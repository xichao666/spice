/*
 * 简化 SPICE 网表 DC 求解器。
 * 支持：R、V、Q（NPN）和 .model ... NPN IS/BF/BR。
 * 默认 Source Stepping；命令行加 --nr 时使用直接 Newton-Raphson。
 */
#include "dc_solver.h"
#include "source_stepping.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 网表解析器使用的固定容量与地节点标记。 */
enum { MAX_ELEMENTS = 128, MAX_MODELS = 16, NAME_SIZE = 32, GROUND = -1 };

/* 电阻、电压源、NPN 模型与 NPN 实例的网表数据。 */
typedef struct { int p, n; double r; } Resistor;
typedef struct { int p, n; double v; char name[NAME_SIZE]; } VoltageSource;
typedef struct {
    char name[NAME_SIZE];
    double is, af, ar;
    int polarity; /* NPN 为 +1，PNP 为 -1。*/
} BjtModel;
typedef struct {
    int c, b, e, model;
    double area;
    char model_name[NAME_SIZE];
} Bjt;
/* 已解析的电路数据与温度设置。 */
typedef struct {
    int node_ids[DC_MAX_UNKNOWNS], node_count;
    Resistor r[MAX_ELEMENTS]; int r_count;
    VoltageSource v[MAX_ELEMENTS]; int v_count;
    Bjt q[MAX_ELEMENTS]; int q_count;
    BjtModel model[MAX_MODELS]; int model_count;
    double temperature_celsius;
} NetlistCircuit;

/* 温度修正后，可直接代入 Ebers-Moll 方程的 BJT 参数。 */
typedef struct { double is, af, ar, vt; } BjtParameters;

/* 限制指数自变量，避免 Newton 的异常猜测使 exp() 溢出。 */
static double safe_exp(double value)
{
    if (value > 40.0) value = 40.0;
    if (value < -40.0) value = -40.0;
    return exp(value);
}

/* 将网表原始节点号映射为求解器连续使用的内部下标。 */
static int node_index(NetlistCircuit *circuit, int node_id)
{
    if (node_id == 0) return GROUND;
    for (int i = 0; i < circuit->node_count; ++i)
        if (circuit->node_ids[i] == node_id) return i;
    if (circuit->node_count >= DC_MAX_UNKNOWNS) return -2;
    circuit->node_ids[circuit->node_count] = node_id;
    return circuit->node_count++;
}

/* 返回指定节点的电压；地节点电压固定为 0 V。 */
static double node_voltage(const double *x, int node)
{
    return node == GROUND ? 0.0 : x[node];
}

/* 将一个数值累加到 Jacobian 指定位置，地节点不参与矩阵。 */
static void add_j(double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
                  int row, int column, double value)
{
    if (row != GROUND && column != GROUND) j[row][column] += value;
}

/* 读取网表数值与常用 SPICE 后缀。 */
static double parse_value(const char *text)
{
    char *end;
    double value = strtod(text, &end);
    if (*end == '\0') return value;
    if (tolower((unsigned char)end[0]) == 'k') return value * 1e3;
    if (tolower((unsigned char)end[0]) == 'm' &&
        tolower((unsigned char)end[1]) == 'e' &&
        tolower((unsigned char)end[2]) == 'g') return value * 1e6;

    if (end[0] == 'm') return value * 1e-3;
    if (tolower((unsigned char)end[0]) == 'u') return value * 1e-6;
    if (tolower((unsigned char)end[0]) == 'n') return value * 1e-9;
    if (tolower((unsigned char)end[0]) == 'p') return value * 1e-12;
    return NAN;
}

/* 读取独立电压源在 DC 工作点时采用的电压值。*/
static double parse_dc_source_value(const char *value, const char *next_value)
{
    if (_stricmp(value, "DC") == 0) return parse_value(next_value);
    if (_strnicmp(value, "SIN(", 4) == 0) return parse_value(value + 4);
    if (_strnicmp(value, "PULSE(", 6) == 0) return parse_value(value + 6);
    if (_strnicmp(value, "PWL(", 4) == 0) return parse_value(next_value);
    return parse_value(value);
}

/* 将字符串转为大写，供参数名比较使用。 */
static void upper_string(char *text)
{
    for (; *text; ++text) *text = (char)toupper((unsigned char)*text);
}

/* 在已读取的模型中按名称查找指定模型。 */
static int find_model(const NetlistCircuit *circuit, const char *name)
{
    for (int i = 0; i < circuit->model_count; ++i)
        if (_stricmp(circuit->model[i].name, name) == 0) return i;
    return -1;
}

/* 将一段模型参数写入当前 BJT 模型。*/
static bool parse_model_parameters(BjtModel *model, char *line, int first_parameter)
{
    char *token[32];
    int count = 0;
    for (char *part = strtok(line, " \t=()\r\n"); part && count < 32;
         part = strtok(NULL, " \t=()\r\n")) {
        token[count++] = part;
    }
    for (int i = first_parameter; i + 1 < count; i += 2) {
        char key[NAME_SIZE];
        const double value = parse_value(token[i + 1]);
        snprintf(key, sizeof(key), "%s", token[i]);
        upper_string(key);
        if (!isfinite(value)) return false;
        if (strcmp(key, "IS") == 0) model->is = value;
        else if (strcmp(key, "BF") == 0) model->af = value / (value + 1.0);
        else if (strcmp(key, "BR") == 0) model->ar = value / (value + 1.0);
    }
    return true;
}

/* 读取一行 .model NPN 定义，并保存 IS、BF、BR 参数。 */
static bool parse_model(NetlistCircuit *circuit, char *line)
{
    char original[512];
    char *token[32]; int count = 0;
    snprintf(original, sizeof(original), "%s", line);
    for (char *part = strtok(line, " \t=()\r\n"); part && count < 32;
         part = strtok(NULL, " \t=()\r\n")) token[count++] = part;
    if (count < 3 || _stricmp(token[0], ".model") != 0 ||
        (_stricmp(token[2], "NPN") != 0 && _stricmp(token[2], "PNP") != 0) ||
        circuit->model_count >= MAX_MODELS)
        return false;
    BjtModel *model = &circuit->model[circuit->model_count++];
    snprintf(model->name, sizeof(model->name), "%s", token[1]);
    model->is = 1e-14; model->af = 100.0 / 101.0; model->ar = 0.5;
    model->polarity = _stricmp(token[2], "NPN") == 0 ? 1 : -1;
    return parse_model_parameters(model, original, 3);
}

/* 读取整个网表，识别 R、V、Q 与 .model 元件行。 */
static bool parse_netlist(const char *path, NetlistCircuit *circuit)
{
    FILE *file = fopen(path, "r");
    char line[512];
    if (!file) return false;
    memset(circuit, 0, sizeof(*circuit));
    circuit->temperature_celsius = 27.0;
    while (fgets(line, sizeof(line), file)) {
        char *cursor = line;
        while (isspace((unsigned char)*cursor)) ++cursor;
        if (*cursor == '\0' || *cursor == '*' || *cursor == ';') continue;
        if (*cursor == '+') {
            char copy[512];
            if (circuit->model_count == 0) { fclose(file); return false; }
            snprintf(copy, sizeof(copy), "%s", cursor);
            if (!parse_model_parameters(&circuit->model[circuit->model_count - 1], copy, 1)) {
                fclose(file); return false;
            }
            continue;
        }
        if (*cursor == '.') {
            char copy[512]; snprintf(copy, sizeof(copy), "%s", cursor);
            if (_strnicmp(copy, ".model", 6) == 0 && !parse_model(circuit, copy)) {
                fclose(file); return false;
            }
            continue;
        }
        char name[NAME_SIZE], a[NAME_SIZE], b[NAME_SIZE], c[NAME_SIZE];
        char d[NAME_SIZE] = "", e[NAME_SIZE] = "", f[NAME_SIZE] = "";
        const int fields = sscanf(cursor, "%31s %31s %31s %31s %31s %31s %31s",
                                  name, a, b, c, d, e, f);
        if (fields < 1) continue;
        const char type = (char)toupper((unsigned char)name[0]);
        if (type == 'R' && fields == 4 && circuit->r_count < MAX_ELEMENTS) {
            Resistor *r = &circuit->r[circuit->r_count++];
            r->p = node_index(circuit, atoi(a)); r->n = node_index(circuit, atoi(b)); r->r = parse_value(c);
            if (r->p == -2 || r->n == -2 || !isfinite(r->r) || r->r <= 0.0) { fclose(file); return false; }
        } else if (type == 'V' && fields >= 4 && circuit->v_count < MAX_ELEMENTS) {
            VoltageSource *v = &circuit->v[circuit->v_count++];
            v->p = node_index(circuit, atoi(a)); v->n = node_index(circuit, atoi(b));
            v->v = parse_dc_source_value(c, d);
            snprintf(v->name, sizeof(v->name), "%s", name);
            if (v->p == -2 || v->n == -2 || !isfinite(v->v)) { fclose(file); return false; }
        } else if (type == 'Q' && fields >= 5 && circuit->q_count < MAX_ELEMENTS) {
            Bjt *q = &circuit->q[circuit->q_count++];
            const bool has_substrate = strspn(d, "+-0123456789") == strlen(d);
            const char *model_name = has_substrate ? e : d;
            const char *area_text = has_substrate ? f : e;
            q->c = node_index(circuit, atoi(a)); q->b = node_index(circuit, atoi(b)); q->e = node_index(circuit, atoi(c)); q->model = -1;
            q->area = 1.0;
            snprintf(q->model_name, sizeof(q->model_name), "%s", model_name);
            if (has_substrate && atoi(d) != 0) { fclose(file); return false; }
            if (_strnicmp(area_text, "area=", 5) == 0)
                q->area = parse_value(area_text + 5);
            if (q->c == -2 || q->b == -2 || q->e == -2) { fclose(file); return false; }
            if (!isfinite(q->area) || q->area <= 0.0) { fclose(file); return false; }
            (void)d;
        }
    }
    fclose(file);
    if (circuit->model_count == 0 || circuit->node_count + circuit->v_count > DC_MAX_UNKNOWNS) return false;
    for (int i = 0; i < circuit->q_count; ++i) {
        circuit->q[i].model = find_model(circuit, circuit->q[i].model_name);
        if (circuit->q[i].model < 0) return false;
    }
    return circuit->r_count > 0 && circuit->v_count > 0 && circuit->q_count > 0;
}

/* 根据温度和模型参数计算本次求解的 BJT 参数。 */
static BjtParameters bjt_parameters(const NetlistCircuit *circuit, int model_index)
{
    const BjtModel *m = &circuit->model[model_index];
    const double t = circuit->temperature_celsius + 273.15;
    const double tnom = 300.15, kq = 8.617333262e-5;
    const double is = m->is * pow(t / tnom, 3.0) *
        exp(-1.11 / kq * (1.0 / t - 1.0 / tnom));
    return (BjtParameters) { is, m->af, m->ar, kq * t };
}

/* 将网表中全部元件的贡献写入 MNA 残差与 Jacobian。 */
static void build_system(const void *context, const double *x, double lambda,
                         double *f, double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
    const NetlistCircuit *circuit = context;
    const int n = circuit->node_count, dimension = n + circuit->v_count;
    memset(f, 0, sizeof(double) * (size_t)dimension);
    memset(j, 0, sizeof(double) * DC_MAX_UNKNOWNS * DC_MAX_UNKNOWNS);
    for (int i = 0; i < circuit->r_count; ++i) {
        const Resistor *r = &circuit->r[i]; const double g = 1.0 / r->r;
        const double current = g * (node_voltage(x, r->p) - node_voltage(x, r->n));
        if (r->p != GROUND) f[r->p] += current;
        if (r->n != GROUND) f[r->n] -= current;
        add_j(j, r->p, r->p, g); add_j(j, r->p, r->n, -g);
        add_j(j, r->n, r->p, -g); add_j(j, r->n, r->n, g);
    }
    for (int i = 0; i < circuit->v_count; ++i) {
        const VoltageSource *v = &circuit->v[i]; const int branch = n + i;
        if (v->p != GROUND) { f[v->p] += x[branch]; j[v->p][branch] += 1.0; j[branch][v->p] += 1.0; }
        if (v->n != GROUND) { f[v->n] -= x[branch]; j[v->n][branch] -= 1.0; j[branch][v->n] -= 1.0; }
        f[branch] += node_voltage(x, v->p) - node_voltage(x, v->n) - lambda * v->v;
    }
    for (int i = 0; i < circuit->q_count; ++i) {
        const Bjt *q = &circuit->q[i]; const BjtParameters p = bjt_parameters(circuit, q->model);
        const int polarity = circuit->model[q->model].polarity;
        const double vc = node_voltage(x, q->c), vb = node_voltage(x, q->b), ve = node_voltage(x, q->e);
        const double ef = safe_exp(polarity * (vb - ve) / p.vt);
        const double er = safe_exp(polarity * (vb - vc) / p.vt);
        const double is = p.is * q->area;
        const double ifc = is * (ef - 1.0), irc = is * (er - 1.0), gf = is * ef / p.vt, gr = is * er / p.vt;
        const double ic = polarity * (ifc - irc / p.ar);
        const double ib = polarity * ((1.0 / p.af - 1.0) * ifc + (1.0 / p.ar - 1.0) * irc);
        const double ie = polarity * (ifc / p.af - irc);
        if (q->c != GROUND) f[q->c] += ic;
        if (q->b != GROUND) f[q->b] += ib;
        if (q->e != GROUND) f[q->e] -= ie;
        add_j(j,q->c,q->c,gr/p.ar); add_j(j,q->c,q->b,gf-gr/p.ar); add_j(j,q->c,q->e,-gf);
        add_j(j,q->b,q->c,-(1.0/p.ar-1.0)*gr); add_j(j,q->b,q->b,(1.0/p.af-1.0)*gf+(1.0/p.ar-1.0)*gr); add_j(j,q->b,q->e,-(1.0/p.af-1.0)*gf);
        add_j(j,q->e,q->c,-gr); add_j(j,q->e,q->b,-gf/p.af+gr); add_j(j,q->e,q->e,gf/p.af);
    }
}

/* 为通用 DC 求解器构造残差向量。 */
static void residual(const void *c, const double *x, double l, double *f)
{
  double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS]; build_system(c,x,l,f,j); }
/* 为通用 DC 求解器构造 Jacobian 矩阵。 */
static void jacobian(const void *c, const double *x, double l, double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
  double f[DC_MAX_UNKNOWNS]; build_system(c,x,l,f,j); }

/*
 * Limit the change of every BJT base-emitter and base-collector voltage.
 * The whole Newton correction is scaled so that its direction is preserved.
 */
static void limit_bjt_junction_steps(
    const void *context,
    const double *x,
    double *delta,
    double maximum_junction_voltage_step)
{
    const NetlistCircuit *circuit = context;
    double scale = 1.0;

    (void)x;
    if (maximum_junction_voltage_step <= 0.0) return;

    for (int i = 0; i < circuit->q_count; ++i) {
        const Bjt *q = &circuit->q[i];
        const int polarity = circuit->model[q->model].polarity;
        const double dc = q->c == GROUND ? 0.0 : delta[q->c];
        const double db = q->b == GROUND ? 0.0 : delta[q->b];
        const double de = q->e == GROUND ? 0.0 : delta[q->e];
        const double d_vbe = polarity * (db - de);
        const double d_vbc = polarity * (db - dc);
        const double junction_steps[2] = { d_vbe, d_vbc };

        for (int k = 0; k < 2; ++k) {
            const double magnitude = fabs(junction_steps[k]);
            if (magnitude > maximum_junction_voltage_step) {
                const double candidate =
                    maximum_junction_voltage_step / magnitude;
                if (candidate < scale) scale = candidate;
            }
        }
    }

    if (scale < 1.0) {
        const int dimension = circuit->node_count + circuit->v_count;
        for (int i = 0; i < dimension; ++i) delta[i] *= scale;
    }
}

/* 输出每个成功 Source Stepping 点的参数和 Newton 次数。 */
static void print_step(double lambda, const double *x, int iterations, void *context)
{
  (void)x; (void)context; printf("lambda = %.4f, Newton = %d\n", lambda, iterations); }

/* 读取网表与命令行选项，执行 DC 求解并输出最终工作点。 */
int main(int argc, char **argv)
{
    const char *path;
    bool direct_nr = false;
    NetlistCircuit circuit;
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS] = {0.0};
    double requested_temperature = NAN;
    DcNewtonReport newton_report;

    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s netlist.txt [--nr] [--secant-predictor] "
                "[--temp Celsius] [--initial-step value] [--max-step value] "
                "[--fast-threshold count] [--slow-threshold count] "
                "[--growth-factor value] [--shrink-factor value] "
                "[--junction-step volts]\n",
                argv[0]);
        return 1;
    }

    path = argv[1];
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--nr") == 0) {
            direct_nr = true;
        } else if (strcmp(argv[i], "--secant-predictor") == 0) {
            options.source_step_policy = DC_SOURCE_STEP_SECANT_PREDICTOR;
        } else if (strcmp(argv[i], "--temp") == 0 && i + 1 < argc) {
            requested_temperature = atof(argv[++i]);
        } else if (strcmp(argv[i], "--initial-step") == 0 && i + 1 < argc) {
            options.initial_lambda_step = atof(argv[++i]);
        } else if (strcmp(argv[i], "--max-step") == 0 && i + 1 < argc) {
            options.maximum_lambda_step = atof(argv[++i]);
        } else if (strcmp(argv[i], "--fast-threshold") == 0 && i + 1 < argc) {
            options.fast_newton_iteration_threshold = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--slow-threshold") == 0 && i + 1 < argc) {
            options.slow_newton_iteration_threshold = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--growth-factor") == 0 && i + 1 < argc) {
            options.lambda_step_growth_factor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--shrink-factor") == 0 && i + 1 < argc) {
            options.lambda_step_shrink_factor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--junction-step") == 0 && i + 1 < argc) {
            options.maximum_junction_voltage_step = atof(argv[++i]);
        } else {
            fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            return 1;
        }
    }

    if (!parse_netlist(path, &circuit)) {
        fprintf(stderr, "Could not parse supported R/V/Q NPN/PNP netlist: %s\n", path);
        return 1;
    }
    if (isfinite(requested_temperature)) {
        circuit.temperature_celsius = requested_temperature;
    }
    problem = (DcProblem) {
        .dimension = circuit.node_count + circuit.v_count,
        .context = &circuit,
        .build_residual = residual,
        .build_jacobian = jacobian,
        .limit_newton_step = limit_bjt_junction_steps
    };
    if (direct_nr) {
        if (!dc_newton_solve_with_report(&problem, &options, 1.0, x,
                                         &newton_report)) {
            fputs("Newton did not converge.\n", stderr);
            return 1;
        }
        printf("Direct Newton converged in %d iterations, %d line-search reductions.\n",
               newton_report.iterations,
               newton_report.line_search_reductions);
    } else {
        int total = 0;
        if (!dc_source_stepping_solve(&problem,&options,x,print_step,NULL,&total)) { fputs("Source stepping failed.\n",stderr); return 1; }
        printf("Total Newton iterations: %d\n", total);
    }
    for (int i = 0; i < circuit.node_count; ++i) printf("V(%d) = % .9f V\n", circuit.node_ids[i], x[i]);
    for (int i = 0; i < circuit.v_count; ++i)
        printf("I(%s) = % .9e A\n", circuit.v[i].name, x[circuit.node_count + i]);
    return 0;
}
