/*
 * 简化 SPICE 网表 DC 求解器。
 * 支持：R、V、Q（NPN）和 .model ... NPN IS/BF/BR。
 * 默认 Source Stepping；命令行加 --nr 时使用直接 Newton-Raphson。
 */
#include "dc_solver.h"
#include "dc_timer.h"
#include "gmin_stepping.h"
#include "source_stepping.h"

#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* 网表解析器使用的固定容量与地节点标记。 */
enum { MAX_ELEMENTS = 2048, MAX_MODELS = 128, NAME_SIZE = 64, GROUND = -1 };

/* 电阻、电压源、NPN 模型与 NPN 实例的网表数据。 */
typedef struct { int p, n; double r; } Resistor;
typedef struct {
    int p, n;
    double v;
    double scale;
    char name[NAME_SIZE];
} VoltageSource;
typedef struct {
    int p, n;
    double i;
    double scale;
    char name[NAME_SIZE];
} CurrentSource;
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
    char node_names[DC_MAX_UNKNOWNS][NAME_SIZE];
    int node_count;
    Resistor r[MAX_ELEMENTS]; int r_count;
    VoltageSource v[MAX_ELEMENTS]; int v_count;
    CurrentSource i[MAX_ELEMENTS]; int i_count;
    Bjt q[MAX_ELEMENTS]; int q_count;
    BjtModel model[MAX_MODELS]; int model_count;
    double temperature_celsius;
    double gmin;
} NetlistCircuit;

/* 温度修正后，可直接代入 Ebers-Moll 方程的 BJT 参数。 */
typedef struct { double is, af, ar, vt; } BjtParameters;

/* 子电路定义只保存引脚表与其在源文件中的行范围；展开时再读取对应主体。 */
enum { MAX_SUBCKTS = 64, MAX_SUBCKT_PINS = 32, INSTANCE_PATH_SIZE = 256 };
typedef struct {
    char name[NAME_SIZE];
    char pins[MAX_SUBCKT_PINS][NAME_SIZE];
    int pin_count;
    int first_body_line;
    int end_line;
} SubcktDefinition;

/* 一次 X 实例展开的引脚映射与内部节点命名空间。 */
typedef struct {
    const SubcktDefinition *definition;
    char actual_pins[MAX_SUBCKT_PINS][NAME_SIZE];
    char instance_path[INSTANCE_PATH_SIZE];
} SubcktInstance;

/* 限制指数自变量，避免 Newton 的异常猜测使 exp() 溢出。 */
static double safe_exp(double value)
{
    if (value > 40.0) value = 40.0;
    if (value < -40.0) value = -40.0;
    return exp(value);
}

/* 将网表原始节点号映射为求解器连续使用的内部下标。 */
static bool is_ground_name(const char *name)
{
    return strcmp(name, "0") == 0 || _stricmp(name, "GND") == 0;
}

/* 将数值或文字节点名映射为连续的内部下标，地节点固定为 GROUND。 */
static int node_index(NetlistCircuit *circuit, const char *name)
{
    if (is_ground_name(name)) return GROUND;
    for (int i = 0; i < circuit->node_count; ++i)
        if (_stricmp(circuit->node_names[i], name) == 0) return i;
    if (circuit->node_count >= DC_MAX_UNKNOWNS) return -2;
    snprintf(circuit->node_names[circuit->node_count], NAME_SIZE, "%s", name);
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
    if (tolower((unsigned char)end[0]) == 'v' ||
        tolower((unsigned char)end[0]) == 'a') return value;
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
        snprintf(key, sizeof(key), "%s", token[i]);
        upper_string(key);
        if (strcmp(key, "IS") == 0 || strcmp(key, "BF") == 0 ||
            strcmp(key, "BR") == 0) {
            const double value = parse_value(token[i + 1]);
            if (!isfinite(value)) return false;
            if (strcmp(key, "IS") == 0) model->is = value;
            else if (strcmp(key, "BF") == 0)
                model->af = value / (value + 1.0);
            else model->ar = value / (value + 1.0);
        }
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
static const SubcktDefinition *find_subckt(
    const SubcktDefinition *definitions, int definition_count,
    const char *name)
{
    for (int i = 0; i < definition_count; ++i) {
        if (_stricmp(definitions[i].name, name) == 0) return &definitions[i];
    }
    return NULL;
}

/* 首次扫描只记录每个 .subckt 的名称、引脚表与文件行范围。 */
static bool collect_subcircuits(
    const char *path, SubcktDefinition *definitions, int *definition_count)
{
    FILE *file = fopen(path, "r");
    char line[512];
    int line_number = 0;
    int active = -1;

    if (file == NULL) return false;
    *definition_count = 0;
    while (fgets(line, sizeof(line), file)) {
        char *cursor = line;
        ++line_number;
        while (isspace((unsigned char)*cursor)) ++cursor;
        if (_strnicmp(cursor, ".subckt", 7) == 0 &&
            isspace((unsigned char)cursor[7])) {
            char copy[512];
            char *part;
            SubcktDefinition *definition;
            if (active >= 0 || *definition_count >= MAX_SUBCKTS) {
                fclose(file);
                return false;
            }
            snprintf(copy, sizeof(copy), "%s", cursor);
            part = strtok(copy, " \t\r\n");
            part = strtok(NULL, " \t\r\n");
            if (part == NULL) { fclose(file); return false; }
            definition = &definitions[*definition_count];
            memset(definition, 0, sizeof(*definition));
            snprintf(definition->name, sizeof(definition->name), "%s", part);
            while ((part = strtok(NULL, " \t\r\n")) != NULL) {
                if (definition->pin_count >= MAX_SUBCKT_PINS) {
                    fclose(file);
                    return false;
                }
                snprintf(definition->pins[definition->pin_count++], NAME_SIZE,
                         "%s", part);
            }
            definition->first_body_line = line_number + 1;
            active = (*definition_count)++;
        } else if (_strnicmp(cursor, ".ends", 5) == 0 &&
                   (isspace((unsigned char)cursor[5]) || cursor[5] == '\0')) {
            if (active < 0) { fclose(file); return false; }
            definitions[active].end_line = line_number;
            active = -1;
        }
    }
    fclose(file);
    return active < 0;
}

/* 将子电路引脚名映射到实例实参；内部节点加实例前缀以保证互不连接。 */
static const char *map_subckt_node(
    const SubcktInstance *instance, const char *name,
    char mapped_name[NAME_SIZE])
{
    if (instance == NULL || is_ground_name(name)) return name;
    for (int i = 0; i < instance->definition->pin_count; ++i) {
        if (_stricmp(name, instance->definition->pins[i]) == 0) {
            return instance->actual_pins[i];
        }
    }
    snprintf(mapped_name, NAME_SIZE, "%s/%s", instance->instance_path, name);
    return mapped_name;
}

static void make_instance_name(
    const SubcktInstance *instance, const char *name,
    char result[NAME_SIZE])
{
    if (instance == NULL) snprintf(result, NAME_SIZE, "%s", name);
    else snprintf(result, NAME_SIZE, "%s/%s", instance->instance_path, name);
}

static bool parse_expanded_line(
    const char *path, NetlistCircuit *circuit, char *line,
    const SubcktDefinition *definitions, int definition_count,
    const SubcktInstance *instance, int depth);

/* 展开一个 X 实例：将其引脚实参绑定到定义的形式引脚，再递归处理主体。 */
static bool expand_subckt_instance(
    const char *path, NetlistCircuit *circuit, char *line,
    const SubcktDefinition *definitions, int definition_count,
    const SubcktInstance *parent, int depth)
{
    char copy[512];
    char *token[64];
    int count = 0;
    const SubcktDefinition *definition;
    SubcktInstance instance;
    FILE *file;
    char body_line[512];
    int line_number = 0;

    if (depth >= 16) return false;
    snprintf(copy, sizeof(copy), "%s", line);
    for (char *part = strtok(copy, " \t\r\n"); part != NULL && count < 64;
         part = strtok(NULL, " \t\r\n")) token[count++] = part;
    if (count < 3) return false;
    definition = find_subckt(definitions, definition_count, token[count - 1]);
    if (definition == NULL || count != definition->pin_count + 2) {
        fprintf(stderr, "Unknown subcircuit or wrong pin count in X instance: %s", line);
        return false;
    }

    memset(&instance, 0, sizeof(instance));
    instance.definition = definition;
    if (parent == NULL) {
        snprintf(instance.instance_path, sizeof(instance.instance_path), "%s", token[0]);
    } else {
        snprintf(instance.instance_path, sizeof(instance.instance_path), "%s/%s",
                 parent->instance_path, token[0]);
    }
    for (int i = 0; i < definition->pin_count; ++i) {
        char mapped[NAME_SIZE];
        snprintf(instance.actual_pins[i], NAME_SIZE, "%s",
                 map_subckt_node(parent, token[i + 1], mapped));
    }

    file = fopen(path, "r");
    if (file == NULL) return false;
    while (fgets(body_line, sizeof(body_line), file)) {
        char *cursor = body_line;
        ++line_number;
        if (line_number < definition->first_body_line) continue;
        if (line_number >= definition->end_line) break;
        while (isspace((unsigned char)*cursor)) ++cursor;
        if (*cursor == '\0' || *cursor == '*' || *cursor == ';' ||
            *cursor == '+' || *cursor == '.') continue;
        if (!parse_expanded_line(path, circuit, cursor, definitions,
                                 definition_count, &instance, depth + 1)) {
            fprintf(stderr, "Could not expand subcircuit line: %s", cursor);
            fclose(file);
            return false;
        }
    }
    fclose(file);
    return true;
}

/* 将一条已处于某个实例层级的原始器件行转换为 MNA 使用的内部元件。 */
static bool parse_expanded_line(
    const char *path, NetlistCircuit *circuit, char *line,
    const SubcktDefinition *definitions, int definition_count,
    const SubcktInstance *instance, int depth)
{
    char name[NAME_SIZE], a[NAME_SIZE], b[NAME_SIZE], c[NAME_SIZE];
    char d[NAME_SIZE] = "", e[NAME_SIZE] = "", f[NAME_SIZE] = "";
    char mapped_a[NAME_SIZE], mapped_b[NAME_SIZE], mapped_c[NAME_SIZE];
    char stored_name[NAME_SIZE];
    const int fields = sscanf(line, "%63s %63s %63s %63s %63s %63s %63s",
                              name, a, b, c, d, e, f);
    const char type = fields > 0 ? (char)toupper((unsigned char)name[0]) : '\0';
    const char *node_a = map_subckt_node(instance, a, mapped_a);
    const char *node_b = map_subckt_node(instance, b, mapped_b);
    const char *node_c = map_subckt_node(instance, c, mapped_c);

    if (fields < 1) return true;
    if (type == 'X') {
        return expand_subckt_instance(path, circuit, line, definitions,
                                      definition_count, instance, depth);
    }
    if (type == 'R' && fields >= 4 && circuit->r_count < MAX_ELEMENTS) {
        Resistor *r = &circuit->r[circuit->r_count++];
        r->p = node_index(circuit, node_a); r->n = node_index(circuit, node_b);
        r->r = parse_value(c);
        return r->p != -2 && r->n != -2 && isfinite(r->r) && r->r > 0.0;
    }
    if (type == 'V' && fields >= 4 && circuit->v_count < MAX_ELEMENTS) {
        VoltageSource *v = &circuit->v[circuit->v_count++];
        v->p = node_index(circuit, node_a); v->n = node_index(circuit, node_b);
        v->v = parse_dc_source_value(c, d); v->scale = 1.0;
        if (!isfinite(v->v) &&
            (_stricmp(c, "SIN") == 0 || _stricmp(c, "PULSE") == 0 ||
             _stricmp(c, "PWL") == 0)) {
            v->v = d[0] == '(' ? parse_value(d + 1) : parse_value(e);
        }
        make_instance_name(instance, name, stored_name);
        snprintf(v->name, sizeof(v->name), "%s", stored_name);
        return v->p != -2 && v->n != -2 && isfinite(v->v);
    }
    if (type == 'I' && fields >= 4 && circuit->i_count < MAX_ELEMENTS) {
        CurrentSource *i = &circuit->i[circuit->i_count++];
        i->p = node_index(circuit, node_a); i->n = node_index(circuit, node_b);
        i->i = parse_dc_source_value(c, d); i->scale = 1.0;
        if (!isfinite(i->i) &&
            (_stricmp(c, "SIN") == 0 || _stricmp(c, "PULSE") == 0 ||
             _stricmp(c, "PWL") == 0)) {
            i->i = d[0] == '(' ? parse_value(d + 1) : parse_value(e);
        }
        make_instance_name(instance, name, stored_name);
        snprintf(i->name, sizeof(i->name), "%s", stored_name);
        return i->p != -2 && i->n != -2 && isfinite(i->i);
    }
    if (type == 'Q' && fields >= 5 && circuit->q_count < MAX_ELEMENTS) {
        Bjt *q = &circuit->q[circuit->q_count++];
        const bool has_substrate = fields >= 6 &&
            _strnicmp(e, "area=", 5) != 0;
        const char *model_name = has_substrate ? e : d;
        const char *area_text = has_substrate ? f : e;
        q->c = node_index(circuit, node_a); q->b = node_index(circuit, node_b);
        q->e = node_index(circuit, node_c); q->model = -1; q->area = 1.0;
        snprintf(q->model_name, sizeof(q->model_name), "%s", model_name);
        /* 当前 Ebers-Moll 模型不使用衬底端；其节点仍由层次展开保留。 */
        if (_strnicmp(area_text, "area=", 5) == 0) q->area = parse_value(area_text + 5);
        return q->c != -2 && q->b != -2 && q->e != -2 &&
               isfinite(q->area) && q->area > 0.0;
    }
    if (type == 'C' && fields >= 4) return true;

    fprintf(stderr,
            "Unsupported element '%s' in %s. Current DC backend supports R, I, V, Q and ignores C.\n",
            name, path);
    return false;
}

static bool parse_netlist(const char *path, NetlistCircuit *circuit)
{
    SubcktDefinition definitions[MAX_SUBCKTS];
    int definition_count;
    FILE *file;
    char line[512];
    int line_number = 0;

    if (!collect_subcircuits(path, definitions, &definition_count)) return false;
    file = fopen(path, "r");
    if (file == NULL) return false;
    memset(circuit, 0, sizeof(*circuit));
    circuit->temperature_celsius = 27.0;
    while (fgets(line, sizeof(line), file)) {
        char *cursor = line;
        bool inside_subckt = false;
        ++line_number;
        if (line_number == 1) continue;
        for (int i = 0; i < definition_count; ++i) {
            if (line_number >= definitions[i].first_body_line - 1 &&
                line_number <= definitions[i].end_line) inside_subckt = true;
        }
        if (inside_subckt) continue;
        while (isspace((unsigned char)*cursor)) ++cursor;
        if (*cursor == '\0' || *cursor == '*' || *cursor == ';') continue;
        if (*cursor == '+') {
            char copy[512];
            if (circuit->model_count == 0) { fclose(file); return false; }
            snprintf(copy, sizeof(copy), "%s", cursor);
            if (!parse_model_parameters(&circuit->model[circuit->model_count - 1], copy, 1)) {
                fclose(file); return false;
            }
        } else if (*cursor == '.') {
            char copy[512];
            snprintf(copy, sizeof(copy), "%s", cursor);
            if (_strnicmp(copy, ".model", 6) == 0 && !parse_model(circuit, copy)) {
                fprintf(stderr, "Unsupported or invalid .model line: %s", cursor);
                fclose(file); return false;
            }
            if (_strnicmp(copy, ".include", 8) == 0 || _strnicmp(copy, ".lib", 4) == 0) {
                fprintf(stderr, "Netlist include/library files are not implemented: %s", cursor);
                fclose(file); return false;
            }
        } else if (!parse_expanded_line(path, circuit, cursor, definitions,
                                         definition_count, NULL, 0)) {
            fprintf(stderr, "Could not parse netlist line: %s", cursor);
            fclose(file);
            return false;
        }
    }
    fclose(file);
    if (circuit->model_count == 0) {
        fputs("No supported BJT .model declaration was found.\n", stderr);
        return false;
    }
    if (circuit->node_count + circuit->v_count > DC_MAX_UNKNOWNS) {
        fprintf(stderr, "Expanded circuit dimension %d exceeds current dense limit %d.\n",
                circuit->node_count + circuit->v_count, DC_MAX_UNKNOWNS);
        return false;
    }
    for (int i = 0; i < circuit->q_count; ++i) {
        circuit->q[i].model = find_model(circuit, circuit->q[i].model_name);
        if (circuit->q[i].model < 0) {
            fprintf(stderr, "BJT model '%s' is not defined.\n",
                    circuit->q[i].model_name);
            return false;
        }
    }
    if (circuit->r_count == 0 || circuit->v_count + circuit->i_count == 0 ||
        circuit->q_count == 0) {
        fputs("The expanded circuit lacks required R, independent-source, or Q elements.\n",
              stderr);
        return false;
    }
    return true;
}

/* Reorder voltage sources according to a comma-separated source-name list. */
static bool reorder_voltage_sources(NetlistCircuit *circuit, const char *text)
{
    char buffer[512];
    VoltageSource ordered[MAX_ELEMENTS];
    bool used[MAX_ELEMENTS] = { false };
    int count = 0;

    snprintf(buffer, sizeof(buffer), "%s", text);
    for (char *item = strtok(buffer, ","); item != NULL;
         item = strtok(NULL, ",")) {
        char *end;
        int found = -1;

        while (isspace((unsigned char)*item)) ++item;
        end = item + strlen(item);
        while (end > item && isspace((unsigned char)end[-1])) --end;
        *end = '\0';

        for (int i = 0; i < circuit->v_count; ++i) {
            if (!used[i] && _stricmp(item, circuit->v[i].name) == 0) {
                found = i;
                break;
            }
        }
        if (found < 0 || count >= circuit->v_count) return false;
        used[found] = true;
        ordered[count++] = circuit->v[found];
    }

    if (count != circuit->v_count) return false;
    memcpy(circuit->v, ordered, sizeof(VoltageSource) * (size_t)count);
    return true;
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
    /* GMIN Stepping：为每个非地节点写入一个到地的附加电导。 */
    if (circuit->gmin > 0.0) {
        for (int node = 0; node < n; ++node) {
            f[node] += circuit->gmin * x[node];
            j[node][node] += circuit->gmin;
        }
    }
    for (int i = 0; i < circuit->v_count; ++i) {
        const VoltageSource *v = &circuit->v[i]; const int branch = n + i;
        if (v->p != GROUND) { f[v->p] += x[branch]; j[v->p][branch] += 1.0; j[branch][v->p] += 1.0; }
        if (v->n != GROUND) { f[v->n] -= x[branch]; j[v->n][branch] -= 1.0; j[branch][v->n] -= 1.0; }
        f[branch] += node_voltage(x, v->p) - node_voltage(x, v->n) -
                     lambda * v->scale * v->v;
    }
    for (int i = 0; i < circuit->i_count; ++i) {
        const CurrentSource *source = &circuit->i[i];
        const double current = lambda * source->scale * source->i;
        if (source->p != GROUND) f[source->p] += current;
        if (source->n != GROUND) f[source->n] -= current;
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
  static double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS]; build_system(c,x,l,f,j); }
/* 为通用 DC 求解器构造 Jacobian 矩阵。 */
static void jacobian(const void *c, const double *x, double l, double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{
  double f[DC_MAX_UNKNOWNS]; build_system(c,x,l,f,j); }

/* 为顺序 Source Stepping 设置网表中一个电压源的缩放比例。 */
static void set_netlist_source_scale(
    const void *context,
    int source_index,
    double scale)
{
    NetlistCircuit *circuit = (NetlistCircuit *)context;
    if (source_index >= 0 && source_index < circuit->v_count) {
        circuit->v[source_index].scale = scale;
    } else if (source_index >= circuit->v_count &&
               source_index < circuit->v_count + circuit->i_count) {
        circuit->i[source_index - circuit->v_count].scale = scale;
    }
}

/* 为 GMIN Stepping 设置当前所有非地节点对地的附加电导。 */
static void set_netlist_gmin(const void *context, double gmin)
{
    NetlistCircuit *circuit = (NetlistCircuit *)context;
    circuit->gmin = gmin;
}

/* 返回顺序 Source Stepping 当前独立源的名称，电压源排在电流源之前。 */
static const char *independent_source_name(
    const NetlistCircuit *circuit,
    int source_index)
{
    if (source_index >= 0 && source_index < circuit->v_count) {
        return circuit->v[source_index].name;
    }
    if (source_index >= circuit->v_count &&
        source_index < circuit->v_count + circuit->i_count) {
        return circuit->i[source_index - circuit->v_count].name;
    }
    return "unknown_source";
}


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

/* 保存一次实验的逐步数据与汇总数据。 */
typedef struct {
    FILE *step_file;
    FILE *summary_file;
    const NetlistCircuit *circuit;
    const char *algorithm;
    double previous_time_ms;
    int accepted_steps;
    int accepted_newton_iterations;
    int accepted_line_search_reductions;
} CsvReporter;

/* 建立 CSV 报告文件；输出目录需由调用者预先创建。 */
static bool csv_report_open(CsvReporter *reporter, const char *directory,
                            const NetlistCircuit *circuit,
                            const char *algorithm)
{
    char step_path[512];
    char summary_path[512];

    memset(reporter, 0, sizeof(*reporter));
    reporter->circuit = circuit;
    reporter->algorithm = algorithm;
    reporter->previous_time_ms = dc_timer_now_milliseconds();
    if (directory == NULL) return true;

    snprintf(step_path, sizeof(step_path), "%s/dc_steps.csv", directory);
    snprintf(summary_path, sizeof(summary_path), "%s/dc_summary.csv", directory);
    reporter->step_file = fopen(step_path, "w");
    reporter->summary_file = fopen(summary_path, "w");
    if (reporter->step_file == NULL || reporter->summary_file == NULL) {
        if (reporter->step_file != NULL) fclose(reporter->step_file);
        if (reporter->summary_file != NULL) fclose(reporter->summary_file);
        memset(reporter, 0, sizeof(*reporter));
        return false;
    }

    fprintf(reporter->step_file,
            "algorithm,step,parameter_kind,parameter_value,source,"
            "newton_iterations,line_search_reductions,residual_inf,"
            "elapsed_since_previous_ms");
    for (int i = 0; i < circuit->node_count + circuit->v_count; ++i) {
        fprintf(reporter->step_file, ",x%d", i + 1);
    }
    fputc('\n', reporter->step_file);
    fprintf(reporter->summary_file,
            "algorithm,converged,total_newton_iterations,accepted_steps,"
            "accepted_newton_iterations,failed_probe_newton_iterations,"
            "accepted_line_search_reductions,total_runtime_ms,"
            "final_residual_inf,dimension,node_count,voltage_source_count,current_source_count,"
            "temperature_celsius,initial_lambda_step,maximum_lambda_step,"
            "initial_gmin,minimum_gmin,gmin_reduction_factor\n");
    return true;
}

/* 将一个成功求解点写入逐步 CSV；每行同时保存该点的未知量。 */
static void csv_report_step(CsvReporter *reporter, const char *parameter_kind,
                            double parameter_value, const char *source,
                            const double *x,
                            const DcNewtonReport *newton_report)
{
    const double now = dc_timer_now_milliseconds();
    const double elapsed_ms = now - reporter->previous_time_ms;
    reporter->previous_time_ms = now;
    ++reporter->accepted_steps;
    reporter->accepted_newton_iterations += newton_report->iterations;
    reporter->accepted_line_search_reductions +=
        newton_report->line_search_reductions;
    if (reporter->step_file == NULL) return;

    fprintf(reporter->step_file, "%s,%d,%s,%.17g,%s,%d,%d,%.17g,%.3f",
            reporter->algorithm, reporter->accepted_steps, parameter_kind,
            parameter_value, source == NULL ? "" : source,
            newton_report->iterations, newton_report->line_search_reductions,
            newton_report->final_residual_norm, elapsed_ms);
    for (int i = 0; i < reporter->circuit->node_count +
                        reporter->circuit->v_count; ++i) {
        fprintf(reporter->step_file, ",%.17g", x[i]);
    }
    fputc('\n', reporter->step_file);
}

/* 写入本次运行的单行汇总数据并关闭两个 CSV 文件。 */
static void csv_report_close(CsvReporter *reporter, bool converged,
                             int total_newton_iterations,
                             double total_runtime_ms,
                             double final_residual_norm,
                             const DcSolverOptions *options)
{
    if (reporter->summary_file != NULL) {
        const NetlistCircuit *circuit = reporter->circuit;
        fprintf(reporter->summary_file,
                "%s,%d,%d,%d,%d,%d,%d,%.3f,%.17g,%d,%d,%d,%d,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%.17g\n",
                reporter->algorithm, converged ? 1 : 0,
                total_newton_iterations, reporter->accepted_steps,
                reporter->accepted_newton_iterations,
                total_newton_iterations - reporter->accepted_newton_iterations,
                reporter->accepted_line_search_reductions, total_runtime_ms,
                final_residual_norm, circuit->node_count + circuit->v_count,
                circuit->node_count, circuit->v_count, circuit->i_count,
                circuit->temperature_celsius, options->initial_lambda_step,
                options->maximum_lambda_step, options->initial_gmin,
                options->minimum_gmin, options->gmin_reduction_factor);
    }
    if (reporter->step_file != NULL) fclose(reporter->step_file);
    if (reporter->summary_file != NULL) fclose(reporter->summary_file);
}

/* 返回残差向量中绝对值最大的分量。 */
static double infinity_norm(const double *values, int count)
{
    double maximum = 0.0;
    for (int i = 0; i < count; ++i) {
        const double value = fabs(values[i]);
        if (value > maximum) maximum = value;
    }
    return maximum;
}

/* 输出并记录每个成功同时 Source Stepping 点。 */
static void print_step(double lambda, const double *x,
                       const DcNewtonReport *report, void *context)
{
    CsvReporter *reporter = context;
    printf("lambda = %.4f, Newton = %d, residual = %.3e\n", lambda,
           report->iterations, report->final_residual_norm);
    csv_report_step(reporter, "lambda", lambda, "", x, report);
}

/* 输出并记录顺序 Source Stepping 中当前电源和局部 lambda。 */
static void print_sequential_step(
    int source_index, double lambda, const double *x,
    const DcNewtonReport *report, void *context)
{
    CsvReporter *reporter = context;
    const NetlistCircuit *circuit = reporter->circuit;
    const char *source_name = source_index < 0 ? "all_sources_off" :
                              independent_source_name(circuit, source_index);
    if (source_index < 0) {
        printf("all sources off, Newton = %d, residual = %.3e\n",
               report->iterations, report->final_residual_norm);
    } else {
        printf("source %s, lambda = %.4f, Newton = %d, residual = %.3e\n",
               source_name, lambda, report->iterations,
               report->final_residual_norm);
    }
    csv_report_step(reporter, "source_lambda", lambda, source_name, x, report);
}

/* 输出并记录每个成功 GMIN 步进点。 */
static void print_gmin_step(double gmin, const double *x,
                            const DcNewtonReport *report, void *context)
{
    CsvReporter *reporter = context;
    printf("gmin = %.3e S, Newton = %d, residual = %.3e\n", gmin,
           report->iterations, report->final_residual_norm);
    csv_report_step(reporter, "gmin_siemens", gmin, "", x, report);
}

/* 读取网表与命令行选项，执行 DC 求解并输出最终工作点。 */
int main(int argc, char **argv)
{
    const char *path;
    bool direct_nr = false;
    bool sequential_sources = false;
    bool gmin_stepping = false;
    const char *source_order = NULL;
    const char *report_directory = NULL;
    NetlistCircuit circuit;
    DcProblem problem;
    DcSolverOptions options = dc_solver_default_options();
    double x[DC_MAX_UNKNOWNS] = {0.0};
    double requested_temperature = NAN;
    DcNewtonReport newton_report;
    CsvReporter reporter;
    const char *algorithm;
    int total_newton_iterations = 0;
    double start_time_ms;

    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s netlist.txt [--nr] [--sequential-sources] [--gmin-stepping] [--source-order names] [--secant-predictor] "
                "[--temp Celsius] [--initial-step value] [--max-step value] "
                "[--fast-threshold count] [--slow-threshold count] "
                "[--growth-factor value] [--shrink-factor value] "
                "[--junction-step volts] [--gmin-initial siemens] "
                "[--gmin-min siemens] [--gmin-factor value] "
                "[--report-dir existing-directory]\n",
                argv[0]);
        return 1;
    }

    path = argv[1];
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--nr") == 0) {
            direct_nr = true;
        } else if (strcmp(argv[i], "--sequential-sources") == 0) {
            sequential_sources = true;
        } else if (strcmp(argv[i], "--gmin-stepping") == 0) {
            gmin_stepping = true;
        } else if (strcmp(argv[i], "--source-order") == 0 && i + 1 < argc) {
            source_order = argv[++i];
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
        } else if (strcmp(argv[i], "--gmin-initial") == 0 && i + 1 < argc) {
            options.initial_gmin = atof(argv[++i]);
        } else if (strcmp(argv[i], "--gmin-min") == 0 && i + 1 < argc) {
            options.minimum_gmin = atof(argv[++i]);
        } else if (strcmp(argv[i], "--gmin-factor") == 0 && i + 1 < argc) {
            options.gmin_reduction_factor = atof(argv[++i]);
        } else if (strcmp(argv[i], "--report-dir") == 0 && i + 1 < argc) {
            report_directory = argv[++i];
        } else {
            fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            return 1;
        }
    }

    if (!parse_netlist(path, &circuit)) {
        fprintf(stderr, "Could not parse supported R/V/Q NPN/PNP netlist: %s\n", path);
        return 1;
    }
    if (source_order != NULL &&
        !reorder_voltage_sources(&circuit, source_order)) {
        fputs("--source-order must list every voltage-source name once.\n", stderr);
        return 1;
    }
    if (isfinite(requested_temperature)) {
        circuit.temperature_celsius = requested_temperature;
    }
    if ((direct_nr ? 1 : 0) + (sequential_sources ? 1 : 0) +
        (gmin_stepping ? 1 : 0) > 1) {
        fputs("--nr, --sequential-sources and --gmin-stepping are mutually exclusive.\n", stderr);
        return 1;
    }
    problem = (DcProblem) {
        .dimension = circuit.node_count + circuit.v_count,
        .context = &circuit,
        .build_residual = residual,
        .build_jacobian = jacobian,
        .limit_newton_step = limit_bjt_junction_steps,
        .independent_source_count = circuit.v_count + circuit.i_count,
        .set_source_scale = set_netlist_source_scale,
        .set_gmin = set_netlist_gmin
    };
    algorithm = direct_nr ? "direct_newton" :
                sequential_sources ? "sequential_source_stepping" :
                gmin_stepping ? "gmin_stepping" : "source_stepping";
    if (!csv_report_open(&reporter, report_directory, &circuit, algorithm)) {
        fputs("Could not create CSV report files; create the report directory first.\n",
              stderr);
        return 1;
    }
    start_time_ms = dc_timer_now_milliseconds();
    if (direct_nr) {
        if (!dc_newton_solve_with_report(&problem, &options, 1.0, x,
                                         &newton_report)) {
            fputs("Newton did not converge.\n", stderr);
            csv_report_close(&reporter, false, 0, 0.0, NAN, &options);
            return 1;
        }
        printf("Direct Newton converged in %d iterations, %d line-search reductions.\n",
               newton_report.iterations,
               newton_report.line_search_reductions);
        total_newton_iterations = newton_report.iterations;
        csv_report_step(&reporter, "none", 1.0, "", x, &newton_report);
    } else if (sequential_sources) {
        if (!dc_sequential_source_stepping_solve(
                &problem, &options, x, print_sequential_step, &reporter,
                &total_newton_iterations)) {
            fputs("Sequential source stepping failed.\n", stderr);
            csv_report_close(&reporter, false, total_newton_iterations,
                             dc_timer_now_milliseconds() - start_time_ms,
                             NAN, &options);
            return 1;
        }
        printf("Total Newton iterations: %d\n", total_newton_iterations);
    } else if (gmin_stepping) {
        if (!dc_gmin_stepping_solve(&problem, &options, x,
                                    print_gmin_step, &reporter,
                                    &total_newton_iterations)) {
            fputs("GMIN stepping failed.\n", stderr);
            csv_report_close(&reporter, false, total_newton_iterations,
                             dc_timer_now_milliseconds() - start_time_ms,
                             NAN, &options);
            return 1;
        }
        printf("Total Newton iterations: %d\n", total_newton_iterations);
    } else {
        if (!dc_source_stepping_solve(&problem, &options, x, print_step,
                                      &reporter, &total_newton_iterations)) {
            fputs("Source stepping failed.\n", stderr);
            csv_report_close(&reporter, false, total_newton_iterations,
                             dc_timer_now_milliseconds() - start_time_ms,
                             NAN, &options);
            return 1;
        }
        printf("Total Newton iterations: %d\n", total_newton_iterations);
    }
    {
        double final_residual[DC_MAX_UNKNOWNS];
        const double runtime_ms = dc_timer_now_milliseconds() - start_time_ms;
        residual(&circuit, x, 1.0, final_residual);
        csv_report_close(&reporter, true, total_newton_iterations, runtime_ms,
                         infinity_norm(final_residual, problem.dimension),
                         &options);
        printf("Runtime: %.3f ms\n", runtime_ms);
    }
    for (int i = 0; i < circuit.node_count; ++i)
        printf("V(%s) = % .9f V\n", circuit.node_names[i], x[i]);
    for (int i = 0; i < circuit.v_count; ++i)
        printf("I(%s) = % .9e A\n", circuit.v[i].name, x[circuit.node_count + i]);
    return 0;
}
