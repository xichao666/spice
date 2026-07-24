/* 根页面三个 modified-nodal 基准的独立实验；不依赖 examples/ 的简化方程。 */
#include "homotopy.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static double safe_exp(double value)
{
    if (value > 40.0) value = 40.0;
    if (value < -40.0) value = -40.0;
    return exp(value);
}

/* 按参考 equation.m 用中心差分生成 Jacobian，避免抄写大型符号 Jacobian。 */
static void numerical_jacobian(const void *context, const double *x,
                               double lambda,
                               double jacobian[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS],
                               DcResidualFunction residual, int n)
{
    double plus[DC_MAX_UNKNOWNS], minus[DC_MAX_UNKNOWNS];
    double fp[DC_MAX_UNKNOWNS], fm[DC_MAX_UNKNOWNS];
    for (int column = 0; column < n; ++column) {
        memcpy(plus, x, sizeof(double) * (size_t)n);
        memcpy(minus, x, sizeof(double) * (size_t)n);
        const double h = 1.0e-6 * fmax(1.0, fabs(x[column]));
        plus[column] += h;
        minus[column] -= h;
        residual(context, plus, lambda, fp);
        residual(context, minus, lambda, fm);
        for (int row = 0; row < n; ++row) jacobian[row][column] = (fp[row] - fm[row]) / (2.0 * h);
    }
}

static void schmitt1_residual(const void *context, const double *x,
                              double lambda, double *f)
{
    (void)context; (void)lambda;
    const double is = 1e-16, af = .99, ar = .5, n = 38.78;
    const double fe1 = -is / af * (safe_exp(-n * (x[1] - x[4])) - 1.0);
    const double fc1 = -is / ar * (safe_exp(-n * (x[0] - x[4])) - 1.0);
    const double fe2 = -is / af * (safe_exp(-n * (x[1] - x[3])) - 1.0);
    const double fc2 = -is / ar * (safe_exp(-n * (x[2] - x[3])) - 1.0);
    const double ie1 = fe1 - ar * fc1, ic1 = fc1 - af * fe1;
    const double ie2 = fe2 - ar * fc2, ic2 = fc2 - af * fe2;
    f[0] = (x[0]-x[3])/10000.0 + (x[0]-x[5])/2000.0 + ic1;
    f[1] = x[1]/100.0 + ie1 + ie2;
    f[2] = (x[2]-x[5])/1000.0 + ic2;
    f[3] = (x[3]-x[0])/10000.0 - ic2 - ie2;
    f[4] = x[4] - 1.5;
    f[5] = x[5] - 10.0;
    f[6] = (x[5]-x[0])/2000.0 + (x[5]-x[2])/1000.0 + x[6];
    f[7] = x[7] - ic1 - ie1;
}
static void schmitt1_jacobian(const void *c,const double*x,double l,double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{ memset(j,0,sizeof(double)*DC_MAX_UNKNOWNS*DC_MAX_UNKNOWNS); numerical_jacobian(c,x,l,j,schmitt1_residual,8); }

static void schmitt2_residual(const void *context, const double *x,
                              double lambda, double *f)
{
    (void)context; (void)lambda;
    const double is = 1e-16, af = .99, ar = .5, n = 38.78;
    const double fe1 = -is / af * (safe_exp(-n * (x[0] - x[4])) - 1.0);
    const double fc1 = -is / ar * (safe_exp(-n * (x[1] - x[4])) - 1.0);
    const double fe2 = -is / af * (safe_exp(-n * (x[0] - x[3])) - 1.0);
    const double fc2 = -is / ar * (safe_exp(-n * (x[2] - x[3])) - 1.0);
    const double ie1 = fe1 - ar * fc1, ic1 = fc1 - af * fe1;
    const double ie2 = fe2 - ar * fc2, ic2 = fc2 - af * fe2;
    f[0] = x[0]/100.0 + ie1 + ie2;
    f[1] = (x[1]-x[3])/10000.0 + (x[1]-x[5])/1500.0 + ic1;
    f[2] = (x[2]-x[5])/1000.0 + ic2;
    f[3] = (x[3]-x[1])/10000.0 + x[3]/1e6 - ic2 - ie2;
    f[4] = (x[4]-x[5])/5000.0 + x[4]/1250.0 - ic1 - ie1;
    f[5] = x[5] - 10.0;
    f[6] = (x[5]-x[1])/1500.0 + (x[5]-x[2])/1000.0 + (x[5]-x[4])/5000.0 + x[6];
}
static void schmitt2_jacobian(const void *c,const double*x,double l,double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{ memset(j,0,sizeof(double)*DC_MAX_UNKNOWNS*DC_MAX_UNKNOWNS); numerical_jacobian(c,x,l,j,schmitt2_residual,7); }

static void chua_currents(const double *x, double *ic, double *ie)
{
    const int c[4]={3,9,4,11}, b[4]={0,5,0,6}, e[4]={7,7,8,8};
    for(int i=0;i<4;++i){ const double ef=safe_exp(-38.7766*(x[e[i]]-x[b[i]]))-1.0;
        const double er=safe_exp(-38.7766*(x[c[i]]-x[b[i]]))-1.0;
        ic[i]=1e-9*ef-2e-9*er; ie[i]=-1.01e-9*ef+1e-9*er; }
}
static void chua_residual(const void *context,const double*x,double lambda,double*f)
{
    (void)context;(void)lambda; double ic[4],ie[4]; chua_currents(x,ic,ie);
    f[0]=x[0]-x[1]-10; f[1]=(x[1]-x[2])/10000-x[16]; f[2]=(x[2]-x[1])/10000+(x[2]-x[9])/30000-x[17];
    f[3]=(x[3]-x[12])/4000+(x[3]-x[5])/30000+ic[0]; f[4]=(x[4]-x[13])/4000+(x[4]-x[6])/30000+ic[2];
    f[5]=(x[5]-x[3])/30000+x[5]/10100-ic[1]-ie[1]; f[6]=(x[6]-x[4])/30000+x[6]/10100-ic[3]-ie[3];
    f[7]=x[7]/500+ie[0]+ie[1]; f[8]=x[8]/500+ie[2]+ie[3]; f[9]=(x[9]-x[12])/4000+(x[9]-x[2])/30000+ic[1];
    f[10]=x[10]-x[2]-2; f[11]=(x[11]-x[13])/4000+(x[11]-x[10])/30000+ic[3]; f[12]=x[12]-12; f[13]=x[13]-12;
    f[14]=(x[12]-x[3])/4000+(x[12]-x[9])/4000+x[14]; f[15]=(x[13]-x[4])/4000+(x[13]-x[11])/4000+x[15];
    f[16]=x[0]/5000+x[16]-ic[0]-ie[0]-ic[2]-ie[2]; f[17]=(x[10]-x[11])/30000+x[17];
}
static void chua_jacobian(const void*c,const double*x,double l,double j[DC_MAX_UNKNOWNS][DC_MAX_UNKNOWNS])
{ memset(j,0,sizeof(double)*DC_MAX_UNKNOWNS*DC_MAX_UNKNOWNS); numerical_jacobian(c,x,l,j,chua_residual,18); }

typedef struct {
    const char *name;
    int dimension;
    int count;
    /* 可选 CSV 输出：路径回调与工作点回调共享同一份实验上下文。 */
    FILE *path_file;
    FILE *solution_file;
} Run;

/* 将每个被接受的 predictor-corrector 点写出，供绘图脚本原样使用。 */
static void record_path(int step, double arc_length, double lambda,
                        double arc_step, double local_error, const double *x,
                        void *data)
{
    Run *run = data;
    if (run->path_file == NULL) return;
    fprintf(run->path_file, "%d,%.17g,%.17g,%.17g,%.17g", step, arc_length,
            lambda, arc_step, local_error);
    for (int i = 0; i < run->dimension; ++i) fprintf(run->path_file, ",%.17g", x[i]);
    fputc('\n', run->path_file);
}
/* 输出每个穿过 lambda=1 后又经原始 Newton 复核的完整 MNA 工作点。 */
static void found(int index,const double*x,double norm,void *data)
{ Run*r=data; ++r->count; printf("%s solution %d: ||F||inf=%.3e\n",r->name,index+1,norm);
  if (r->solution_file != NULL) {
      fprintf(r->solution_file, "%d,%.17g", index + 1, norm);
      for (int i = 0; i < r->dimension; ++i) fprintf(r->solution_file, ",%.17g", x[i]);
      fputc('\n', r->solution_file);
  }
  for(int i=0;i<r->dimension;++i) printf("  x%-2d = % .10e\n",i+1,x[i]); }
static int run(const char*name,DcProblem*p,const double*a,int steps,const char *output_directory)
{ DcSolverOptions no=dc_solver_default_options(); DcHomotopyOptions ho=dc_homotopy_default_options(a); Run r={name,p->dimension,0,NULL,NULL}; int n=0;
  char path_name[512], solution_name[512];
  if (output_directory != NULL) {
      snprintf(path_name, sizeof(path_name), "%s/%s_path.csv", output_directory, name);
      snprintf(solution_name, sizeof(solution_name), "%s/%s_solutions.csv", output_directory, name);
      r.path_file = fopen(path_name, "w"); r.solution_file = fopen(solution_name, "w");
      if (r.path_file == NULL || r.solution_file == NULL) {
          fprintf(stderr, "Cannot create Homotopy CSV output for %s.\n", name);
          if (r.path_file != NULL) fclose(r.path_file);
          if (r.solution_file != NULL) fclose(r.solution_file);
          return -1;
      }
      fprintf(r.path_file, "step,arc_length,lambda,arc_step,local_error");
      fprintf(r.solution_file, "solution,residual_norm");
      for (int i = 0; i < p->dimension; ++i) {
          fprintf(r.path_file, ",x%d", i + 1);
          fprintf(r.solution_file, ",x%d", i + 1);
      }
      fputc('\n', r.path_file); fputc('\n', r.solution_file);
  }
  no.maximum_newton_iterations=120; ho.initial_arc_step=.01; ho.maximum_arc_step=.1; ho.maximum_path_steps=steps; ho.local_error_tolerance=1e-4;
  const bool ok=dc_fixed_point_homotopy_solve(p,&no,&ho,record_path,found,&r,&n);
  if (r.path_file != NULL) fclose(r.path_file);
  if (r.solution_file != NULL) fclose(r.solution_file);
  printf("%s: path=%s, solutions=%d\n",name,ok?"complete":"failed",n); return ok?n:-1; }
int main(int argc, char **argv)
{ const double a1[]={.5,.4799,.9047,.6099,.6177,.8594,.8055,.5767}; const double a2[]={.5828,.5862,.9258,.5751,.01,.8094,.6088}; const double a3[]={.0759,.054,.5308,.7792,.934,.1299,.5688,.4694,.0119,.3371,.1622,.7943,.3112,.5285,.1656,.602,.263,.6541};
  DcProblem p1={8,NULL,schmitt1_residual,schmitt1_jacobian,NULL},p2={7,NULL,schmitt2_residual,schmitt2_jacobian,NULL},p3={18,NULL,chua_residual,chua_jacobian,NULL};
  const char *output_directory = argc == 2 ? argv[1] : NULL;
  if (argc > 2) { fprintf(stderr, "Usage: %s [csv-output-directory]\n", argv[0]); return 1; }
  const int n1=run("Schmitt1",&p1,a1,10000,output_directory), n2=run("Schmitt2",&p2,a2,10000,output_directory), n3=run("Chua",&p3,a3,10000,output_directory);
  return n1==3 && n2==3 && n3==9 ? 0 : 1; }
