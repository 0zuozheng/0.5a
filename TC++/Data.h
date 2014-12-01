#include <fstream>
using namespace std;

//声明输入数据外部变量
//控制参数信息Control
extern int nodenum_c, enrich_nodenum_c, pipe_nodenum_c, elementnum_c, materialnum_c, loadnum_c, fixnum_c, elsetnum_c, nsetnum_c, sectionnum_c, initialnum_c, surfsetnum_c, surfnum_c,tablenum_c;

//材料信息Material
extern float *density_m, **diffusivity_m, *c_m, *rise_m, *m_m, **conduction_m;
extern char **name_m;

//结点信息Node
extern float **xyz_n, *t_n;
extern int *enrichorder_n;

//单元信息Element
extern int *material_e, **node_e, *NodeNum_e, *plan_e, *PointNum_e, *pipe_e, **surf_e;
extern float *t_e, ***LMN_Pipe_e;

//边界条件信息Boundary
extern int boundarynum_c,*num_b,**list_b,*table_b;
extern float *fix_b;

//工况信息 Loadcase
extern int incnum_c,inc_l, *deactive_l, maxelement_l;
extern float timeinc_l, totaltime_l;
extern char name_l[500];

//分组定义信息Set
extern char **name_elset,**name_nset;
extern int *type_elset, **list_elset, *num_elset, *type_nset, **list_nset, *num_nset;

//材料部件信息Section
extern char **elset_section,**material_section;

//面信息
extern char **name_surfset,**name_surf;
extern int *num_surfset, *type_surfset, **list_surfset, *num_surf, **elelist_surf, **facelist_surf;

//表格信息
extern char **name_tb;
extern int  *num_tb;
extern float **time_tb,**value_tb;
extern float **Table_step;

//Debug文件
extern FILE *DEBUG;

//水管信息
extern int ppnum_c, pipenum_c, *pts_p, **pId_p, *TbInlet_p, *TbOutlet_p, *TbFlow_p; // 水管条件 Pipe
extern float **xyz_pp, **tem_pp;

// iniArray.cpp
extern void Realloc2DArray_float(float ***arr, int k0, int k1, int k2);
extern void Alloc2DArray_float(float*** arr,int k1,int k2);
extern void Alloc2DArray_int  (int*** arr,int k1,int k2);
extern void Alloc2DArray_char  (char*** arr,int k1,int k2);
extern void	Alloc3DArray_float(float**** arr,int k1,int k2,int k3);
extern void Alloc3DArray_int  (int**** arr,int k1,int k2,int k3);
extern void FREE2D_int(int** arr,int num);
extern void FREE2D_float(float** arr,int num);
extern void FREE2D_char(char** arr,int num);
extern void FREE3D_int(int*** arr,int num1,int num2);
extern void FREE3D_float(float*** arr,int num1,int num2);