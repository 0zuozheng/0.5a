#include<stdlib.h>
// 所有全局变量, 如欲引用，请引用Data.h
extern void FREE2D_int(int** arr,int num);
extern void FREE2D_float(float** arr,int num);
extern void FREE2D_char(char** arr,int num);
extern void FREE3D_int(int*** arr,int num1,int num2);
extern void FREE3D_float(float*** arr,int num1,int num2);
extern int *list;

int nodenum_c, elementnum_c, materialnum_c, loadnum_c, fixnum_c, elsetnum_c, nsetnum_c, sectionnum_c, initialnum_c, surfsetnum_c, surfnum_c;//控制参数信息Control

float *density_m, **diffusivity_m, *c_m, *rise_m, *m_m, **conduction_m;//材料信息Material
char **name_m;

float **xyz_n, *t_n;//结点信息Node

int *material_e, **node_e, *class_e, *plan_e;//单元信息Element
float *t_e;
int **surf_e;//单元面信息ElementFace
	
int boundarynum_c=0,*num_b,**list_b,*table_b;	//边界条件 Boundary
float *fix_b;
//int *node_b;
//float *fix_b;

int incnum_c,inc_l, *deactive_l, maxelement_l;//工况信息 Loadcase
float timeinc_l, totaltime_l; 
char name_l[500];

int tablenum_c;

int pipenum_c,*pts_p;
float ***xyz_p;

char **name_elset, **name_nset;		// 组信息 Set
int *type_elset, **list_elset, *num_elset, *type_nset, **list_nset, *num_nset;

char **elset_section, **material_section;	// 材料部件信息 Section

char **name_surfset,**name_surf;		// 面信息
int *num_surfset, *type_surfset, **list_surfset, *num_surf, **elelist_surf, **facelist_surf;

char **name_tb;	// 表格信息
int *num_tb;
float **time_tb,**value_tb;

float **Table_step;//step

void Free_Read(){	// 读入之后需要删除掉的变量
	FREE2D_char(name_m,materialnum_c);
	FREE2D_char(name_surfset,surfsetnum_c);
	free(num_surfset);
	free(type_surfset);
	FREE2D_int(list_surfset,surfsetnum_c);
}

void FreeVariables(){	// 程序运行完需要删除掉的变量

	FREE2D_float(xyz_n,nodenum_c);			free(t_n);
	free(material_e);   free(t_e);			FREE2D_int(node_e,elementnum_c);
	free(plan_e);		free(class_e);		FREE2D_int(surf_e,elementnum_c);
	free(deactive_l);
	
	// 材料信息
	FREE2D_float(diffusivity_m,materialnum_c);	free(density_m);	free(c_m);			
	FREE2D_float(conduction_m,materialnum_c);	free(rise_m);		free(m_m);			

	// 组信息(单元、结点)
	FREE2D_char(name_elset,elsetnum_c);		free(type_elset);	
	FREE2D_int(list_elset,elsetnum_c);		free(num_elset);
	FREE2D_char(name_nset,nsetnum_c);		free(type_nset);	
	FREE2D_int(list_nset,nsetnum_c);		free(num_nset);

	// Section信息
	FREE2D_char(elset_section,sectionnum_c);
	FREE2D_char(material_section,sectionnum_c);

	// 单元面信息		
	FREE2D_char(name_surf,surfnum_c);		free(num_surf);
	FREE2D_int(elelist_surf,surfnum_c);		FREE2D_int(facelist_surf,surfnum_c);

	// 表格信息
	FREE2D_char(name_tb,tablenum_c);		free(num_tb);
	FREE2D_float(time_tb,tablenum_c);		FREE2D_float(value_tb,tablenum_c);

	// 边界信息
	free(num_b);	free(table_b);		free(fix_b);
	FREE2D_int(list_b,boundarynum_c);

	// step
	FREE2D_float(Table_step,tablenum_c);

	// free list
	free(list);
}