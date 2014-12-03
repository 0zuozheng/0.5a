/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *   Global Variables Declaration
 *
 * (all the global vars can be called by Data.h)
 *
 *	 Author:          Z Zuo
 *   Last Modified :  Dec. 2014
 *  
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include<stdlib.h>
extern void FREE2D_int(int** arr,int num);
extern void FREE2D_float(float** arr,int num);
extern void FREE2D_char(char** arr,int num);
extern void FREE3D_int(int*** arr,int num1,int num2);
extern void FREE3D_float(float*** arr,int num1,int num2);
extern int* list;

/* Controal Vars, with _c as postfix... */
int nodenum_c, pipe_nodenum_c, enrich_nodenum_c, elementnum_c, materialnum_c, loadnum_c, fixnum_c, elsetnum_c, nsetnum_c, sectionnum_c, initialnum_c, surfsetnum_c, surfnum_c, tablenum_c;

/* Material Vars, with _m as postfix... */
float *density_m, **diffusivity_m, *c_m, *rise_m, *m_m, **conduction_m;
char **name_m;

/* Node Vars, with _n as postfix... */
float **xyz_n, *t_n;
int *enrichorder_n;

/* Element Vars, with _e as postfix... */
int *material_e, **node_e, *class_e, *NodeNum_e, *plan_e, *PointNum_e, *pipe_e, **surf_e;
float *t_e, ***LMN_Pipe_e;

/* Boundary Vars, with _b as postfix... */
int boundarynum_c=0,*num_b,**list_b,*table_b;
float *fix_b;

/* Loadcase Vars, with _l as postfix... */
int inc_l, *deactive_l, maxelement_l;
float timeinc_l, totaltime_l; 
char name_l[500];

/* Pipe Vars, with _p as postfix; Pipe Point Vars, with _pp as postfix... */
int ppnum_c, pipenum_c;
float **xyz_pp, **tem_pp;
int *pts_p, **pId_p, *TbInlet_p, *TbOutlet_p, *TbFlow_p;

/* Set Vars, with _nset or _elset as postfix... */
char **name_elset, **name_nset;
int *type_elset, **list_elset, *num_elset, *type_nset, **list_nset, *num_nset;

/* Section Vars, with _section as postfix... */
char **elset_section, **material_section;

/* Surface Vars, with _surf and _surfset as postfix... */
char **name_surfset,**name_surf;
int *num_surfset, *type_surfset, **list_surfset, *num_surf, **elelist_surf, **facelist_surf;

/* Table Vars, with _tb as postfix... */
int *num_tb;	char **name_tb;
float **time_tb,**value_tb;

/* Step Vars, with _step as postfix... */
float **Table_step;

/* Free some variables which are used only in ReadData(); called when data reading is completed...*/
void Free_Read(){
	FREE2D_char(name_m,materialnum_c);
	free(num_surfset);
	free(type_surfset);
	FREE2D_char(name_surfset,surfsetnum_c);
	FREE2D_int(list_surfset,surfsetnum_c);
}

/* Free all global variables; called when computation is over (end of 40.Assembly.cpp)...*/
void FreeVariables()
{	// Node
	FREE2D_float(xyz_n, enrich_nodenum_c);			free(t_n);

	// Element
	free(material_e);   free(t_e);			FREE2D_int(node_e,elementnum_c);
	free(plan_e);		free(class_e);		FREE2D_int(surf_e,elementnum_c);

	// Loadcase
	free(deactive_l);
	
	// Material
	FREE2D_float(diffusivity_m,materialnum_c);	free(density_m);	free(c_m);			
	FREE2D_float(conduction_m,materialnum_c);	free(rise_m);		free(m_m);			

	// Set
	FREE2D_char(name_elset,elsetnum_c);		free(type_elset);	
	FREE2D_int(list_elset,elsetnum_c);		free(num_elset);
	FREE2D_char(name_nset,nsetnum_c);		free(type_nset);	
	FREE2D_int(list_nset,nsetnum_c);		free(num_nset);

	// Section
	FREE2D_char(elset_section,sectionnum_c);	FREE2D_char(material_section,sectionnum_c);

	// Surface
	FREE2D_char(name_surf,surfnum_c);		free(num_surf);
	FREE2D_int(elelist_surf,surfnum_c);		FREE2D_int(facelist_surf,surfnum_c);

	// Table
	FREE2D_char(name_tb,tablenum_c);		free(num_tb);
	FREE2D_float(time_tb,tablenum_c);		FREE2D_float(value_tb,tablenum_c);

	// Boundary
	free(num_b);	free(table_b);		free(fix_b);
	FREE2D_int(list_b,boundarynum_c);

	// Step
	FREE2D_float(Table_step,tablenum_c);

	// List (A temporary array used in 12.Assign.cpp)
	free(list);
}