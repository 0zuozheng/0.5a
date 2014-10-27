#include <stdio.h>
#include <fstream>
using namespace std;

FILE *log_answer, *log_solver, *log_element, *log_N, *log_check, *log_transform;	//数据输出校核文件流
FILE *inp_material, *inp_element, *inp_node, *inp_pipe, *inp_loadcase, *inp_table, *inp_set, *inp_section, *inp_initial, *inp_surface, *inp_step;	//数据输入文件流
FILE *out;
FILE *DEBUG;

extern char* DTFS;
extern char *Prog_Version,*Solver_Version,*Dev_Version,*Copy_Right,*DTFS;

void File(){
	/* Log File */
	log_solver = fopen("Solver.log","w");
	fprintf(log_solver,"%s\n*** This log file includes the A,B and X of the linear equations in each increment.\n\n",DTFS);

	log_element = fopen("Element.log","w");
	fprintf(log_element,"%s\n*** This log file includes the H,P,Q and Q3 of each element in each loadcase.\n\n",DTFS);

	log_transform = fopen("Transform.log","w");
	fprintf(log_transform,"%s\n*** This log file includes the spacial transforming matrix.\n\n",DTFS);

	log_N = fopen("N_shape.log","w");
	fprintf(log_N,"%s\n*** This log file includes the shape function matrix.\n\n",DTFS);

	log_check = fopen("CheckDAT.log","w");
	fprintf(log_check,"%s\n\n",DTFS);
	fprintf(log_check,"%s\n%s\n%s\n%s\n",Prog_Version,Solver_Version,Dev_Version,Copy_Right);

	/* Input File */
	inp_material    = fopen("Material.inp","r");
	inp_node        = fopen("Node.inp","r");
	inp_element     = fopen("Element.inp","r");
	inp_loadcase	= fopen("Loadcase.inp","r");
	inp_table		= fopen("Table.inp","r");
	inp_set			= fopen("Set.inp","r");
	inp_section		= fopen("Section.inp","r");
	inp_surface		= fopen("Surface.inp","r");
	inp_step		= fopen("Step.inp","r");
	inp_initial     = fopen("Initial.inp","r");
	inp_pipe		= fopen("Pipe.inp","r");
	fclose(inp_pipe);

	/* Output Job.inp file for abaqus */
	out				= fopen("Job.inp","w");

	/* Debug File */
	DEBUG = fopen("Debug.log","w");
}

void restart_files(){
	fclose(out);				out	= fopen("Job.inp","a");
	fclose(log_check);			log_check = fopen("CheckDAT.log","a");
	fclose(DEBUG);				DEBUG = fopen("Debug.log","a");
	fclose(log_solver);			log_solver = fopen("Solver.log","a");
	fclose(log_element);		log_element = fopen("Element.log","a");
	fclose(log_transform);		log_transform = fopen("Transform.log","a");
}

void close_files(){
	printf("\n***************\n*\n*    Mission Complete!\n*\n***************\n");
	fclose(inp_step);
	fclose(log_check);
	fclose(DEBUG);
	fclose(log_solver);
	fclose(log_element);
	fclose(log_transform);
	fclose(out);
}