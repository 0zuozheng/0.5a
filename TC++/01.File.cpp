/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *		FILE Pointer Declaration
 *
 *	 Author:          Z Zuo
 *   Last Modified :  Dec. 2014
 *  
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <stdio.h>
extern char *DTFS, *Prog_Version, *Solver_Version, *Dev_Version, *Copy_Right;

/* Log FILE Pointer */
FILE *log_answer, *log_solver, *log_element, *log_N, *log_check, *log_transform, *log_bound, *log_enrich;

/* Input FILE Pointer */
FILE *inp_material, *inp_element, *inp_node, *inp_pipe, *inp_loadcase, *inp_table, *inp_set, *inp_section, *inp_initial, *inp_surface, *inp_step, *inp_cooling;

/* Output FILE Pointer */
FILE *out;			FILE *DEBUG;

/* Open FILE Pointer, called by main()... */
void File(){
	log_solver = fopen("Solver.log","w");
	fprintf(log_solver,"%s\n*** This log file includes the A,B and X of the linear equations in each increment.\n\n",DTFS);

	log_element = fopen("Element.log","w");
	fprintf(log_element,"%s\n*** This log file includes the H,P,Q and Q3 of each element in each loadcase.\n\n",DTFS);

	log_transform = fopen("Transform.log","w");
	fprintf(log_transform,"%s\n*** This log file includes the spacial transforming matrix.\n\n",DTFS);

	log_N = fopen("N_shape.log","w");
	fprintf(log_N,"%s\n*** This log file includes the shape function matrix.\n\n",DTFS);

	log_enrich = fopen("Enrich.inp","w");
	fprintf(log_N,"*** This log file includes the model that contains pipe nodes and integration points.\n\n");

	log_bound = fopen("Bound.log","w");
	fprintf(log_bound,"%s\n*** This log file includes the boundary condition.\n\n",DTFS);

	log_check = fopen("CheckDAT.log","w");
	fprintf(log_check,"%s\n\n%s\n%s\n%s\n%s\n",DTFS,Prog_Version,Solver_Version,Dev_Version,Copy_Right);

	/* Input File Pointer */
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
	inp_cooling		= fopen("Cooling.inp","r");

	/* Output Job.inp file for abaqus */
	out				= fopen("Job.inp","w");

	/* Debug File, which can be called everywhere in this program... */
	DEBUG = fopen("Debug.log","w");
}

/* Restart FILE Pointer when another calc step begins... */
void restart_files(){
	fclose(out);				out	= fopen("Job.inp","a");	// Format = "Additional"
	fclose(log_check);			log_check = fopen("CheckDAT.log","a");
	fclose(DEBUG);				DEBUG = fopen("Debug.log","a");
	fclose(log_solver);			log_solver = fopen("Solver.log","a");
	fclose(log_element);		log_element = fopen("Element.log","a");
	fclose(log_transform);		log_transform = fopen("Transform.log","a");
}

/* Close all FILE Pointers ... */
void close_files(){
	printf("\n***************\n*\n*    Mission Complete!\n*\n***************\n");	
	fclose(log_check);	fclose(log_solver);	fclose(log_element);	fclose(log_transform);
	fclose(inp_step);	fclose(DEBUG);		fclose(out);
}