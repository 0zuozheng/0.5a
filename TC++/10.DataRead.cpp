/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *		Input Reading Functions
 *
 *	 Author:          Z Zuo
 *   Last Modified :  Dec. 2014
 *
 *  All the file pointer(*fp) are already defined in 01.File.cpp.
 *  
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include "DataRead.h"
extern void Free_Read();

/* Read and Write the line from *in into *out, and extract the string into *tmp... */
void RnWinp(FILE *in,FILE *out,char *tmp){
	fgets(tmp,500,in);
	fprintf(out,"%s",tmp);
	if (tmp[strlen(tmp)-1] != 10)	fprintf(out,"\n");
}

/* Output the file head into the *out... */
void OutHead(FILE *out, char *head){
	fprintf(out,"************\n");
	fprintf(out,"** %s\n",head);
	fprintf(out,"************\n");
	fprintf(log_check,"\n#### %s Information\n",head);
}

/* Read Node.inp and extract the nodes info (variables with _n as postfix, i.e, xyz_n...) */
void Read_Node(FILE *in)	
{	int id;	char tmp[500];	char *p;
	OutHead(out,"PARTS (Nodes)");
	fprintf(log_check,"%10s%15s%15s%15s\n","NodeID","X","Y","Z");
	fscanf (in,"**Node Num =%d\n",&nodenum_c);	RnWinp(in,out,tmp);
	enrich_nodenum_c = pipe_nodenum_c = nodenum_c;

	Alloc2DArray_float(&xyz_n,nodenum_c,3);	enrichorder_n = (int*)calloc(nodenum_c,sizeof(float));

	for(int i=0;i<nodenum_c;i++){
		RnWinp(in,out,tmp);		p = strtok(tmp,",");	sscanf(p,"%d",&id);		id--;	// Node Id - 1
		for(int j=0;j<3;j++){
			p = strtok(NULL,",");	sscanf(p,"%f",&xyz_n[id][j]);
		}
		fprintf(log_check,"%10d%15.3f%15.3f%15.3f\n",id+1,xyz_n[id][0],xyz_n[id][1],xyz_n[id][2]);
	}fclose(in);
}

/* Read Element.inp and extract the ele info (variables with _e as postfix, i.e, material_e...) */
void Read_Element(FILE *in) 
{	int num,id;	char cmd[10][500], tmp[500], *p;
	OutHead(out,"PARTS (Elements)");
	fprintf(log_check,"  EleID| * * * * * * * *  C o n n e c t i v i y * * * * * * * * |\n");      
	fscanf (in,"**Element Num =%d\n",&elementnum_c);
	LineExtract(in,cmd,&num);	// Extract the keywords
	KeywordElement(cmd,num,0);	// Parse the keywords

	material_e=(int*)calloc(elementnum_c,sizeof(int));		t_e       =(float*)calloc(elementnum_c,sizeof(float));
	Alloc2DArray_int(&node_e,elementnum_c,8);				Alloc3DArray_float(&LMN_Pipe_e,elementnum_c,2,3);
	NodeNum_e =(int*)calloc(elementnum_c,sizeof(int));		PointNum_e = (int*)calloc(elementnum_c,sizeof(int));
	plan_e = (int*)calloc(elementnum_c,sizeof(int));		pipe_e = (int*)calloc(elementnum_c,sizeof(int));

	for(int i=0;i<elementnum_c;i++){
		RnWinp(in,out,tmp);		p = strtok(tmp,",");	sscanf(p,"%d,",&id);
		fprintf(log_check,"%5d",id);	id--;	// Element Id - 1

		NodeNum_e[id] = 8;		// Number of nodes per ele
		plan_e[id]  = 5;		// Scheme of integration point layout

		if(i==0) {pipe_e[id]  = 1;	NodeNum_e[id] = 16;}	// ZUOZUO, for test, only pipe in ele 0

		if(pipe_e[id])	node_e[id] = (int*)realloc(node_e[id],16*sizeof(int));

		for(int j=0;j<8;j++) {	// Output the element topo/connectivity
			p = strtok(NULL,",");	sscanf(p,"%d",&node_e[id][j]);
			fprintf(log_check,"%7d",node_e[id][j]);
			node_e[id][j]--;
		}	fprintf(log_check,"|\n");
	}	fclose(in);
}

/* Read Set.inp and extract the set info (variables with _set as postfix, i.e, list_elset...) */
void Read_Set(FILE *in){	
	int num, ElsetId = -1, NsetId = -1;	// If first set exists, then ElsetId++ -> ElsetId=0
	char cmd[20][500];	// Attention: 20 keywords as maximum
	OutHead(out,"PARTS (Element sets & Node sets)");
	fscanf (in,"**Element Set Num =%d\n",&elsetnum_c);
	fscanf (in,"**Node Set Num =%d\n",&nsetnum_c);

	Alloc2DArray_char(&name_elset,elsetnum_c,500);			type_elset = (int*)calloc(elsetnum_c,sizeof(int));
	list_elset = (int**)calloc(elsetnum_c,sizeof(int*));	num_elset  = (int*)calloc(elsetnum_c,sizeof(int));
	Alloc2DArray_char(&name_nset,nsetnum_c,500);			type_nset = (int*)calloc(nsetnum_c,sizeof(int));
	list_nset = (int**)calloc(nsetnum_c,sizeof(int*));		num_nset  = (int*)calloc(nsetnum_c,sizeof(int));

	while(!feof(in)){
		LineExtract(in,cmd,&num);
		KeywordSet(cmd,num,&ElsetId,&NsetId);
	}

	for(int i=0;i<elsetnum_c;i++){
		if (type_elset[i]==1)	fprintf(log_check,"  Element Set Name: %s , Element from %d to %d , in total of %d .\n",name_elset[i],list_elset[i][0]+1,list_elset[i][1]+1,num_elset[i]);
		if (type_elset[i]==0){
			fprintf(log_check,"  Element Set Name: %s\n",name_elset[i]);
			for(int j=0;j<num_elset[i];j++)	fprintf(log_check,"\t%d",list_elset[i][j]+1);
			fprintf(log_check,"\n");
		}
	}
	
	for(int i=0;i<nsetnum_c;i++){
		if (type_nset[i]==1)	fprintf(log_check,"  Node Set Name: %s , Node from %d to %d , in total of %d .\n",name_nset[i],list_nset[i][0]+1,list_nset[i][1]+1,num_nset[i]);
		if (type_nset[i]==2)	fprintf(log_check,"  Node Set Name: %s , is generated by the element set %s\n",name_nset[i],name_elset[num_nset[i]]);
		if (type_nset[i]==0){
			fprintf(log_check,"  Element Set Name: %s\n",name_nset[i]);
			for(int j=0;j<num_nset[i];j++)	fprintf(log_check,"\t%d",list_nset[i][j]+1);
			fprintf(log_check,"\n");
		}
	}fclose(in);
}

/* Read Section.inp and extract the set info (variables with _section as postfix, i.e, material_section...) */
void Read_Section(FILE *in){
	int num;	char cmd[20][500];	// Attention: 20 keywords as maximum
	OutHead(out,"PARTS (sections)");
	fscanf (in,"**Section Num = %d\n",&sectionnum_c);

	Alloc2DArray_char(&elset_section,sectionnum_c,500);		
	Alloc2DArray_char(&material_section,sectionnum_c,500);
	
	for (int i=0;i<sectionnum_c;i++){
		LineExtract(in,cmd,&num);	// Extract keywords
		KeywordSection(cmd,num,i);	// Parse keywords
		fprintf(log_check,"  Section %2d Element Set Name = '%s' Material Name = '%s'\n",i+1,elset_section[i],material_section[i]);
	}	fclose(in);
}

/* Read Material.inp and extract the mat info (variables with _m as postfix, i.e, density_m...) */
void Read_Material(FILE *in)
{	char cmd[20][500];	int num,id = -1;
	OutHead(out,"MATERIAL");
	fscanf (in,"**Material Num = %d\n",&materialnum_c);

	Alloc2DArray_char(&name_m,materialnum_c,500);	
	density_m = (float*)calloc(materialnum_c,sizeof(float));		Alloc2DArray_float(&diffusivity_m,materialnum_c,3);
	rise_m    = (float*)calloc(materialnum_c,sizeof(float));		Alloc2DArray_float(&conduction_m,materialnum_c,3);
	c_m = (float*)calloc(materialnum_c,sizeof(float));				m_m = (float*)calloc(materialnum_c,sizeof(float));

	while(!feof(in)){
		LineExtract(in,cmd,&num);
		KeywordMaterial(cmd,num,&id);
	}
	for(int i=0;i<materialnum_c;i++){
		rise_m[i] = 0.0;	m_m[i]=0.16;	// ZUOZUO No consideration of hydration heat
		conduction_m[i][2] = conduction_m[i][1] = conduction_m[i][0];	// Isotropy. But can be developed into Anisotropy...
		diffusivity_m[i][0] = diffusivity_m[i][1] = diffusivity_m[i][2] = conduction_m[i][0]/(density_m[i]*c_m[i]);

		fprintf(log_check," Material Name = \t%s\n Density (kg*m^-3):%10.3f\n Specific heat(kJ/(kg*oC)):%10.3f\n Max adiabatic temperature rise (oC):%10.3f\n Adiavatic rate (-):%10.3f\n",name_m[i],density_m[i],c_m[i],rise_m[i],m_m[i]);
		fprintf(log_check," Diffusivity coefficient X (m^2/d):%10.3f\n Diffusivity coefficient Y (m^2/d):%10.3f\n Diffusivity coefficient Z (m^2/d):%10.3f\n",diffusivity_m[i][0],diffusivity_m[i][1],diffusivity_m[i][2]);
		fprintf(log_check," Conductivity X (kJ/(m*d*oC)):%10.3f\n Conductivity Y (kJ/(m*d*oC)):%10.3f\n Conductivity Z (kJ/(m*d*oC)):%10.3f\n",conduction_m[i][0],conduction_m[i][1],conduction_m[i][2]);
	}	fclose(in);
}

/* Assign the material onto the section... */
void Assign_Section(){
	OutHead(out,"Assigning Material");
	for (int i=0;i<sectionnum_c;i++)
		AssignMaterial(i,material_section[i]);	

	fprintf(log_check," MatID");
	for (int i=0;i<elementnum_c;i++){
		if(i%10==0)	fprintf(log_check,"\n%8d:",i+1);
		fprintf(log_check,"%8d",material_e[i]+1);
	}fprintf(log_check,"\n");
}

/* Read Surface.inp and extract the surf info (variables with _surf as postfix, i.e, list_surfset...) */
void Read_Surface(FILE *in){
	char cmd[20][500]; // Attention: 20 keywords as maximum
	int num, SetId = -1, SurfId = -1;	// If first surf exists, then SurfId++ -> SurfId=0
	OutHead(out,"Surface");
	fscanf (in,"**Set Num = %d Surface Num = %d\n",&surfsetnum_c,&surfnum_c);

	Alloc2DArray_char(&name_surfset,surfsetnum_c,500);		// Surfset Name, e.g., _Surf-1_S4
	num_surfset  = (int*)calloc(surfsetnum_c,sizeof(int));
	type_surfset = (int*)calloc(surfsetnum_c,sizeof(int));
	list_surfset = (int**)calloc(surfnum_c,sizeof(int*));	// Only Ele Id, Face Id in the name

	Alloc2DArray_char(&name_surf,surfnum_c,500);
	num_surf      = (int*)calloc(surfnum_c,sizeof(int));
	elelist_surf  = (int**)calloc(surfnum_c,sizeof(int**));	// Ele Id
	facelist_surf = (int**)calloc(surfnum_c,sizeof(int**));	// Face Id

	while(!feof(in)){
		LineExtract(in,cmd,&num);
		KeywordSurface(cmd,num,&SetId,&SurfId);
	}

	fprintf(log_check,"** Attention: The face number definition is different with Abaqus.\n 0 - > 4;\n 1 - > 6;\n 2 - > 5;\n 3 - > 3;\n 4 - > 2;\n 5 - > 1;\n");
	for(int i=0;i<surfnum_c;i++){
		fprintf(log_check,"  Surface Name = %s: Total Num = %d",name_surf[i],num_surf[i]);
		for(int j=0;j<num_surf[i];j++){
			if(j%8==0)	fprintf(log_check,"\n");
			fprintf(log_check,"\t%d:%d",elelist_surf[i][j]+1,facelist_surf[i][j]);
		}fprintf(log_check,"\n");
	}	fclose(in);
}

/* Output control data into CheckDAT.log ...) */
void welcome(){
	fprintf(log_check,"\n#### Control Information\n");
	fprintf(log_check," Number of nodes:%d\n",nodenum_c);
	fprintf(log_check," Number of elements:%d\n",elementnum_c);
	fprintf(log_check," Number of materials:%d\n",materialnum_c);
}

/* Read Initial.inp and extract the initial condition... */
void Read_Initial(FILE *in){
	int num;	char cmd[10][500];		float iniT;
	OutHead(out,"INITIAL CONDITION");
	fscanf (in,"**Initial Condition Num =%d\n",&initialnum_c);

	for(int i=0;i<initialnum_c;i++){
		LineExtract(in,cmd,&num);	// 提取关键字
		KeywordInitial(cmd,num,i);	// 解析关键字

		LineExtract(in,cmd,&num);	// 提取关键字
		sscanf(cmd[1],"%f",&iniT);	// cmd[0]:组的名称; cmd[1]:初始温度
		fprintf(out,"%s, %f\n",cmd[0],iniT);
		AssignInitial_T(cmd[0],iniT);
	}

	fprintf(log_check,"Elements Temperature");
	for(int i=0;i<elementnum_c;i++){
		if (i%10==0)	fprintf(log_check,"\n%8d:",i+1);
		fprintf(log_check,"%15.3f",t_e[i]);
	}fclose(in);
}

/* Check the table is given ordered by time... */
void Check_Table(){
	for (int i=0;i<tablenum_c;i++){
		if (time_tb[i][0]!=0)
			Warning("The initial time of a tabular amplitude curve is not zero.");
		for(int j=0;j<num_tb[i]-1;j++){
			if (time_tb[i][j+1] <= time_tb[i][j])
				Warning("A tabular amplitude curve is defined with time not monotonically increasing");
		}
	}
}

/* Read Table.inp and extract the table (variables with _tb as postfix, i.e, time_tb...) */
void Read_Table(FILE* in){
	int num,id = -1;	char cmd[20][500];
	OutHead(out,"Table/Amplitude");
	fscanf (in,"**Table Num = %d\n",&tablenum_c);

	Alloc2DArray_char(&name_tb,tablenum_c,500);		time_tb  = (float**)calloc(tablenum_c,sizeof(float*));
	num_tb  = (int*)calloc(tablenum_c,sizeof(int));	value_tb = (float**)calloc(tablenum_c,sizeof(float*));

	while(!feof(in)){
		LineExtract(in,cmd,&num);
		KeywordTable(cmd,num,&id);
	}
	for (int i=0;i<tablenum_c;i++){
		fprintf(log_check,"Table Name : %s\n",name_tb[i]);
		for(int j=0;j<num_tb[i];j++){
			fprintf(log_check,"%10.2f%10.2f\n",time_tb[i][j],value_tb[i][j]);
		}
	}
	fclose(in);
	Check_Table();
}

/* Read Pipe.inp and extract the pipe point info (variables with _pp as postfix, i.e, xyz_pp...)
   and extract the pipe topo info (variables with _p as postfix, i.e, pId_p...) */
void Read_Pipe(FILE *fp){
	char cmd[20][500];		int num;
	fscanf(fp,"**Pipe Point Num =%d\n*Pipe Point\n",&ppnum_c);
	Alloc2DArray_float(&xyz_pp, ppnum_c, 3);	Alloc2DArray_float(&tem_pp, ppnum_c, 1);
	
	fprintf(log_check,"\n\n#### Pipe Geometry\n Total Pipe Points Number = %d\n",ppnum_c);
	for (int ij,i=0; i<ppnum_c; i++){		
		fscanf(fp,"%d %f %f %f\n",&ij,&xyz_pp[i][0],&xyz_pp[i][1],&xyz_pp[i][2]);
		fprintf(log_check,"\t%.3f\t%.3f\t%.3f\n",xyz_pp[i][0],xyz_pp[i][1],xyz_pp[i][2]);
	}

	// Pipe topo (pipe points included in each pipe segment) is parsed by the prog below...
	fscanf(fp,"**Pipe Segment Num =%d\n*Pipe Segment\n",&pipenum_c);
	fprintf(log_check,"\n Total Pipe Segments Number = %d\n",pipenum_c);

	pts_p = (int*)calloc(pipenum_c,sizeof(int));		pId_p = (int**)calloc(pipenum_c,sizeof(int*));		
	TbInlet_p = (int*)calloc(pipenum_c,sizeof(int*));	TbOutlet_p = (int*)calloc(pipenum_c,sizeof(int*));	TbFlow_p  = (int*)calloc(pipenum_c,sizeof(int*));

	for (int i=0;i<pipenum_c;i++){		
		LineExtract(fp,cmd,&num);
		pts_p[i] = num-1;
		pId_p[i] = (int*)calloc(num-1,sizeof(int));
		fprintf(log_check," Pipe Segment No.%d has a total of %d points.\n",i+1,num-1);
		for(int j=1;j<num;j++){
			pId_p[i][j-1] = atoi(cmd[j])-1;
			fprintf(log_check,"\t%d", pId_p[i][j-1]+1);
		}fprintf(log_check,"\n");
	}
	fclose(fp);
}

/* Read Cooling.inp and extract the cooling condition */
void Read_Cooling(FILE *fp){	char cmd[20][500];	int num;
	OutHead(out,"Cooling Condition");
	while(!feof(fp)){
		LineExtract(fp,cmd,&num);
		KeywordCooling(cmd,num);
	}
}

/* ReadData() calls all the read sub above, which is called by main()... */
void ReadData(){
	printf("Reading model information form INP...\n");
	printf(" *Reading node ...\n");		Read_Node(inp_node);
	printf(" *Reading element ...\n");	Read_Element(inp_element);
	printf(" *Reading set ...\n");		Read_Set(inp_set);
	printf(" *Reading section ...\n");	Read_Section(inp_section);
	printf(" *Reading material ...\n");	Read_Material(inp_material);
	printf(" *Assigning material ...\n");	Assign_Section();
	welcome();
	printf(" *Reading surface...\n");	Read_Surface(inp_surface);
	printf(" *Reading table ...\n");	Read_Table(inp_table);
	printf(" *Reading initial condition ...\n");	Read_Initial(inp_initial);
	printf(" *Reading pipe information ...\n");		Read_Pipe(inp_pipe);
	printf(" *Reading cooling condition ...\n");	Read_Cooling(inp_cooling);
	Free_Read();
	printf("Reading model data Finished!\n\n");
}