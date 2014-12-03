/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *   Assign Functions
 *
 *	 Author:          Z Zuo
 *   Last Modified :  Dec. 2014
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#include <string.h>
#include "Data.h"
extern FILE* log_check;
extern void Warning(char*,char* tmp2="",char* tmp3="");

int* list;	int len; /* Integer list and length of the list, temporary saved, should only be used here. */

/* Search a Set named setname, returns the Set Id ... */
int SearchSetName(char* setname,int type){
	int id = 0;
	switch(type){
		case 1:		// Type = 1 (node set): name_nset, e.g., Nset-1
			for (int i=0;i<nsetnum_c;i++)
				if (strcmp(setname,name_nset[i])==0)	id=i;
			break;
		case 2:		// Type = 2 (element set): name_elset, e.g., Elset-1
			for (int i=0;i<elsetnum_c;i++)
				if (strcmp(setname,name_elset[i])==0)	id=i;
			break;
		case 3:		// Type = 3 (surface set): name_surfset, e.g., _Surf-1_S4
			for (int i=0;i<surfsetnum_c;i++)
				if (strcmp(setname,name_surfset[i])==0)	id=i;
			break;
		case 4:		// Type = 4 (surface): name_surf, e.g., Surf-1
			for (int i=0;i<surfnum_c;i++)
				if (strcmp(setname,name_surf[i])==0)	id=i;
			break;
		case 10:	// Type = 10 (material set): name_m, e.g., MAT001
			for (int i=0;i<materialnum_c;i++)
				if (strcmp(setname,name_m[i])==0)		id=i;
			break;
		case 20:	// Type = 20 (table set): name_tb, e.g., AMP-1
			for (int i=0;i<tablenum_c;i++)
				if (strcmp(setname,name_tb[i])==0)		id=i;
			break;
	}
	return id;
}

/* Search the list of the set by the set id, returns to the list ... */
void SearchSetList(int id,int type,char* cmd="Default"){
	free(list);	// clear the list, to receive new list obtained in this func.
	if (strcmp(cmd,"Default")!=0)	id = SearchSetName(cmd,type);	// if a set name is given, first search the set id of the set name
	
	/* Type = 1 (Node set)  */
	if (type==1){	// All type_nset can be checked in AbqPDF/NSET.pdf...
		if (type_nset[id]==1){ // GENERATE : Each data line should give a first node, a last node and the increment in node numbers between these nodes. 
			len  = num_nset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_nset[id];i++){
				list[i] = list_nset[id][0] + i*list_nset[id][2];
			}
		}else if(type_nset[id]==2){ // ELSET : The nodes that define the elements that belong to this element set at the time this option is encountered will be assigned to the node set specified.
			len  = num_nset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int j=0;j<num_nset[id];j++){
				list[j] = list_nset[id][j];
			}
		}else if (type_nset[id]==0){ // NO GENERATE : List of nodes or node set labels to be assigned to this node set. Only previously defined node sets can be assigned to another node set.
			len  = num_nset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_nset[id];i++){
				list[i] = list_nset[id][i];
			}
		}
	}
	/* Type = 2 (Element set)  */
	if (type==2){	// All type_elset can be checked in AbqPDF/ELSET.pdf...
		if (type_elset[id]==1){	// Each data line should give a first element, a last element, and the increment in element numbers between these elements.
			len  = num_elset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_elset[id];i++){
				list[i] = list_elset[id][0] + i*list_elset[id][2];
			}
		}else if (type_elset[id]==0){ // List of elements or element set labels to be assigned to this element set. Only previously defined element sets can be assigned to another element set.
			len  = num_elset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_elset[id];i++){
				list[i] = list_elset[id][i];
			}
		}
	}
	/* Type = 3 (Surface set)  */
	if (type==3){	// All type_surfset can be checked in AbqPDF/ELSET.pdf...
		if (type_surfset[id]==0){ // Each data line should give a first element, a last element, and the increment in element numbers between these elements.
			len  = num_surfset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_surfset[id];i++){
				list[i] = list_surfset[id][i];
			}
		}if (type_surfset[id]==1){ // List of elements or element set labels to be assigned to this element set. Only previously defined element sets can be assigned to another element set.
			len  = num_surfset[id];
			list = (int*)calloc(len,sizeof(int));
			for (int i=0;i<num_surfset[id];i++){
				list[i] = list_surfset[id][0] + i*list_surfset[id][2];
			}
		}
	}
}

/* Assign initial temperatures to the node set ... */
void AssignInitial_T(char* nsetname, float iniT){
	int ele;
	int id = SearchSetName(nsetname,1);
	if (type_nset[id]==2){			// type_nset : ELSET
		ele = num_nset[id];			// 获得单元组编号
		SearchSetList(ele,2);		// 获得单元列表
		for (int i=0;i<len;i++)		// 为单元赋初值
			t_e[list[i]] = iniT;
	}else if(type_nset[id]==1){		// type_nset : GENERATE
		SearchSetList(id,1);
		for (int i=0;i<len;i++)		// 为结点赋初值
			t_n[list[i]] = iniT;
	}else		
		Warning("The node set used in Initial Condition should be the elset format");	
}

/* Assign material id onto the section elements  ... */
void AssignMaterial(int id,char* MaterialName){
	int MatID = SearchSetName(MaterialName,10);		//查询材料组号
	int SetID = SearchSetName(elset_section[id],2);	//查询单元组号
	SearchSetList(SetID,2);

	for (int i=0;i<len;i++){
		material_e[list[i]] = MatID;
	}
}

/* Store element list and face list into the surface ... */
void AssignSurface(int SurfId,char* SetName,int face){
	int SetID = SearchSetName(SetName,3);
	SearchSetList(SetID,3);

	int num_old = num_surf[SurfId];

	elelist_surf[SurfId]  = (int*)realloc(elelist_surf[SurfId],(num_old+len)*sizeof(int));
	facelist_surf[SurfId] = (int*)realloc(facelist_surf[SurfId],(num_old+len)*sizeof(int));

	for (int i=0;i<len;i++){
		elelist_surf[SurfId][num_old+i] = list[i];
		facelist_surf[SurfId][num_old+i] = face;
	}num_surf[SurfId] += len;
}