#include "Data.h"
extern FILE *out,*log_check;
extern void LineExtract(FILE *in,char (*cmd)[500],int* num);
extern void Warning(char* tmp1,char *tmp2="",char* tmp3="");
extern int SearchSetName(char* setname,int type);
extern void SearchSetList(int id,int type,char* cmd="Default");
extern int *list, len;
float *BetaValue={0};	int *SinkTable={0},*SurfLoadType={0};	int SurfLoadId=1;	// 第三类边界条件

void deactivate(char* elset){
	int setid = SearchSetName(elset,2);
	SearchSetList(setid,2);
	fprintf(DEBUG,"  A total of %d elements in '%s' are DEactivated:\n",setid,elset);
	for (int i=0;i<len;i++){
		deactive_l[list[i]] = 1;
		fprintf(DEBUG,"\t%d",list[i]+1);
	}fprintf(DEBUG,"\n");	
}

void activate(char* elset){
	int setid = SearchSetName(elset,2);
	SearchSetList(setid,2);
	fprintf(DEBUG,"  A total of %d elements in '%s' are activated:\n",setid,elset);
	for (int i=0;i<len;i++){
		deactive_l[list[i]] = 0;
		fprintf(DEBUG,"\t%d",list[i]+1);
	}fprintf(DEBUG,"\n");	
}

void Sfilm(int setid,int load,float beta,int TableId){
	
	BetaValue = (float*)realloc(BetaValue,(load+1)*sizeof(float));
	SinkTable = (int*)realloc(SinkTable,(load+1)*sizeof(int));
	SurfLoadType = (int*)realloc(SurfLoadType,(load+1)*sizeof(int));
	
	SurfLoadType[load] = 3; //第三类边界条件	
	BetaValue[load] = beta;
	SinkTable[load] = TableId;

	for(int i=0;i<num_surf[setid];i++){
		int e = elelist_surf[setid][i];
		int f = facelist_surf[setid][i];
		surf_e[e][f] = load;
	}
}

void FixBoundary(int FixId, int nodeset, float fixT, int TableId){

	num_b   = (int*)realloc(num_b,  (FixId+1)*sizeof(int));
	fix_b   = (float*)realloc(fix_b,(FixId+1)*sizeof(float));
	table_b = (int*)realloc(table_b,(FixId+1)*sizeof(int));
	list_b  = (int**)realloc(list_b,(FixId+1)*sizeof(int*));

	fix_b[FixId] = fixT;
	num_b[FixId] = num_nset[nodeset];
	table_b[FixId] = TableId;

	list_b[FixId] = (int*)calloc(num_b[FixId],sizeof(int));	
	SearchSetList(nodeset,1);
	for(int i=0;i<num_b[FixId];i++){
		list_b[FixId][i] = list[i];	// 数组粘帖
	}

}

void EmptySfilm(){
	for (int i=0;i<elementnum_c;i++)
		for (int j=0;j<6;j++)
			surf_e[i][j] = 0;
}

void Read_step(FILE *in){
	char KeywordData[500] = "";
	char cmd[10][500], tmp[500];
	int end = 0, output = 1, num;

	while(!end){
		LineExtract(in,cmd,&num);

		if (cmd[0][1]=='*'){	// **注释
			for (int i=0;i<num;i++)	fprintf(out,"%s",cmd[i]);
		}
		else if (strstr(cmd[0],"*Step")){	// *Step
			fprintf(out,"%s",cmd[0]);
			for(int i=1;i<num;i++){
				if (strstr(cmd[i],"name")){	// 参数 name
					sscanf(cmd[i],"name=%s",&name_l);
					printf(" Running Step : %s\n",name_l);
				}else if (strstr(cmd[i],"inc")){	// 参数 name					
					// inc=xxxxx，设置最大步长，这里不管
				}else{
					output = 0;
					Warning(cmd[i]," is not a parameter in *Step");
				}if(output) fprintf(out,", %s",cmd[i]);
			}
		}
		else if(strstr(cmd[0],"*Heat Transfer")){// *Heat Transfer
			fprintf(out,"%s",cmd[0]);
			for (int i=1;i<num;i++)		fprintf(out,", %s",cmd[i]);
			fprintf(out,"\n");
			LineExtract(in,cmd,&num);	// 再读取一行数据
			sscanf(cmd[0],"%f",&timeinc_l);	sscanf(cmd[1],"%f",&totaltime_l);
			inc_l = (int)(totaltime_l/timeinc_l);
			fprintf(out,"%s, %s, %s, %s,",cmd[0],cmd[1],cmd[2],cmd[3]);
		}
		else if(strstr(cmd[0],"*Output")){		// *Output
			fprintf(out,"%s",cmd[0]);
			for (int i=1;i<num;i++)	fprintf(out," ,%s",cmd[i]);
		}
		else if(strstr(cmd[0],"*Sfilm")){		// *Sfilm
			fprintf(out,"%s",cmd[0]);
			for (int i=1;i<num;i++)	{
				if(strstr(cmd[i],"amplitude")){
					sscanf(cmd[i],"amplitude=%s",&tmp);
					fprintf(DEBUG," *Sfilm Amplitude name = %s\n",tmp);
					int k=SearchSetName(tmp,20);	sprintf(tmp,"*Sfilm=%d",k);
				}if(strstr(cmd[i],"op")){
					sscanf(cmd[i],"op=%s",&tmp);
					if(strstr(tmp,"NEW")){
						EmptySfilm();
					}
				}fprintf(out,", %s",cmd[i]);
			}strcpy(KeywordData,tmp);
		}
		else if(strstr(cmd[0],"*Boundary")){		// *Sfilm
			fprintf(out,"%s",cmd[0]);
			for (int i=1;i<num;i++){
				if(strstr(cmd[i],"amplitude")){
					sscanf(cmd[i],"amplitude=%s",&tmp);
					fprintf(DEBUG," *Boundary Amplitude name = %s\n",tmp);
					int k=SearchSetName(tmp,20);	sprintf(tmp,"*Boundary=%d",k);
				}if(strstr(cmd[i],"op")){
					sscanf(cmd[i],"op=%s",&tmp);
					if(strstr(tmp,"NEW")){
						// EmptySfilm();
					}
				}fprintf(out,", %s",cmd[i]);
			}strcpy(KeywordData,tmp);
		}
		else if(strstr(cmd[0],"*model change")){	// *model change
			sprintf(tmp,"%s, %s",cmd[0],cmd[1]);	// cmd[1] : remove/add
			fprintf(out,"%s",tmp);			
			strcpy(KeywordData,tmp);
		}else if(strstr(cmd[0],"*Restart")){
			fprintf(out,"%s",cmd[0]);
			for(int i=1;i<num;i++)	fprintf(out,", %s",cmd[i]);
		}
		else if(strstr(cmd[0],"*End Step")){		// *End Step
			fprintf(out,"%s",cmd[0]);	end = 1;
		}

		/* * * * * * * * * * * * * * * * 
		* *							 * *
		* *			DATA BLOCK		 * *
		* *							 * *
		* * * * * * * * * * * * * * * */

		else{	// 数据段
			if (!strcmp(KeywordData,"*model change, remove")){
				for(int i=0;i<num;i++){
					deactivate(cmd[i]);
					fprintf(out,"%s,",cmd[i]);
				}
			}else if (!strcmp(KeywordData,"*model change, add")){
				for(int i=0;i<num;i++){
					activate(cmd[i]);
					fprintf(out,"%s,",cmd[i]);
				}
			}else if (strstr(KeywordData,"*Sfilm")){
				int k,surfset;	float beta;	// k临时变量，存储气温表号
				sscanf(KeywordData,"*Sfilm=%d",&k);
				sscanf(cmd[0],"%s",&tmp);
				surfset = SearchSetName(tmp,3);	// 搜索单元面组
				sscanf(cmd[1],"%s",&tmp);
				if(!strcmp(tmp,"F")){
					SurfLoadId ++;
					sscanf(cmd[3],"%f",&beta);	// beta 值
					Sfilm(surfset,SurfLoadId-1,beta,k);
					fprintf(out,"%s, %s, %s, %s",cmd[0],cmd[1],cmd[2],cmd[3]);
				}else if(!strcmp(tmp,"FNU"))
					Warning(" Opiton FNU is not supported in this version.");
			}else if (strstr(KeywordData,"*Boundary")){
				int k,nodeset;	float fixT;	// k临时变量，存储气温表号
				boundarynum_c ++;
				nodeset = SearchSetName(cmd[0],1);	// 搜索结点组
				sscanf(cmd[1],"%d",&k);		// 获得约束自由度编号
				if(k==11){	// 温度自由度
					sscanf(cmd[3],"%f",&fixT);
					if (KeywordData[9]=='='){
						sscanf(KeywordData,"*Boundary=%d",&k);
						FixBoundary(boundarynum_c-1,nodeset,fixT,k);
					}else{
						FixBoundary(boundarynum_c-1,nodeset,fixT,-1);	// 没有变化表
					}
					fprintf(out,"%s, %s, %s, %s",cmd[0],cmd[1],cmd[2],cmd[3]);
					
				}else{
					Warning(" Only Freedom 11 (Temperature) is supported in this version.");
				}
			}else
				Warning("Keyword ",cmd[0]," should not be included in step definition");
		}
		fprintf(out,"\n");
	}
}