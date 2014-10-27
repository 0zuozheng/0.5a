#include <stdio.h>
#include <fstream>
#include <iomanip>
#include "Data.h"
using namespace std;

//声明文件指针
extern FILE *DEBUG;
extern FILE *inp_node,*inp_element,*inp_material,*inp_loadcase,*inp_table,*inp_set,*inp_pipe,*inp_section,*inp_initial,*inp_surface;
extern FILE *log_check;
extern FILE *out;

//声明用到的函数

void Assign_Boundary(int,int);
void Free_Variable();
void welcome(); //欢迎信息

extern void Read_Result();

//声明Assign函数
extern int SearchSetName(char* setname,int type);
extern void SearchSetList(int id,int type,char* cmd="Default");
extern void AssignInitial_T(char*,float);
extern void AssignMaterial(int,char*);

// 声明关键字解析相关变量、函数
extern void Warning(char*,char* tmp2="",char* tmp3="");
extern void LineExtract(FILE *in,char (*cmd)[500],int* num);
extern void KeywordElement(char (*cmd)[500],int num,int id);
extern void KeywordSet(char (*cmd)[500],int num,int* id, int* id2);
extern void KeywordInitial(char (*cmd)[500],int num,int id);
extern void KeywordSection(char (*cmd)[500],int num,int id);
extern void KeywordMaterial(char (*cmd)[500],int num,int* id);
extern void KeywordSurface(char (*cmd)[500],int num,int* SetId,int* SurfId);
extern void KeywordTable(char (*cmd)[500],int num,int* id);