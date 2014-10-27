#include "CONSTRUCT.h"
#include <iomanip>
#include <cmath>
#include "imsl.h"

#pragma comment(lib, "imslcmath_imsl_dll.lib")
#pragma comment(linker, "\"/manifestdependency:type='Win32' name='Microsoft.VC80.CRT' version='8.0.50608.0' processorArchitecture='X86' publicKeyToken='1fc8b3b9a1e18e3b' language='*'\"")

//声明文件指针
extern FILE *log_solver,*log_answer;

//声明计算函数
void Insert(int row, int col, float value);
void SolveWrite(int ln,int pn);

extern void Tecplot(int,int,int);

extern void third_film(float*,int,int,int,int,int,int);	//第三类边界条件

extern void Rst_wrt_node(int,float*);

extern int rst_flag;		//标识，是否写出结点结果