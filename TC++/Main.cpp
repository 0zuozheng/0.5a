#include <stdio.h>
extern void File();
extern void ReadData();
extern void ElementAnalysis();
extern void Assembly();	

void main()
{
	File();
	ReadData();
	ElementAnalysis();
	Assembly();
	getchar();
}