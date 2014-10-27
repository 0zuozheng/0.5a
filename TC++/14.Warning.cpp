#include <stdio.h>
extern FILE *log_check;

void Warning(char* tmp1,char *tmp2="",char* tmp3=""){
	fprintf(log_check,"\n!!!!  Warning:\n\t\t\t%s%s\n",tmp1,tmp2,tmp3);
}