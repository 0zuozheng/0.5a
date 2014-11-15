#include "Tecplot.h"

void Tecplot(int ln,int pn,int ptn) //ln -> loadcase_now;  ptn -> period_totalnow; pn -> period_now
{
	if(ptn==0){
		tecplot<<"TITLE = \" Job Result \" "<<endl;
		tecplot<<"VARIABLES = \"X\",\"Y\",\"Z\", "<<endl;
		tecplot<<"\"Temperature\" "<<endl;
	}	tecplot<<endl;
	
	tecplot<<"ZONE"<<endl<<"DATAPACKING = BLOCK, ZONETYPE = FEBRICK, T = \"Step="<<ln<<" Inc="<<pn+1<<"\""<<endl;
	tecplot<<"STRANDID = 1,"<<"SOLUTIONTIME = "<<time_totalnow<<endl;
	tecplot<<"NODES= "<<nodenum_c<<endl;
	tecplot<<"ELEMENTS= "<<maxelement_l<<endl;
	if(ptn>0.01)	tecplot<<"VARSHARELIST = ([1,2,3]=1)"<<endl;
//	tecplot<<"VARLOCATION=([4]=CELLCENTERED)"<<endl;
	if(pn>0)	tecplot<<"CONNECTIVITYSHAREZONE = "<<ptn<<endl;
	tecplot<<"#### Increment:"<<ptn+1<<endl;
	
	// ZUOZUO  这里有待改进
	//tecplot<<"#### Set Definition"<<endl;
	//for(int i=0;i<maxelement_l[ln];i++){
	//	tecplot<<setw(8)<<set_e[i];
	//	if((i+1)%16==0)	tecplot<<endl;
	//}	tecplot<<endl;
	
	if(ptn==0){
		tecplot<<"#### X"<<endl;
		for(int i=0;i<nodenum_c;i++){	//Temperature Output
			tecplot<<setw(16)<<setprecision(3)<<xyz_n[i][0];
			if((i+1)%8==0)	tecplot<<endl;
		}	tecplot<<endl;
		tecplot<<"#### Y"<<endl;
		for(int i=0;i<nodenum_c;i++){	//Temperature Output
			tecplot<<setw(16)<<setprecision(3)<<xyz_n[i][1];
			if((i+1)%8==0)	tecplot<<endl;
		}	tecplot<<endl;
		tecplot<<"#### Z"<<endl;
		for(int i=0;i<nodenum_c;i++){	//Temperature Output
			tecplot<<setw(16)<<setprecision(3)<<xyz_n[i][2];
			if((i+1)%8==0)	tecplot<<endl;
		}	tecplot<<endl;
	}

	tecplot<<"#### Temperature"<<endl;
	for(int i=0;i<nodenum_c;i++){	//Temperature Output
		tecplot<<setw(16)<<setprecision(3)<<t_n[i];
		if((i+1)%8==0)	tecplot<<endl;
	}	tecplot<<endl;

	if(pn==0){
		tecplot<<"#### Connectivity"<<endl;
		for(int i=0;i<maxelement_l;i++){
			for(int j=0;j<8;j++){	tecplot<<setw(8)<<node_e[i][j]+1;}
			tecplot<<endl;
		}	tecplot<<endl;
	}
}