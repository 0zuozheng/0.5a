#include "Construct.h"

extern float i3[3],w3[3];
float face[6]={1,-1,1,-1,1,-1};


void Tface(int ln,int en,int facenum,int LoadId) //ln -> loadcase_now; en -> element_now
{
	int j,k,l,r,k1,k2,k3,m;
	float xyz[3],st[3][2],*sf_,**dsf_,intergral_x,intergral_y,intergral_z,tmp;
	m=material_e[en];
	GetEleXYZ(8,en);

	sf_ = new float[8];		 Alloc2DArray_float(&dsf_,3,8);

	if (SurfLoadType[LoadId]==3){	// 第三类边界条件
		k1=facenum/2;		k2=ijk[k1];		k3=ijk[k2];
		xyz[k1]=face[facenum];
		for(j=0;j<3;j++)
		{
			xyz[k2]=i3[j];
			for(k=0;k<3;k++)
			{
				xyz[k3]=i3[k];
				GetSf(xyz[0],xyz[1],xyz[2],sf_);	// 获得 (偏N/偏epsilon)
				GetDsf(xyz[0],xyz[1],xyz[2],dsf_);	// 获得 (偏N/偏epsilon)
				for(l=0;l<3;l++)
					for(r=0;r<2;r++) st[l][r]=0;
				for(l=0;l<8;l++)
				{	// (偏x/偏epsilon)=(偏N/偏epsilon)*x0+(偏N/偏epsilon)*x1+...(偏N/偏epsilon)*x7
					st[0][0]+=dsf_[k2][l]*nodeXYZ_e[l][0];  st[0][1]+=dsf_[k3][l]*nodeXYZ_e[l][0];
					st[1][0]+=dsf_[k2][l]*nodeXYZ_e[l][1];  st[1][1]+=dsf_[k3][l]*nodeXYZ_e[l][1];
					st[2][0]+=dsf_[k2][l]*nodeXYZ_e[l][2];  st[2][1]+=dsf_[k3][l]*nodeXYZ_e[l][2];
				}
				intergral_x=st[1][0]*st[2][1]-st[1][1]*st[2][0];
				intergral_y=st[2][0]*st[0][1]-st[0][0]*st[2][1];
				intergral_z=st[0][0]*st[1][1]-st[0][1]*st[1][0];
				tmp=sqrt(intergral_x*intergral_x+intergral_y*intergral_y+intergral_z*intergral_z);	// 见《有限单元法》P134

				fprintf(log_N,"\n积分点位置 x=%f, y=%f, z=%f, A=%f\n",xyz[0],xyz[1],xyz[2],tmp);

				for(l=0;l<8;l++)
				{
					//Q3[en][l]+=sf_[l]*tmp*w3[j]*w3[k]*beta_m[m][kind_e[en][i]]*diffusivity_m[m][k1]/conduction_m[m][k1];//i - face
					Q3[en][l]+=sf_[l]*tmp*w3[j]*w3[k]*diffusivity_m[m][k1]/conduction_m[m][k1]; // 插值函数 N * 面积微元 / (c * rou)
					for(r=0;r<8;r++)
						H[en][l][r]+=sf_[r]*sf_[l]*tmp*w3[j]*w3[k]*BetaValue[LoadId]*diffusivity_m[m][k1]/conduction_m[m][k1];
				}
			}
		}
	}
}
/*  形状函数	*/