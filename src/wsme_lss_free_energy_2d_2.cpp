#include "wsme_lss_free_energy_2d.hpp"

int FreeEnergy2D::get_contact_number(const double cutoff){
	int cn=0;
	for(const auto&val1:contact_map){
    for(const auto&val2:val1){
      if(val2<cutoff)cn++;
    }
  }
	return cn;	
}

double FreeEnergy2D::get_absolute_max_contact_map(){
  double absmax=0;
  for(const auto&val1:contact_map){
    for(const auto&val2:val1){
      if(fabs(val2)>absmax)absmax=fabs(val2);
    }
  }
  return absmax;
}

void FreeEnergy2D::calculate_weight(const double energy,const double temperature){
  double effective_energy[n+1][n+1];
  std::fill(effective_energy[0],effective_energy[n+1],0);

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      effective_energy[i][j]=effective_energy[i][j-1];
      for(int k=i;k<=j;k++){
      	effective_energy[i][j]+=(energy*contact_map[k][j]/absmax_contact_map)/(boltzmann_const*temperature);
      }
    }
  }

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      for(int k=i;k<=j;k++){
      	effective_energy[i][j]-=entropic_cost[k]/boltzmann_const/1000;
      }
    }
  }

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      weight[i][j]=exp(-effective_energy[i][j]);
    }
  }
}

void FreeEnergy2D::calculate_weight_linker(const double energy,const double temperature,int uu,int vv){
  
	//Contact map without contacts around linker
	vector<vector<double>> contact_map_=contact_map;
	if(1<=uu  &&  uu<=n&&1<=vv  &&  vv<=n) contact_map_[uu ][vv ]=0;
	if(1<=uu-1&&uu-1<=n&&1<=vv  &&  vv<=n) contact_map_[uu-1][vv]=0;
	if(1<=uu+1&&uu+1<=n&&1<=vv  &&  vv<=n) contact_map_[uu+1][vv]=0;
	if(1<=uu  &&  uu<=n&&1<=vv-1&&vv-1<=n) contact_map_[uu][vv-1]=0;
	if(1<=uu  &&  uu<=n&&1<=vv+1&&vv+1<=n) contact_map_[uu][vv+1]=0;

	//Boltzmann weight
  double effective_energy[n+1][n+1];
  std::fill(effective_energy[0],effective_energy[n+1],0);

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      effective_energy[i][j]=effective_energy[i][j-1];
      for(int k=i;k<=j;k++){
      	effective_energy[i][j]+=(energy*contact_map_[k][j]/absmax_contact_map)/(boltzmann_const*temperature);
      }
    }
  }

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      for(int k=i;k<=j;k++){
      	effective_energy[i][j]-=entropic_cost[k]/boltzmann_const/1000;
      }
    }
  }

  for(int i=1;i<=n;i++){
    for(int j=i;j<=n;j++){
      weight[i][j]=exp(-effective_energy[i][j]);
    }
  }

	//Boltzmann weight of linker
  double *wl000=weight_linker[0][0][0];
  for(int k=3;k<=n;k++){
    for(int l=k;l<=n;l++){
      for(int i=1;i<=n-2;i++){
        wl000[i]=0;
        for(int j=k;j<=l;j++){
          wl000[i]+=(energy*contact_map_[i][j]/absmax_contact_map)/(boltzmann_const*temperature);
        }
      }
      double **wlkl=weight_linker[k][l];
      for(int i=1;i<=n-2;i++){
        for(int j=i;j<=n-2;j++){
          if(j==i){
            wlkl[i][j]=wl000[j];
          }else{
            wlkl[i][j]=wl000[j]+wlkl[i][j-1];
          }
        }
      }
    }
  } 
  for(int k=3;k<=n;k++){
    for(int l=k;l<=n;l++){
      for(int i=1;i<=n-2;i++){
        for(int j=i;j<=n-2;j++){
          weight_linker[k][l][i][j]=exp(-weight_linker[k][l][i][j]);
        }
      }
    }
  }
}

void FreeEnergy2D::allocate(double **&p){
  p=new double*[na+1];
  for(int i=0;i<na+1;i++){
    p[i]=new double[n2+1]();
  }
}

void FreeEnergy2D::allocate(double ***&p,const int size){
	p=new double**[size];
  for(int i=0;i<size;i++){
    p[i]=new double*[na+1];
    for(int j=0;j<na+1;j++){
      p[i][j]=new double[n2+1]();
    }
  }
}

void FreeEnergy2D::allocate_r(){
  r=new double***[n+2];
  for(int i=0;i<n+2;i++){
    r[i]=new double**[n+2];
    for(int j=i;j<n+2;j++){
      const int op1=(j==0)?0:order_parameter1[i+1][j-1];
      const int op2=(j==0)?0:order_parameter2[i+1][j-1];
      r[i][j]=new double*[op1+1];
      for(int k=0;k<=op1;k++){
        r[i][j][k]=new double[op2+1];
      }
    }
  }
}

void FreeEnergy2D::free(double **&p){
  for(int i=0;i<na+1;i++){
    delete[] p[i];
  }
  delete[] p;
}

void FreeEnergy2D::free(double ***&p,const int size){
  for(int i=0;i<size;i++){
    for(int j=0;j<na+1;j++){
      delete[] p[i][j];
    }
    delete[] p[i];
  }  
  delete[] p;
}

void FreeEnergy2D::free_r(){
  for(int i=0;i<n+2;i++){
    for(int j=i;j<n+2;j++){
      const int op1=(j==0)?0:order_parameter1[i+1][j-1];
      for(int k=0;k<=op1;k++){
        delete[] r[i][j][k];
      }
      delete[] r[i][j];
    }
    delete[] r[i];
  }
  delete[] r;
}


vector<vector<double>> FreeEnergy2D::to_vector(double **&p){
  vector<vector<double>> ret(na+1,vector<double>(n2+1));
	for(int i=0;i<=na;i++){
		for(int j=0;j<=n2;j++){
			ret[i][j]=p[i][j];
		}
	}
  return ret;
}

void FreeEnergy2D::calculate_z(double ***z,int o){
	if(o==-1){
		z[1][0][0]=1;
		return;
	}

  int op1,op2,op3,op4;
  double ww;

  for(int i=1;i<=o+2;i++){
    op1=order_parameter1[i][o];
    op2=order_parameter2[i][o];
    for(int j=0;j<=op1;j++){
			std::fill_n(z[i][j],op2+1,0);
    }
  }

  z[o+2][0][0]=1;
  for(int j=o;j>=0;j--){
    op1=order_parameter1[j+2][o];
    op2=order_parameter2[j+2][o];
    for(int i=1;i<=j+1;i++){
      op3=order_parameter1[i][j];
      op4=order_parameter2[i][j];
      ww=weight[i][j];
      for(int k=0;k<=op1;k++){
      	for(int l=0;l<=op2;l++){
      		z[i][k+op3][l+op4]+=z[j+2][k][l]*ww;
      	}
      }
    }
  }
}

void FreeEnergy2D::calculate_r(){
	int op1,op2;
  double ***rr;
  rr=new double**[n+3];
  for(int i=0;i<n+3;i++){
    op1=order_parameter1[i][n];
    op2=order_parameter2[i][n];
    rr[i]=new double*[op1+1];
    for(int j=0;j<op1+1;j++){
      rr[i][j]=new double[op2+1]();
    }
  }

  for(int j=0;j<=n+1;j++){
    calculate_z(rr,j-1);
		for(int i=0;i<=j;i++){
      op1=(j==0)?0:order_parameter1[i+1][j-1];
      op2=(j==0)?0:order_parameter2[i+1][j-1];
      for(int k=0;k<=op1;k++){
      	for(int l=0;l<=op2;l++){
       		r[i][j][k][l]=rr[i+1][k][l];
       	}
      }
    }
  }

  for(int i=0;i<n+3;i++){
    op1=order_parameter1[i][n];
    for(int j=0;j<op1+1;j++){
      delete[] rr[i][j];
    }
    delete[] rr[i];
  }
	delete[] rr;
}



void FreeEnergy2D::calculate_zl(int i1,int j1,double ***z,int o){
  if(o==-1){
		z[1][0][0]=1;
		return;
	}

  bool top1=(i1<=v1&&v1<=j1);

  if(!top1){
    int a=o+2;
    int op1=(a==1)?0:order_parameter1[1][a-2],op2=(a==1)?0:order_parameter2[1][a-2];
    auto rr=r[0][a-1];
    for(int i=0;i<=op1;i++){
      for(int j=0;j<=op2;j++){
        z[1][i][j]=rr[i][j];
      }
    }
    return;
  }

  int op1,op2,op3,op4;
  double ww;
  auto wl1=weight_linker[i1][j1];
  for(int i=1;i<=o+2;i++){
		op1=order_parameter1[i][o];
		op2=order_parameter2[i][o];
		for(int j=0;j<=op1;j++){
			std::fill_n(z[i][j],op2+1,0);
		}
	}

  z[o+2][0][0]=1;
  for(int j=o;j>=0;j--){
    op1=order_parameter1[j+2][o];
    op2=order_parameter2[j+2][o];
    for(int i=1;i<=j+1;i++){
      op3=order_parameter1[i][j];
      op4=order_parameter2[i][j];
      bool low1=(i<=u1&&u1<=j);
      ww=weight[i][j];
      if(top1&&low1) ww*=wl1[i][j];
			for(int k=0;k<=op1;k++){
      	for(int l=0;l<=op2;l++){
      		z[i][k+op3][l+op4]+=z[j+2][k][l]*ww;
      	}
      }
    }
  }
}

void FreeEnergy2D::calculate_zl(int i1,int j1,int i2,int j2,double ***z,int o){
  if(o==-1){
		z[1][0][0]=1;
		return;
	}

  bool top1=(i1<=v1&&v1<=j1);
  bool top2=(i2<=v2&&v2<=j2);

  if(!top1&&!top2){
    int a=o+2;
    int op1=(a==1)?0:order_parameter1[1][a-2],op2=(a==1)?0:order_parameter2[1][a-2];
    auto rr=r[0][a-1];
    for(int i=0;i<=op1;i++){
      for(int j=0;j<=op2;j++){
        z[1][i][j]=rr[i][j];
      }
    }
    return;
  }

  bool dp12=(i1==i2&&j1==j2)&&top1&&top2;
  int op1,op2,op3,op4;
  double ww;
  auto wl1=weight_linker[(i1!=n+1)?i1:0][j1],wl2=weight_linker[i2][j2];

  for(int i=1;i<=o+2;i++){
		op1=order_parameter1[i][o];
		op2=order_parameter2[i][o];
		for(int j=0;j<=op1;j++){
			std::fill_n(z[i][j],op2+1,0);
		}
	}

  z[o+2][0][0]=1;
  for(int j=o;j>=0;j--){
    op1=order_parameter1[j+2][o];
    op2=order_parameter2[j+2][o];
    for(int i=1;i<=j+1;i++){
      op3=order_parameter1[i][j];
      op4=order_parameter2[i][j];
      bool low1=(i<=u1&&u1<=j);
      bool low2=(i<=u2&&u2<=j);
      ww=weight[i][j];
      if(dp12){
        if(low1||low2) ww*=wl1[i][j];
      }else{
        if(top1&&low1) ww*=wl1[i][j];
        if(top2&&low2) ww*=wl2[i][j];
      }
			for(int k=0;k<=op1;k++){
      	for(int l=0;l<=op2;l++){
      		z[i][k+op3][l+op4]+=z[j+2][k][l]*ww;
      	}
      }
    }
  }

}

void FreeEnergy2D::calculate_zl(int i1,int j1,int i2,int j2,int i3,int j3,double ***z,int o){
  if(o==-1){
		z[1][0][0]=1;
		return;
	}

  bool top3=(i3<=v3&&v3<=j3);

  if(!top3){
    int a=o+2;
    int op1=(a==1)?0:order_parameter1[1][a-2],op2=(a==1)?0:order_parameter2[1][a-2];
    auto rr=r[0][a-1];
    for(int i=0;i<=op1;i++){
      for(int j=0;j<=op2;j++){
        z[1][i][j]=rr[i][j];
      }
    }
    return;
  }

  int op1,op2,op3,op4;
  double ww;
  auto wl3=weight_linker[i3][j3];
  for(int i=1;i<=o+2;i++){
		op1=order_parameter1[i][o];
		op2=order_parameter2[i][o];
		for(int j=0;j<=op1;j++){
			std::fill_n(z[i][j],op2+1,0);
		}
	}

  z[o+2][0][0]=1;
  for(int j=o;j>=0;j--){
    op1=order_parameter1[j+2][o];
    op2=order_parameter2[j+2][o];
    for(int i=1;i<=j+1;i++){
      op3=order_parameter1[i][j];
      op4=order_parameter2[i][j];
      bool low3=(i<=u3&&u3<=j);
      ww=weight[i][j];
      if(top3&&low3) ww*=wl3[i][j];
			for(int k=0;k<=op1;k++){
      	for(int l=0;l<=op2;l++){
      		z[i][k+op3][l+op4]+=z[j+2][k][l]*ww;
      	}
      }
    }
  }
}

void FreeEnergy2D::calculate_r_2(double **r1,int i1,int j1,double ***z){
  double **srw;
  allocate(srw);

  clear(r1,i1-2);
  for(int j2=v2;j2<=i1-2;j2++){

    clear(srw,j2);
    for(int i2=1;i2<=v2+1;i2++){
      
    	calculate_zl(i1,j1,i2,j2,z,i2-2);

    	sum_w(srw,z[1],weight[i2][j2]*(i1<=v1&&v1<=j1&&i2<=u1&&u1<=j2?weight_linker[i1][j1][i2][j2]:1),i2,j2);
    }
    sum_r(r1,srw,j2+1,i1-1);
  }
  free(srw);
}

void FreeEnergy2D::calculate_rr_2(double **r1,int i1,int j1,double ***z){
  double **srw,**r2;
  allocate(srw);
  allocate(r2);

  clear(r1,i1-2);
  for(int j2=v2;j2<=i1-2;j2++){

    clear(srw,j2);
    for(int i2=1;i2<=v3+1;i2++){
    	calculate_zl(i1,j1,i2,j2,i2,j2,z,i2-2);
    	sum_w(srw,z[1],weight[i2][j2]*(i1<=v1&&v1<=j1&&i2<=v2&&v2<=j2?weight_linker[i1][j1][i2][j2]:1),i2,j2);
    }
    for(int i2=v3+2;i2<=v2+1;i2++){
      calculate_rr_3(r2,i1,j1,i2,j2,z);
    	sum_w(srw,r2,weight[i2][j2]*(i1<=v1&&v1<=j1&&i2<=v2&&v2<=j2?weight_linker[i1][j1][i2][j2]:1),i2,j2);
    }
    sum_r(r1,srw,j2+1,i1-1);
  }
  free(r2);
  free(srw);
}

void FreeEnergy2D::calculate_rr_3(double **r2,int i1,int j1,int i2,int j2,double ***zp){
  double **srw;
  allocate(srw);

  clear(r2,i2-2);
  for(int j3=v3;j3<=i2-2;j3++){
    clear(srw,j3);
    for(int i3=1;i3<=v3+1;i3++){
    	calculate_zl(i1,j1,i2,j2,i3,j3,zp,i3-2);
    	sum_w(srw,zp[1],weight[i3][j3],i3,j3);
    }
    sum_r(r2,srw,j3+1,i2-1);
  }
  free(srw);
}

void FreeEnergy2D::clear(double **zz,int a){
  int op1=order_parameter1[1][a],op2=order_parameter2[1][a];
  for(int i=0;i<=op1;i++){
		std::fill_n(zz[i],op2+1,0);
  }
}

void FreeEnergy2D::sum_w(double **rw,double **rr,double ww,int a,int b){
  int op1=(a==1)?0:order_parameter1[1][a-2],op2=(a==1)?0:order_parameter2[1][a-2];
	int op3=order_parameter1[a][b],op4=order_parameter2[a][b];
  for(int i=0;i<=op1;i++){
    for(int j=0;j<=op2;j++){
      rw[i+op3][j+op4]+=rr[i][j]*ww;
    }
  }
}

void FreeEnergy2D::sum_w_p(double **rw,double **rr,double ww,int a,int b){
  int op1=(a==1)?0:order_parameter1[1][a-2],op2=(a==1)?0:order_parameter2[1][a-2];
	int op3=order_parameter1[a][b],op4=order_parameter2[a][b];
  for(int i=0;i<=op1;i++){
    for(int j=0;j<=op2;j++){
      #pragma omp atomic
      rw[i+op3][j+op4]+=rr[i][j]*ww;
    }
  }
}

void FreeEnergy2D::sum_r(double **rr,double **srw,int a,int b){
  int op1=order_parameter1[1][a-1],op2=order_parameter2[1][a-1];
	int m3=order_parameter1[a+1][b-1],m4=order_parameter2[a+1][b-1];
  auto rrr=r[a][b];
  for(int i=0;i<=op1;i++){
    for(int j=0;j<=op2;j++){
      for(int k=0;k<=m3;k++){
        for(int l=0;l<=m4;l++){
          rr[i+k][j+l]+=srw[i][j]*rrr[k][l];
        }
      }
    }
  }
}

double FreeEnergy2D::calculate_wld(const double energy,const double temperature,const int uu,const int vv){
	double el=0;
  if(1<=uu  &&  uu<=n&&1<=vv  &&  vv<=n) el+=contact_map[uu ][vv ];
	if(1<=uu-1&&uu-1<=n&&1<=vv  &&  vv<=n) el+=contact_map[uu-1][vv];
	if(1<=uu+1&&uu+1<=n&&1<=vv  &&  vv<=n) el+=contact_map[uu+1][vv];
	if(1<=uu  &&  uu<=n&&1<=vv-1&&vv-1<=n) el+=contact_map[uu][vv-1];
	if(1<=uu  &&  uu<=n&&1<=vv+1&&vv+1<=n) el+=contact_map[uu][vv+1];
  return std::exp(-(energy*el/absmax_contact_map)/(boltzmann_const*temperature));
}

double FreeEnergy2D::calculate_sr(const int i,const double coef){
	return (i>0)?-3/2.0*(std::log(i)+coef/i):0;
}

double FreeEnergy2D::combination_fraction(const int i,const int j,const int m,const int n){
  if(m<j||n-m<i-j||i<j) return 0;
  
  double cf=1;
  for(int k=1;k<=i;k++){
    double c0=(n-i+k)*1.0/k;
    double c1=(k<=j)?(m-j+k)*1.0/k:1;
    double c2=(k<=i-j)?(n-m-i+j+k)*1.0/k:1;
    cf*=c1*c2/c0;
  }
  return cf;
}