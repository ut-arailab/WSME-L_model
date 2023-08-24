#include "wsme_lss_free_energy_2d.hpp"

FreeEnergy2D::FreeEnergy2D(int residue_number,int residue_number_domain1,int residue_number_domain2,int residue_number_domain3):
n(residue_number),n1(residue_number_domain1),n2(residue_number_domain2),n3(residue_number_domain3),na(n1+n3)
{
  distance_map.resize(n+1,vector<double>(n+1,0));
  contact_map.resize(n+1,vector<double>(n+1,0));
  entropic_cost.resize(n+1,0);

  order_parameter1=new int*[n+3];
  order_parameter2=new int*[n+3];
  for(int i=0;i<n+3;i++){
    order_parameter1[i]=new int[n+1]();
    order_parameter2[i]=new int[n+1]();
  }

  weight=new double*[n+2];
  for(int i=0;i<n+2;i++){
    weight[i]=new double[n+1];
    std::fill_n(weight[i],n+1,1);
  }

  weight_linker=new double***[n+1];
  for(int k=0;k<n+1;k++){
    weight_linker[k]=new double**[n+1];
    for(int l=0;l<n+1;l++){
      weight_linker[k][l]=new double*[n-1];
      for(int i=0;i<n-1;i++){
        weight_linker[k][l][i]=new double[n-1]();
      }
    }
  }
}

FreeEnergy2D::~FreeEnergy2D(){
  for(int i=0;i<n+3;i++){
  	delete[] order_parameter1[i];
  	delete[] order_parameter2[i];
  }
  delete[] order_parameter1;
  delete[] order_parameter2;

  for(int i=0;i<n+2;i++){
  	delete[] weight[i];
  }
  delete[] weight;

  for(int k=0;k<n+1;k++){
    for(int l=0;l<n+1;l++){
      for(int i=0;i<n-1;i++){
        delete[] weight_linker[k][l][i];
      }
      delete[] weight_linker[k][l];
    }
    delete[] weight_linker[k];
  }
  delete[] weight_linker;
}

void FreeEnergy2D::import_distance_map(const string fn_distance_map){
  std::ifstream in(fn_distance_map);
  if(!in){
    std::cout<<"Unable to open file: "+fn_distance_map<<std::endl;
    std::exit(1);
  }
  for(int i=1;i<=n;i++){
    for(int j=1;j<=n;j++){
      in>>distance_map[i][j];
    }
  }

  /*for(const auto&val1:distance_map){
    for(const auto&val2:val1){
      std::cout<<val2<<" ";
    }
    std::cout<<std::endl;
  }*/
}


void FreeEnergy2D::import_contact_map(const string fn_contact_map,const double cutoff,const vector<vector<int>> ss_residue){
  std::ifstream in(fn_contact_map);
  if(!in){
    std::cout<<"Unable to open file: "+fn_contact_map<<std::endl;
    std::exit(1);
  }
  for(int i=1;i<=n;i++){
    for(int j=1;j<=n;j++){
      in>>contact_map[i][j];
    }
  }

  for(const auto&val:ss_residue){
    int uu=val[0],vv=val[1];
 	  if(1<=uu  &&  uu<=n&&1<=vv  &&  vv<=n) contact_map[uu ][vv ]=0;
	  if(1<=uu-1&&uu-1<=n&&1<=vv  &&  vv<=n) contact_map[uu-1][vv]=0;
	  if(1<=uu+1&&uu+1<=n&&1<=vv  &&  vv<=n) contact_map[uu+1][vv]=0;
	  if(1<=uu  &&  uu<=n&&1<=vv-1&&vv-1<=n) contact_map[uu][vv-1]=0;
	  if(1<=uu  &&  uu<=n&&1<=vv+1&&vv+1<=n) contact_map[uu][vv+1]=0;
  }

  absmax_contact_map=get_absolute_max_contact_map();
  contact_number=ss_residue.size()*get_contact_number(cutoff);

  /*for(const auto&val1:contact_map){
    for(const auto&val2:val1){
      std::cout<<val2<<" ";
    }
    std::cout<<std::endl;
  }*/
}

void FreeEnergy2D::import_entropic_cost(const string fn_entropic_cost){
  std::ifstream in(fn_entropic_cost);
  if(!in){
    std::cout<<"Unable to open file: "+fn_entropic_cost<<std::endl;
    std::exit(1);
  }
  for(int i=1;i<=n;i++){
    in>>entropic_cost[i];
  }

  /*for(const auto&val:entropic_cost){
    std::cout<<val<<std::endl;
  }*/
}

void FreeEnergy2D::set_order_parameter(){
	for(int i=1;i<=n+2;i++){
		for(int j=i;j<=n;j++){
			if(j<=n1){
				order_parameter1[i][j]=j-i+1;
				order_parameter2[i][j]=0;
			}else if(j<=n1+n2){
				if(i<=n1){
					order_parameter1[i][j]=n1-i+1;
					order_parameter2[i][j]=j-n1;
				}else{
					order_parameter1[i][j]=0;
					order_parameter2[i][j]=j-i+1;
				}
			}else{
				if(i<=n1){
					order_parameter1[i][j]=j-i-n2+1;
					order_parameter2[i][j]=n2;
				}else if(i<=n1+n2){
					order_parameter1[i][j]=j-n1-n2;
					order_parameter2[i][j]=n1+n2-i+1;
				}else{
					order_parameter1[i][j]=j-i+1;
					order_parameter2[i][j]=0;
				}
			}
		}
	}
}

vector<vector<double>> FreeEnergy2D::initialize_partition_function(){
  vector<vector<double>> ret(na+1,vector<double>(n2+1,0));
  return ret;
}


double FreeEnergy2D::get_contact_map(const int i,const int j){
  return contact_map[i][j];
}

void FreeEnergy2D::print_progress(){
	static int i=0;
	i++;
	std::cout<<"\r";
  std::cout<<"Progress: "<<std::setfill(' ')<<std::setw(3)<<std::right<<i<<"/"<<contact_number<<std::flush;
  if(i==contact_number)std::cout<<std::endl;
}

vector<vector<double>> FreeEnergy2D::calculate_partition_function(const double energy,const double temperature){
  calculate_weight(energy,temperature);
	double ***z;
	allocate(z,n+3);
  calculate_z(z,n);  
	auto ret=to_vector(z[1]);
  free(z,n+3);
	return ret;
}

vector<vector<double>> FreeEnergy2D::calculate_partition_function_linker(const double energy,const double temperature,const int u1_,const int v1_,const int openmp_number){
  calculate_weight_linker(energy,temperature,u1_,v1_);

  u1=u1_;
  v1=v1_;

  allocate_r();
  calculate_r();

	double **r0,***srwp,***z;
  allocate(r0);
	allocate(srwp,n-v1+1);

  #pragma omp parallel for private(z) schedule(dynamic) num_threads(openmp_number)
	for(int ij=0;ij<(n-v1+1)*(v1+1);ij++){
    int j1=ij/(v1+1)+v1;
    int i1=ij%(v1+1)+1;
		allocate(z,n+3);
	  calculate_zl(i1,j1,z,i1-2);
    sum_w_p(srwp[j1-v1],z[1],weight[i1][j1],i1,j1);
    free(z,n+3);
	}
  for(int j1=v1;j1<=n;j1++){
    sum_r(r0,srwp[j1-v1],j1+1,n+1);
  }
  
	auto ret=to_vector(r0);
  free(srwp,n-v1+1);
  free(r0);
  free_r();
	return ret;
}

vector<vector<double>> FreeEnergy2D::calculate_partition_function_linker(const double energy,const double temperature,int u1_,int v1_,int u2_,int v2_,const int openmp_number){
  calculate_weight_linker(energy,temperature,u1_,v1_);
  double wld=calculate_wld(energy,temperature,u1_,v1_);
  
  if(v2_<=v1_){
    u1=u1_;
    v1=v1_;
    u2=u2_;
    v2=v2_;
  }else{
    u1=u2_;
    v1=v2_;
    u2=u1_;
    v2=v1_;  
  }

  allocate_r();
  calculate_r();

  double **r0,***srwp,**r1,***z;
  allocate(r0);
  allocate(srwp,n-v1+1);

  if(u1<=v2){
    #pragma omp parallel for private(r1,z) schedule(dynamic) num_threads(openmp_number)
    for(int ij=0;ij<(n-v1+1)*(v1+1);ij++){
      int j1=ij/(v1+1)+v1;
      int i1=ij%(v1+1)+1;
      
      allocate(r1);
      allocate(z,n+3);
      if(i1<=v2+1){
        calculate_zl(i1,j1,i1,j1,z,i1-2);
        sum_w_p(srwp[j1-v1],z[1],weight[i1][j1],i1,j1);
      }
      else if(i1<=v1+1){
        calculate_r_2(r1,i1,j1,z);
        sum_w_p(srwp[j1-v1],r1,weight[i1][j1],i1,j1);
      }
      free(z,n+3);
      free(r1);
    }
    for(int j1=v1;j1<=n;j1++){
      sum_r(r0,srwp[j1-v1],j1+1,n+1);
    }
  }
  else{
    v3=v2;
    u3=u2;
    v2=u1;
    u2=0;
    v1=v1;
    u1=0;

    #pragma omp parallel for private(r1,z) schedule(dynamic) num_threads(openmp_number)
    for(int ij=0;ij<(n-v1+1)*(v1+1);ij++){
      int j1=ij/(v1+1)+v1;
      int i1=ij%(v1+1)+1;

      allocate(r1);
      allocate(z,n+3);

      if(i1<=v3+1){
        calculate_zl(i1,j1,i1,j1,i1,j1,z,i1-2);
        sum_w_p(srwp[j1-v1],z[1],weight[i1][j1],i1,j1);
      } 
      else if(i1<=v2+1){
        calculate_rr_3(r1,i1,j1,i1,j1,z);
        sum_w_p(srwp[j1-v1],r1,weight[i1][j1],i1,j1);
      }
      else if(i1<=v1+1){
        calculate_rr_2(r1,i1,j1,z);
        sum_w_p(srwp[j1-v1],r1,weight[i1][j1],i1,j1);
      }
      free(z,n+3);
      free(r1);
    }

    for(int j1=v1;j1<=n;j1++){
      sum_r(r0,srwp[j1-v1],j1+1,n+1);
    }
    v1=v1;
    u1=v2;
    u2=u3;
    v2=v3;
    u3=0;
    v3=0;
  }

  auto ret=to_vector(r0);
  free(r0);
  free(srwp,n-v1+1);
  free_r();

  for(auto&val1:ret){
    for(auto&val2:val1){
      val2*=wld;
    }
  }
	return ret;
}

vector<vector<double>> FreeEnergy2D::calculate_weight_ring_entropy(const double ring_entropy_scaling,const int uu,const int vv){
 int nl=vv-uu+1;
 double coefficient=(std::pow(distance_map[uu][vv],2)-std::pow(3.8,2))/(2*20*3.8);
 vector<vector<double>> ret(na+1,vector<double>(n2+1,0));
 vector<double>sr(n+1,0);

  for(int i=0;i<=n;i++){
  double tmp=0;
    for(int j=0;j<=i;j++){
      tmp+=calculate_sr(nl-j,coefficient)*combination_fraction(i,j,nl,n);
    }
    sr[i]=tmp;
  }

	for(int i=0;i<=na;i++){
		for(int j=0;j<=n2;j++){
			ret[i][j]+=exp(ring_entropy_scaling*sr[i+j]);
		}
	}
  
  /*for(const auto&val1:ret){
    for(const auto&val2:val1){
      std::cout<<std::log(val2)<<" ";
    }
    std::cout<<std::endl;
  }*/

  return ret;
}

void FreeEnergy2D::export_free_energy(const string fn_free_energy,const vector<vector<double>>& partition_function){
	std::ofstream out(fn_free_energy);
  out<<std::scientific<<std::setprecision(6);
  for(const auto&val1:partition_function){
    for(const auto&val2:val1){
      out<<-log(val2)<<" ";
    }
    out<<std::endl;
  }
}
