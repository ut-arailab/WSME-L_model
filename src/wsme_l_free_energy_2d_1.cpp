#include "wsme_l_free_energy_2d.hpp"

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

  /*for(const auto&v1:distance_map){
    for(const auto&v2:v1){
      std::cout<<v2<<" ";
    }
    std::cout<<std::endl;
  }*/
}


void FreeEnergy2D::import_contact_map(const string fn_contact_map,const double cutoff){
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
  absmax_contact_map=get_absolute_max_contact_map();
  contact_number=get_contact_number(cutoff);

  /*for(const auto&v1:contact_map){
    for(const auto&v2:v1){
      std::cout<<v2<<" ";
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

  /*for(const auto&v:entropic_cost){
    std::cout<<v<<std::endl;
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

vector<vector<double>> FreeEnergy2D::calculate_partition_function_linker(const double energy,const double temperature,const int u1,const int v1,const int openmp_number){
  calculate_weight_linker(energy,temperature,u1,v1);
  double wld=calculate_wld(energy,temperature,u1,v1);

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
	  calculate_zl(i1,j1,z,i1-2,u1,v1);
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

  for(auto&val1:ret){
    for(auto&val2:val1){
      val2*=wld;
    }
  }
	return ret;
}

vector<vector<double>> FreeEnergy2D::calculate_weight_ring_entropy(const double ring_entropy_scaling,const int u1,const int v1){
 int nl=v1-u1+1;
 double coefficient=(std::pow(distance_map[u1][v1],2)-std::pow(3.8,2))/(2*20*3.8);
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
  
  /*for(const auto&v1:ret){
    for(const auto&v2:v1){
      std::cout<<std::log(v2)<<" ";
    }
    std::cout<<std::endl;
  }*/

  return ret;
}

void FreeEnergy2D::export_free_energy(const string fn_free_energy,const vector<vector<double>>& partition_function){
	std::ofstream out(fn_free_energy);
  out<<std::scientific<<std::setprecision(6);
  for(const auto&v1:partition_function){
    for(const auto&v2:v1){
      out<<-log(v2)<<" ";
    }
    out<<std::endl;
  }
}
