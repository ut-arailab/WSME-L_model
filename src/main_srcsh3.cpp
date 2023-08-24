#include "wsme_l_free_energy_2d.hpp"

int main(){
  const int residue_number=60;
  const int residue_number_domain1=32;
  const int residue_number_domain2=28;

  const string fn_distance_map="../data/srcsh3_dm.dat";
  const string fn_contact_map="../data/srcsh3_cm.dat";
  const string fn_entropic_cost="../data/srcsh3_ec.dat";

  const double energy=2.5064236176;
  const double temperature=298;
  const double ring_entropy_scaling=0.5;
  const double cutoff=-0.6;
  const int openmp_number=128;
  const string fn_free_energy="../output/srcsh3_fe2d_001.dat";


  FreeEnergy2D protein(residue_number,residue_number_domain1,residue_number_domain2);
  protein.import_distance_map(fn_distance_map);
  protein.import_contact_map(fn_contact_map,cutoff);
  protein.import_entropic_cost(fn_entropic_cost);
  protein.set_order_parameter();

  auto z0=protein.calculate_partition_function(energy,temperature);
  vector<vector<double>> partition_function=z0;
  for(int i=1;i<=residue_number;i++){
    for(int j=i+2;j<=residue_number;j++){
      if(protein.get_contact_map(i,j)<cutoff){
        protein.print_progress();
        auto zl=protein.calculate_partition_function_linker(energy,temperature,i,j,openmp_number);
        auto ws=protein.calculate_weight_ring_entropy(ring_entropy_scaling,i,j);
        for(unsigned int k=0;k<z0.size();k++){
          for(unsigned int l=0;l<z0[0].size();l++){
            partition_function[k][l]+=(zl[k][l]-z0[k][l])*ws[k][l];
          }
        }
      }
    }
  }
  protein.export_free_energy(fn_free_energy,partition_function);
}