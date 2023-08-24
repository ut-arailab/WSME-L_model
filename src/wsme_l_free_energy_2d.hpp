#pragma once

#include <algorithm>
#include<cmath>
#include<fstream>
#include<iomanip>
#include<iostream>
#include<string>
#include<vector>
using string=std::string;
template<typename T> using vector=std::vector<T>;

class FreeEnergy2D{
public:
  FreeEnergy2D(int,int,int,int=0);
  ~FreeEnergy2D();

  void import_distance_map(const string);
  void import_contact_map(const string,const double);
  void import_entropic_cost(const string);
  void set_order_parameter();
  double get_contact_map(const int,const int);
  void print_progress();
  vector<vector<double>> calculate_partition_function(const double,const double);
  vector<vector<double>> calculate_partition_function_linker(const double,const double, const int,const int,const int);
  vector<vector<double>> calculate_weight_ring_entropy(const double,const int,const int);
  void export_free_energy(const string,const vector<vector<double>>&);

private:
  FreeEnergy2D(const FreeEnergy2D&);
  FreeEnergy2D& operator=(const FreeEnergy2D&);

  const double boltzmann_const=1.987/1000;
  const int n;   //residue_number
  const int n1;  //residue_number_domain1
  const int n2;  //residue_number_domain2
  const int n3;  //residue_number_domain3
  const int na;  //residue_number_domain1+3

  vector<vector<double>> distance_map;
  vector<vector<double>> contact_map;
  vector<double> entropic_cost;
  double absmax_contact_map;
  int contact_number;
  int **order_parameter1,**order_parameter2;
  double **weight;
  double ****weight_linker;
  double ****r;//partition function fragment

  double get_absolute_max_contact_map();
  int get_contact_number(const double);
  void calculate_weight(const double,const double);
  void calculate_weight_linker(const double,const double,int,int);

  void allocate(double**&);
  void allocate(double***&,const int);
  void allocate_r();
  void free(double**&);
  void free(double***&,const int);
  void free_r();

  vector<vector<double>> to_vector(double**&);
  void calculate_z(double***,int);
  void calculate_r();
  void calculate_zl(int,int,double***,int,int,int);
  void sum_w_p(double**,double**,double,int,int);
  void sum_r(double**,double**,int,int);
  double calculate_wld(const double,const double, const int,const int);
  double calculate_sr(const int,const double);
  double combination_fraction(const int,const int,const int,const int);
};
