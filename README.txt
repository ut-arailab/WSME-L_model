"Accurate prediction of protein folding mechanisms by simple structure-based statistical mechanical models"

Koji Ooka(1,2) and Munehito Arai(1,2,3)
1 Department of Physics, Graduate School of Science, The University of Tokyo, 3-8-1 Komaba, Meguro, Tokyo 153-8902, Japan
2 Komaba Organization for Educational Excellence, College of Arts and Sciences, The University of Tokyo, 3-8-1 Komaba, Meguro, Tokyo 153-8902, Japan
3 Department of Life Sciences, Graduate School of Arts and Sciences, The University of Tokyo, 3-8-1 Komaba, Meguro, Tokyo 153-8902, Japan
E-mail: arai@bio.c.u-tokyo.ac.jp



# Calculation program for WSME-L model and WSME-L(SS intact) model

## 1. WSME-L model
This program was used to calculate two-dimensional free energy landscape of a model protein, src SH3 domain, based on the WSME-L model.
The source code and input data are stored in "src" and "data" directories, respectively.
Do not change the location of the files when running the calculation.

The source code is divided into five files: 
- main_srcsh3.cpp
- wsme_l_free_energy_2d.hpp
- wsme_l_free_energy_2d_1.cpp
- wsme_l_free_energy_2d_2.cpp

All parameters for calculating the free energy of src SH3 domain were already described in the main code.

The input data consist of two files:
- srcsh3_dm.dat (Cα-Cα distance map)
- srcsh3_cm.dat (contact map)
- srcsh3_ec.dat (entropic cost)


### Compilation
This program does not require any installation and can be run by compiling the source code.

Compilation was tested on gcc version 7.5.0 (OS: Ubuntu 18.0.4.6 LTS):

g++ main_srcsh3.cpp wsme_l_free_energy_2d_1.cpp wsme_l_free_energy_2d_2.cpp -std=c++17 -O3 -fopenmp


### Execution

./a.out

16GB RAM or more is recommended.
Expected run time is about 13 sec on a standard desktop computer.


### WSME-L(SS) model 
The two-dimensional free energy landscape of WSME-L(SS) model can be calculated with the same program as in the WSME-L model, using a contact map in which the contact energy of disulfide bonds is taken into account.


### Φ-value analysis
To obtain theoretical Φ-values, it is necessary to calculate a partition function with the perturbations due to pseudo-mutation. This can be calcualted with the same program as in the WSME-L model, using a contact map that takes into account perturbations due to pseudo-mutation.



## 2. WSME-L(SS intact) model
This program was used to calculate two-dimensional free energy landscape of a model protein, BPTI, based on the WSME-L model.
The source code and input data are stored in "src" and "data" directories, respectively.
Do not change the location of the files when running the calculation.

The source code is divided into five files: 
- main_bpti.cpp
- wsme_lss_free_energy_2d.hpp
- wsme_lss_free_energy_2d_1.cpp
- wsme_lss_free_energy_2d_2.cpp

All parameters for calculating the free energy of BPTI were already described in the main code.

The input data consist of two files:
- bpti_dm.dat (Cα-Cα distance map)
- bpti_cm.dat (contact map)
- bpti_ec.dat (entropic cost)


### Compilation
This program does not require any installation and can be run by compiling the source code.

Compilation was tested on gcc version 7.5.0 (OS: Ubuntu 18.0.4.6 LTS):

g++ main_bpti.cpp wsme_lss_free_energy_2d_1.cpp wsme_lss_free_energy_2d_2.cpp -std=c++17 -O3 -fopenmp


### Execution

./a.out

16GB RAM or more is recommended.
Expected run time is about 74 sec on a standard desktop computer.
