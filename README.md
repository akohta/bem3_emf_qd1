# bem3_emf_qd1
This is the three-dimensional electromagnetic field analysis program for one-dimensional periodic arrangement objects irradiated by a plane wave. 
This is based on boundary element method, the own developed numerical solution is used. 
This is the full vector field three-dimensional analysis, the corner problem free.
Intel Math Kernel Library is required. Gmsh is used for create a mesh data of object. 
The calculation program of quasi-periodic Green's function "d3_qpgf_d1" is used.  

![analysis model](model_qpbc1.png "analysis model (model_qpbc1.png)")  


## Usage of example code  
1. type 'make' command to compile.  
   The executable d3qd1_bv_solver, example1.out, example2.out are created. 
   The d3qd1_bv_solver is the main solver of boundary integral equations.
   The example1.out is the executable of souce code example1.c, it shows a simplest example using "bem3_emf_qd1". 
   The example2.out is the executable of souce code example2.c, it shows a example of electromagnetic field intensity analysis.  
   
2. type './d3qd1_bv_solver' with arguments of plane wave datafile name, periodicity datafile name, medium datafile name, mesh datafile name, output datafile name, rotation and translation settings (optional).  
   For example, './d3qd1_bv_solver ipw.txt periodicity_data.txt medium_data.txt sphere_m2.msh ex.dat'. 
   The ipw.txt is the sample of incident field datafile, a plane wave is defined in it. 
   The periodicity_data.txt is the sample of periodicity datafile, periodic boundary condition and lattice constant are defined in it. 
   The medium_data.txt is the sample of medium datafile, two mediums are defined in it. The domain number is assinged to the medium from 1 in order. 
   The sphere_m2.msh is the sample of mesh datafile, it is a two layered sphere object. 
   It was created by Gmsh geometry file sphere_m2.geo in the mesh_sample folder.
   The sphere_m2_image.png is the visualization result of the sphere_m2.msh. 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this case, the file ex.particles is output and the visualization result is ex_particle.png (using ParaView).  
   
3. type './example1.out' with an argument of datafile name output by d3qd1_bv_solver.  
   For example, './example1.out ex.dat'. This executable calculates electromagnetic field, radiaton force and torque.  

4. type './example2.out' with an artument of datafile name output by d3qd1_bv_solver.  
   For example, './example2.out ex.dat'. This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of intensity distributions, created by Gnuplot script gscritp_example2.plt
   (using ImageMagick to convert eps to png).  

Please see d3qd1_src/bem3_emf_qd1.h for detail of functions.
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS.  

![mesh image 0](sphere_m2_image.png "mesh image of the object (sphere_m2_image.png)") 
![point cloud data 0](ex_particles.png "nodes for surface integral (ex_particles.png)")  
![intensity distributions 0](I_example2.png "intensity distributions (I_example2.png)")  


## Analysis sample 3 (in the analysis_sample3)  

This is the analysis result of plane wave scattering by cone objects. 

![mesh image 3](analysis_sample3/cone_m1_image.png "mesh image of the object (analysis_sample3/cone_m1_image.png)") 
![point cloud data 3](analysis_sample3/ex3_particles.png "nodes for surface integral (analysis_sample3/ex3_particles.png)")  
![intensity_distributions 3](analysis_sample3/I_example2_logcb.png "intensity distributions (analysis_sample3/I_example2_logcb.png)")  


## Verification  
The verification results are in the folder verification. 
The analysis result of arrangement of five spheres using "emf_mie_mmls" (in the folder emf_mie_mmls_result) and the analysis result of periodic arrangement of spheres are shown. 

![point cloud data](verification/emf_mie_mmls_result/v1_particles.png "nodes for surface integral (verification/emf_mie_mmls/v1_particles.png)")  


## About mesh file 
This code can use quadrangular ( bi-linear ) and triangular ( linear triangular ) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh file created by Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line ( xxxx.geo is a geometry file). 
The domain number ( Physical Surface ) 99 is assigned to the open region in Gmsh geometry file, becase Gmsh can't use the number 0 (assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.  


## Reference

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)
2. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
3. Command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
4. The calculation program of quasi-periodic Green's function [d3_qpgf_d1](https://github.com/akohta/d3_qpgf_d1/)
5. The electromagnetic field analysis program [emf_mie_mmls](https://github.com/akohta/emf_mie_mmls/)  
6. The data analysis and visualization application [ParaView](https://www.paraview.org/)  
