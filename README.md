2D_Voronoi_Cell_Radial_Alignment (MATLAB CODEs) (Author: Nonthakorn Olaranont and Sarah Pierre)

Voronoi Cell Model for elife paper "Condensation tendency  with planar isotropic actin gradient induces cell radial alignment in the connected contractile tissue".
Here are the steps to run the simulation in the paper:
Open Main.m:
1. Get the inner and outer cells for each Initial Condition:

      a. Each initial condition is labeled IC_625_1 etc. (there are 5)
          
      b. You might need to do some trial and error to get Router right (lines 76-77 in main.m) to get the number of cells in the model to match the experiment (average of 121 cells in each pattern). 
      
      c. In main.m, change K_b to 0.1, lambda1 to 40, and all other lambdas to zero.
      
      d. In solver.m comment out lines 251-4. We want to allow the cells to move.
      
      e. In diff_e_omega.m, uncomment lines 123-135
      
      f. Run main.m, and record/save the lists of inside and outside cells as well as the cdat final configuration. 

2. Get the results of varying contractility and stiffness
      
      a. Uncomment lines 33-36 in main.m and use the inner/outer cell lists from part 1 along with the cdat (line 25).
      
      b. Change K_b to 0, K to 0, lambda1 to 0, lambda 2 &4 to 15. 
      
      c. In create_voronoi, line 164 will be varied from 1, 0.75, 0.5 (as seen in figure 4 in the paper)
      
      d. In solver.m make sure to uncomment lines 251-4 because we don’t want these cells to move now.
      
      e. In diff_e_omega, comment lines 123-135 out. Line 17 will be varied from 0.4, 0.6, 0.8, 1 (as in the paper).
