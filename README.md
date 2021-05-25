# blender_project

Description:  MATLAB code solving for the flow inside a cylindrical domain with rotating floor and stationary walls (a simplified blender). The dimensions of the problem are reduced by assuming that the flow is centrosymmetric about the axis through the center of the cylinder, allowing us to reduce the 3D domain to a 2D rectangular cross-section inside the cylinder extending from the central axis to the edge. Using the streamfunction-vorticity formulation of the incompressible Navier-Stokes equations, we numerically solve (employing similarity transformations on the discretized equations) within the 2D domain [0,r]x[0,\Gamma] over time (with Heun's method) and use centrosymmetry to extrapolate the flow back to the full 3D domain.  


Helpful info for running solver:  

   At line 6, Reynolds number can be changed (1, 100, 1000).
   At line 18, stable timestep values are provided for varying choices of Reynolds number.
   At lines 100-104 and 108-112, the momentum equations can be solved with/without the non-linear terms by commenting in/out appropriately.
   At line 166, there is a block of commented code for saving PNGs at every 50th timestep so that you can make movies.
   At lines 183/186, either error term works.
   At lines 202/203 are the streamfunction and vorticity solutions, if desired.
