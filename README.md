# pencil_tools
For Max Planck Astronomy. Tools to process simulations produced by the PENCIL CODE.
A hydrodynamics solver.

This routine should print out the radial difusivity. For now it assumes a 3d simulation with periodic boundries. 

To be called in run directory after the simulation:

```console
python rad_diff.py s st n

```

parameters:
    s:     int pvar file to start from
    st:    int pvar file to stop at 
    n:     int number of particles to track
 
The routine tracks the variance of the particle displacment over time and 
fits it linearly assuming fickian diffusion. the diffusivity is then given by:
    
H * c_s * delta_x = 0.5 * d(sig^2)/dt
    
where sig^2 is the variance of the particle displacment at a given time and 
d(sig^2)/dt is for fickian diffusion the coefficient of the linear fit. A plot 
is preduced of the fit and the particle distrebution as it spreads out with time.
The routine needs suffciently many PVARs.
