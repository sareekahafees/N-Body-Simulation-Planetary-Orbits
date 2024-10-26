# N-Body-Simulation-Planetary-Orbits
This is a simulation of the orbits of 4 solar system bodies around the sun (5 bodies including the sun). It can be modified to simulate the orbits of all 8 planets of the solar system. This was my project for my computational physics course during the last year of my undergraduate degree. The initial task was to simulate the orbits of the inner solar system planets around the sun, but I wrote the code such that it can easily be modified to simulate the trajectories of all 8 planets. 

This code is written in Fortran 90 using the Verlet algorithm and the simulation was tested on Origin.

INPUT FILE: 
There is an input file named 'input2.dat' which includes the positions, velocities and masses of the planet in an array in standard units. The code will convert the standard units to astronomical units.

OUTPUT FILE: 
The data generated by the code is saved in an external txt file named 'output5.dat'. This file can be run on Origin to simulate the trajectories of the planets. The txt file found within this repository already contains enough data to simulate complete orbits of the first 5 planets. If you would like to run the code yourself, please make sure to save the original output5 file separately and rename it so you have the original data for reference. Let the code run for about 1 minute so there is enough data generated to simulate complete trajectories.

EXTENDING THE CODE TO THE GAS GIANTS: 
If you include the gas giants in your code, let the code run for longer than 2 minutes, as the gas giants are massive and take longer to orbit the sun. If the code is terminated early, the simulation (on Origin) will still run, but the trajectories of the gas giants will not be complete.
