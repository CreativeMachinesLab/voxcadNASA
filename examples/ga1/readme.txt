Naive GA 1
====================

Created: June 16, 2014



Description

Naive GA (direct-encoding single hillclimber),
to be used as basic example for future GAs
which wish to interface with Voxelyze for
soft robot evolution.


Library dependencies:

o  voxelyze.0.9.92 
(note: must be build into file named: voxelize,
and moved to, or linked to in, current directory)



Useage:

Parameters are defined within "main.cpp" and "Individual.cpp"
A GA intended for actual use should create a "config.txt" file 
that contains all relevant parameters and their meanings for your GA.

$ ./naiveGA



To Make Executables:

$ make



Outputs:

o  fitness over time recorded in "champFitness.txt"

o  VXA (input type for Voxelyze) files kept as indicated by parameter in "main.cpp"
   labeled as "VoxelyzeGA_Gen=[generation]_Fit=[fitness].vxa"

   These VXA files can be used in VoxCad GUI to visualize evolved creatures
   FILE -> IMPORT -> SIMULATION -> [vxaFileName].vxa



Questions:

Contact Nick Cheney at nac93@cornell.edu



