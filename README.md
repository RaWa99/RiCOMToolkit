# RiCOMToolkit
This repository contains processing programs to support the hydrodynamic model RiCOM,
but can be used with almost any other model by changing inputs.

PROGRAMS
decomp: Grid decomposition using the METIS library.
partition: Grid partitioning for MPI using the METIS library.
reorder: Reorder an element list to minimize bandwidth.
refinex2: Refine a grid by a factor of 2 by subdividing on element midsides.
grid2nc: Convert grid file to netCDF.
rcmha: Harmonic analysis of hourly model output.
rcm2tec: Convert binary model output file to Tecplot ascii or binary input file.
rcm2vtk: Convert binary model output file to VTK ascii or binary input file.
rcm2nc: Convert binary output file to netCDF format.
zoom_sigma: Refine upper and/or lower parts of sigma coordinate profile.
