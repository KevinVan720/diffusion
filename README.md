
# HQ diffusion 

========= This is the Langevin diffusion code which is left by Shanshan. This code can read in OSU hydro information 
and a set of initial position of binary collisions (in which case, we consider it as the initial position of heavy quarks), 
and generate the informaiton of heavy quarks in OSCAR format.

to run this code, extra file: parameter_df.dat (in which some parameter are set)

JetCtl.dat JetData.dat (OSU hydro information) 

initial_xy.dat (initial position of particles)

export ftn10=dNg_over_dt_bC6.dat   // radiative energy loss table

export ftn20=HQ_result.dat

export ftn30=initial_xy.dat

./diffusion
