# EcologicalHRL
A Computational Model of SensoriMotor Learning under Risk of Injury (Hierarchical RL under Ecological Control)

The Matlab code here is meant to test the EcologicalHRL model through an experimental paradigm used in
Babic et al. (expected-to-appear). 
The model was developed on Matlab R2024b; but it should work on most versions.
The model can be tested either
(1) With an individual run emulating a single human participant's sensorimotor learning. 
or
(2) As multiple runs emulating a set of participants going through the experiments.


HOW TO RUN:

(1) For INVIDIAL PARTICIPANT emulation

From the Matlab console go to the directory of the code, and type HIREL_model.
This will apply Ecological Hierarchical learning on an inverted pendulum model 
simulating the Center of Mass (COM) of a standardized human body. The code will 
draw some figures showing the progress and stages of the simulated experiment.
The data will be stored in a global variable which can be used by another 
script to generate multi-participant emulated experiments (see below). 

(2) For MULTI-PARTICIPANT experiment emulation

From the Matlab consolude  go to the directory of the code, and type repeatEXP.
This scripts calls HIREL_model.m to make a virtual experiment over "repeats" many
times (change repeats to a desired value in the code). 
Each repeat is considered analogous to an experiment by a human participant.
All the virtual participant data is stored in a folder with a name constructed from the date 
and time of the execution. The folder can then be analysed, for example by 
show_stats_disc_succnorm_sem.m which generates several plots about the experiments.

HOW TO PLOT MULTI-PARTICIPANT RESULTS:

After running repeatEXP note the name of the data folder name created.
(The folder name will be reported to the user). Then from the Matlab console type
 show_stats_disc_succnorm_sem(datafoldername, num_of_repeats) here num_of_repeats 
must be a number less than or equal to the repeats used in repeatEXP.m
The data featured in the associated paper is based on the multi-participant 
experiment simulation that is stored under keep1_REPS2021_4_21_15_0.
The data plots and results can be regenerated with entering this command from
the Matlab console: show_stats_disc_succnorm_sem("keep1_REPS2021_4_21_15_0",60)


 