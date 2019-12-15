In the top most folder, you will find our submission for task 4.
The program is callable as requested by:
./stencilomp resolution numthreads e.g.
./stencilomp 64 4

The sub folder "/scale_analysis/" contains everything that was used to conduct the scale and parallelization analysis.
We used a version of our program that uses a slightly different calling signature:
./stencilomp resolution schedule e.g.
./stencilomp 64 dynamic
The "schedule" does not actually affect the schedule used for the parallel for loop but is just used to be output to the console afterwards.
The output of this version is also slightly different, so that processing the output into the plots (found in the subfolder "/plots/") via our python script ("plot_byResolution.py") was easier.