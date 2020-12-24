# Figure 16

The R code to run the methods (and the respective outputs) is in the following folders:

1. Accept/Reject Algorithm Raw Data: "~/AR-Raw"
2. Accept/Reject Algorithm Derivative Data: "~/AR-Derivative"

For the Accept/Reject Algorithm, we ran the code with the following 10 seeds (seed changed within foreach loop):
y, y+10, y+20,..., y+90

where "y" is the study number that was held out. 

The output files of the AR algorithm are in "~/AR-Raw/AR Output Files" and "~/AR-Deriv/AR Output Files"

## Figure Generation

The code that takes that output and generates the Figure is in "FSCV-Between Seed Variability Plot Final.R". Working directories will have to be changed. Methods were run on Harvard Odyssey Cluster ("general" partition).
