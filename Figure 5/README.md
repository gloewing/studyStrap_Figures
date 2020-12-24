# Figure 5

The R code to run the methods (and the respective outputs) for the raw covariates is in the folder "~/Raw" and for the derivative of the covariates is in the folder "~/Diff"

Within each of these folders are the following folders that each contain code and output for the corresponding algorithms.

1. Accept/Reject Algorithm: "~/AR"
2. OSE:  "~/Ensemble"
3. TOM: "~/Merged"
4. SSE: "~/Study Strap"

For the Accept/Reject Algorithm, we ran the code with the following 10 seeds (seed changed within foreach loop):
y, y+10, y+20,..., y+90

where "y" is the study number that was held out. 

The output files of the AR algorithm are in "~/AR/AR Output Files"

## Figure Generation

The code that takes that output and generates the Figure is in "Box and Whisker Plot of Raw FullLog.R"
