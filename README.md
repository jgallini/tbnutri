# Tuberculosis Nutritional Intervention Cost-Effectiveness Analysis
This repository contains the code used to simulate a hypothetical nutritional intervention in persons with tuberculosis in India. 

*Author:* Julia Gallini

## Code Overview
The top level program for this project is `Cost_Effectiveness_Results.Rmd`. This
creates an HTML file with all of the main results for this manuscript. It calls
a custom function `cost_eff.R` several times to produce these results. The 
`cost_eff.R` function calls `Transition_function_wide.R`, another custom 
function. Descriptions of each of these functions can be found at the top of
each respective R script.

The `qsub_rmd.sh` file is what was used to submit `Cost_Effectiveness_Results.Rmd`
to the BU Shared Computing Cluster. Please note this code requires significant
time to run; `qsub_rmd.sh` requests 108 hours from the cluster for 10,000 
iterations to be safe. Typically actual run time is closer to 72 hours for 
10,000 iterations.

## Session Requirements
This code was tested under `R/4.2.3`. The code requires the following packages:

- mc2d       
- mvtnorm     
- lhs        
- iterators 
- here   

To install, execute: 
```
install.packages(c("mc2d","mvtnorm","lhs","iterators","here"))
```