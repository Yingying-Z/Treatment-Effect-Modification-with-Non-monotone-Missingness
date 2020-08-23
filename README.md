# Treatment-Effect-Modification-with-Non-monotone-Missingness
This repository contains the analysis code used in paper 'Treatment-Effect-Modification-with-Non-monotone-Missingness'. We provide the code to generate and analyze a simulated data with similar design as the real dengue data used in the paper in the repository and a step-by-step demostration on how to use the code to obtain results of this simulated data set and make the the VE curves as shwon in Figure 1 in the supplementary material. 

## Step 1: Simulate Data
The generateData function provided in the myFunctions.R generated a simulated data with similar design as the real dengue data used in the paper.

## Step 2: Calculate the marginal, baseline seropositive and baseline seronegative VE curves and their corresponding 95% CIs and CBs based on 500 perturbation iterations
The CalculateVE.R demostrates how to calculate the marginal, baseline seropositive and baseline seronegative VE curves and their corresponding95% CIs and CBs based on 500 perturbation iterations, and save the results into a RData file. 

## Step 3: plot the individual VE curves and their  corresponding 95% CIs and CBs as shown in Figure 1 in the supplementary material
The VEPlot.R reads the results from the RData.file and generates the VE plots as shown in Figure 1 in the supplementary material.
