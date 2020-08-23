# Treatment-Effect-Modification-with-Non-monotone-Missingness
This repository contains the analysis code used in paper 'Treatment-Effect-Modification-with-Non-monotone-Missingness'. We provide the code to generate and analyze a simulated data with similar design as the real dengue data used in the paper in this repository and a step-by-step demostration on how to use the code to obtain results of this simulated data set and make the VE curves as shown in Figure 1 in the supplementary material. 

## Step 1: Simulate Data
The generateData function provided in myFunctions.R generates a simulated data with a similar design as the real dengue data used in the paper.

## Step 2: Calculate VE curves and their CIs and CBs
The CalculateVE.R demostrates how to calculate the marginal, baseline seropositive and baseline seronegative VE curves and their corresponding 95% confidence intervals and confidences bands based on 500 perturbation iterations, and save the results into a RData file. 

## Step 3: Plot the VE Curves
The VEPlot.R reads the results from the RData.file and plot the individual VE curves and their corresponding 95% confidence intervals and confidences bands as shown in Figure 1 in the supplementary material.
