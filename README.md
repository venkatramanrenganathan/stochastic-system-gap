# Distance Between Stochastic SISO Linear Systems 
This repository contains the code to simulate and demonstrate distance measure between stochastic linear systems both in the frequency domain and time domain settings. 

 **Associated Paper:** Venkatraman Renganathan and Sei Zhen Khong, `Distance Between Stochastic SISO Linear Systems`, Submitted to the IEEE Transactions on Automatic Control, 2026.

## Dependencies
- Matlab

## Information about code
There are three main coding files:
1. `frequentDomainSimulation.m` 
2. `timeDomainSimulation.m`
3. `comparisonSimulation.m`

## Steps to Run the Code 
1. While running any of the first two coding files, set the parameters of the nominal plant accordingly.
2. Similarly set the distribution information (mean, standard deviation, covariance) about the parameter theta.
3. Find the interested quantity and compare it against the lower/upper bound from the theory in the paper.

## Caution
1. Increasing the sample size (N - number of perturbed plant models) will increase the computation time significantly larger in the order of O(N^2). Hence, choose the sample size according to your available computation resources. 

## Author Information
1. Venkatraman Renganathan, IIT Hyderabad, India. Email: venkatraman@ai.iith.ac.in
2. Sei Zhen Khong, National Sun Yat-Sen University, Taiwan. Email: szkhong@mail.nsysu.edu.tw
   
# Affiliation
1. Venkatraman Renganathan is associated with the department of artificial intelligence at Indian Institute of Technology Hyderabad, India. 
2. Sei Zhen Khong is associated with the department of electrical engineering at the National Sun Yat-Sen University, Taiwan.


