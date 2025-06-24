# seir-vaccine
 Matlab codes for the sensitivity analysis of an SEIR model with a mass vaccination program
 
 This folder contains the codes for the one-city model, with the following constant parameters:
 (1) N = 5,000,000
 (2) I0 = 5,000, S0 = 4,995,000
 (3) initial values of other compartments are 0
 (4) death rate of vaccinated individuals = 0

  Here is a brief summary of what each file does:
  - onecity.m: (function) right-hand side of the ODE system model
  - eEuler.m: (function) one-step explicit Euler scheme implementation of the model
  - modelsimulation.m: generating plots for the cumulative infected and the cumulative death count for the model, t \in [0,140]
  - lhsu.m: (function) generating the LHS matrix using Latin hypercube sampling
  - modeloutputDt.m: (function) computing the cumulative death count for t = 7,14,....,140 given a matrix of input parameters
  - modeloutputWt.m: (function) computing the cumulative infected count for t = 7,14,...,140 given a matrix of input parameters
  - monotonicitycheckDt.m: checking that D is monotonic with respect to the input parameters
  - monotonicitycheckWt.m: checking that W is monotonic with respect to the input parameters
  - modelstore.m: generating the LHS matrix, computing the output matrices, and storing for use in residualplots.m and sanalysistime.m 
  - residualplots.m: generating the scatter plots of the residuals for the output variables vs the input variables
  - sanalysistime.m: sensitivity analysis implementation using PRCC