# MatrixAnalyticalSamplingFormula
## A matrix-analytical sampling formula for time-homogeneous coalescent processes under the infinite sites mutation model
This README file explains how the R code for calculating the probability of a sample in a population genetic model is organized.

The sampling formula in Hobolth et al (2024) is implemented in step 1-3 below.<br>
Step 1 describes how an evolutionary model is encoded as an object.<br>
Step 2 converts the evolutionary model into a probability model.<br>
Step 3 calculates the probability of a sample for the evolutionary model


In example A-C we give three examples of the sampling formla.<br>
Example A: Probability for a site frequency spectrum (SFS) with sample size 3.<br>
Example B: Probability for a SFS with arbitrary sample size for the Kingman coalescent.<br>
Example C: Probability for a SFS with arbitrary sample size for the Beta-coalescent.

#### Step 1. Evolutionary Model 
The basic object is an evolutionary model.<br>
We encode an evolutionary model as an R object determined by a list of four components:<br>
$initVec: Initial probability vector <br>
$SMat: Sub-intensity matrix <br>
$RMat: Reward matrix <br>
$mutaRate: Mutation rate <br>
An evolutionary model can be coded directly (as described in 'Example A. A simple example' below), but we also provide models for analysis of the site frequency spectrum (as described in 'Example B. Application: Analysis of SFS' below) or Beta-coalescent (as described in 'Example C. Application: Analysis of Beta-coalescent' below).

#### Step 2. Probability Model 
The evolutionary model (object) is then converted into a probability model (object).<br>
We encode the probability model as an R object determined by a list of three components:<br>
$initVec: Initial probability vector<br>
$PrbArray: Array of sub-transition probability matrices<br> 
$exitVec: Exit probability vector<br>
The conversion follows the main Theorem in Hobolth et al (2024).<br>
File: EvoModelToPrbModel.R<br>
Function: EvoMdlToPrbMdl( EvolMdl )

#### Step 3. Sampling Probability
We use the main Theorem in Hobolth et al (2024) to calculate the probability of a sample for a given model.<br>
File: PrbOfSample.R<br>
Function: SmplPrb( Smpl,PrbMdl ) 

#### Example A. A simple example 
As an illustration of the set-up described above we provide a simple example where we analyze the SFS for a sample of size n=3. <br>
Please see Example 3 in Section 2 in Hobolth et al (2024) for a detailed description of the example.<br>
File: SFSThreeAnalysis.R 

#### Example B. Application: Analysis of SFS 
A general application of the SFS is provided. We first define the rate matrix and state space for the block counting process (SFSCoal.R) and second turn this information into an evolutionary model (SFSModel.R). Finally, the sampling probabilities of the site frequency spectra are analyzed.<br>
Files: SFSCoal.R, SFSModel.R, SFSAnalysis.R   

#### Example C. Application: Analysis of Beta-coalescent
A general applicattion of the SFS for the Beta-coalescent is provided. We first define the rate matrix and state space for the Beta-coalescent (BetaCoal.R) and second turn this information into an evolutionary model (BetaModel.R). Finally, the sampling probabilities of the Beta-coalescent are analyzed.<br>
Files: BetaCoal.R, BetaModel.R, BetaAnalysis.R
  
### References 
Hobolth, A, Futschik, A, Boitard, S, and Leblois, R (2024).
A matrix-analytical sampling formula for time-homogeneous coalescent processes under the infinite sites mutation model.
