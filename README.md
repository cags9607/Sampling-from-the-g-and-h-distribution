# Sampling-from-the-g-and-h-distribution
A method to simulate data from Tukey's g and h distribution. 

The ghrootsv2.R function computes the g and h values required for a desired combination of skewness and kurtosis for Tukey's g and h distribution. By using homotopy continuation, the function also computes several solutions (with different values of variance and mean).

The multgh function simulates data from Tukey's g and h distribution, for specific values of skewness and kurtosis.
