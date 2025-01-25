# Darts Empirical Bayesian Hypothesis Testing

## Introduction

The [original Darts-BHT](https://www.nature.com/articles/s41592-019-0351-9) was designed as a first step in a 
deep-learning augmented RNA-seq analysis of splicing variation. The original Darts-BHT framework offers two 
options to specify prior: an informative prior which is from the output of the second step Darts-DNN model,
and an uninformative prior (Darts-BHT-flat). 

In this version, we provide an updated feature borrowing the idea from empirical Bayesian, where it can 
iteratively update the prior starting from uninformative prior. Such updated priors are data-driven and usually
lead to more accurate results. In this implementation, it will repeat until convergence, which usually happens
within five rounds. 

## Install

You can install directly through GitHub.

		library(devtools)
		install_github("quansun98/Darts")

## Instructions

The easiest way to run this version of Darts (v0.2.0) is through calling the two main pipeline wrapper functions,
`Darts()`, for data without replicate, and `Darts_replicate()` for data with replicates.
The input and output look similar for these two functions.

### Input file

Input files are counts from [rMATS-turbo](https://github.com/xinglab/rmats-turbo).
Basically, it requires columns ID,I1,S1,I2,S2,inc\_len,skp\_len. Additional columns are allowed.
An example input file for data without replicate is `test_rmats.txt`, which looks like the following

		ID	I1	S1	I2	S2	inc_len	skp_len
		101	651	51	667	14	2	1
		102	161	113	202	124	2	1

An example input file for data with replicate is `test_rmats_rep.txt`, which looks like the following,
where different replicates are separated by comma:

		ID	I1	S1	I2	S2	inc_len	skp_len
		1	1200,1135,1214	81,180,106	1212,1167,1207	89,69,31	2	1
		2	19,22,25	36,26,19	11,17,15	24,32,35	2	1

### Perform analysis

Detailed documentation of the two functions, `Darts()` and `Darts_replicate()`,
are available in R. In brief, it can start with one simple line of command:

		out = Darts(in_fn = "test_rmats.txt", out_fn = "test_darts.txt", iter=T)

Or for data with replcates:

		out2 = Darts(in_fn = "test_rmats_rep.txt", out_fn = "test_darts_rep.txt", iter=T)


### Output file

The above two functions will write output files as specified, and return the results.
They will return a dataframe of posterior probability. The following columns will be added appending to the original input file.

1. `rho`: The final prior P1, proportion of differential splicing events, from the last iteration to generate this output.

2. `psi1`: Exon inclusion rate (psi values) for the first group.

3. `psi2`: Exon inclusion rate (psi values) for the second group.

4. `mu.mle`: Maximum likelihood estimator of mean psi value for group 1.

5. `delta.mle`: Maximum likelihood estimator of mean difference in psi values in group 2 compared to group 1.

6. `post_pr`: Posterior probability of this event being differential splicing. 
We recommend using `post_pr > 0.9` to define significant differential splicing, 
`post_pr < 0.2` to define no differential splicing, and in between as unsure.

7. `covg`: Mean coverage for this event.



