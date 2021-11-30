# Compressed sensing:

This project deals with the paper *"Deterministic Sesning Matrices Arising From Near Orthogonal Systems - Shuxing Li and Gennian Ge"*. This paper deals with the concept of compressed sensing. Compressed sensing is a concept that says that a sparse signal can be reconstructed in much fewer samples than what the Nyquist-Shannon sampling theorems suggest. This is possible with the help of sensing matrices. 

Traditionally, the sensing matrices used for compressed sensing are random matrices. This paper introduces several deterministic matrices and compares the recovery percentage of these matrices with that of the traditional matrices.

This project recreates 2 figures from the paper- Singer and MacFarland matrix

## Files included

```
● main.m
● compressed_sensing.m
● sensing_matrix_method.m
● normalized_value.m
● algo_omp.m
● gen_bch_matrix.m
● generate_singer.m
● generate_macfarland.m
```
## Execution

The entire project is run by running `main` without any arguments. The execution time is roughly 10 minutes.

## Structure

* #### main.m
This function is the main function and once called, will call the rest of the defined function to output the required graphs
    
* #### compressed_sensing.m
Function to perform compressed sensing and to plot out all the required graphs
    
* #### sensing_matrix_method.m
Function to recover the output matrix from the input signal and the sensing matrix.
    
* #### normalized_value.m
Function to normalize the graph outputs. This function is used to compensate for the lack of tests (for computational and time reasons). It normalizes the output and ensure that the plot does not deviate a lot (which would usually be normalized with more tests)
    
* #### algo_omp.m
Function algo_omp solves y  Ax, takes the input parameters y, A, k where y is the output field, A is the dataset field and k is the sparsity. It returns the solution of x.
    
* #### gen_bch_matrix.
Function to generage bch matrix
    
* #### generate_singer.m
This function returns a singer matrix of order pxp^ 2 that satisfies the RTP with specified values for k, t, delta.
    
* #### generate_macfarland.m
Function to generate the Macfarland matrix

