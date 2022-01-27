# Intro
This is a collection of my MatLab codes related to scientific computing and probability. 

# 1. Sampling Probability Distributions in MatLab
For the first part I'm defining a couple functions for to draw random samples from probability distributions and plot their histograms. The main function is called `hst` and the function that draws the samples is `sampler`.

## [Histogram function `hst`](https://github.com/erik-dali/Scientific-Computing/blob/d1930d14214fdd64134b8b4d33cd3f77a295084c/histogram-function.m)
I have defined a function `hst` which creates a histogram for a given probability distribution by generating a large number of i.i.d. samples using a sampler function called `sampler`. The histogram function returns a set of parameters that we are interested in. It can also plot the confidence intervals for our bins along with the pdf itself. The function is as follows:
```
[bins, frq, fx, cip, mu, var_dist, var_mean] = hst(@sampler, bins, samples, xmin, xmax, graph)
```
where on the right-hand side the functions takes as input:
* `@sampler` is the function handle of the sampler of a pdf function
* `bins` is the number of bins for the histogram
* `samples` is the number of samples we want to be generated for the histogram
* The argument `graph` takes `y` and `n` to determenie if it should plot the histogram,

and on the left-hand side the function returns:
* `bins` is the vector of bins
* `frq` is the vector of frequencies for each bin
* `fx` is the vector of empirical (numerical) probability distribution i.e. probability of each bin
* `cip` is the vector of probabilities for each bin
* `mu` is the empirical mean
* `var_dist` is the empirical variance
* `var_mean` is the variance of the mean estimator

Note that the last bin will contain all observations larger than the upper bound of the interval.

### Testing the histogram function `hst`
I use the `hst` function to plot the histogram of `f_E(x) = exp(−x)` in the interval `[0, 10]` using 50 bins and
2,500 samples. The green curve and right y-axis show the graph of f_E.
<p align="center">
<img width="732" alt="Histogram of the exponential function." src="https://user-images.githubusercontent.com/65843134/150830562-93975576-13ea-4a4c-8ca2-efab777261b9.png">
</p>

## [Sampling Sum of Two i.i.d. Random Variables](https://github.com/erik-dali/Scientific-Computing/blob/d1930d14214fdd64134b8b4d33cd3f77a295084c/sum-of-random-variables.m)
Let `X_1` and `X_2` be two i.i.d. random variables with pdf functions `f(x) = exp(−x)`. Then the sum of these two random variables is also a random variable `Y = X_1 + X_2` and the pdf for `Y` is the convolution of the pdfs for `X_1` and `X_2`. Therefore we can define a new function to sample from `Y`. Here's how I derived the pdf for `Y` using convolution:
<p align="center">
  <img width="828" alt="Sum of two i.i.d. random variables using convolution" src="https://user-images.githubusercontent.com/65843134/150833408-f2c2d040-7657-43f8-b84b-caddfe13917c.png">
</p>

### Histogram for the new pdf
Using the new sample function derived above, 100 bins, and N = 10^5 samples on the interval `[0,10]` we get the following:
<p align="center">
  <img width="822" alt="Histogram of sum of two i.i.d. random variables" src="https://user-images.githubusercontent.com/65843134/150833810-90d24c6a-2918-4875-8f78-2343466ac295.png">
</p>

# 2. Nonlinear Least-Squares Fitting
Here I fit randomly generated data to a nonlinear function with parameters `y = f (x; c)`. The least-squares fit is the objective function to be minimized. This problem can thought of as an overdetermined system of nonlinear equations for the parameter `c`.
<p align="center">
  <img width="659" alt="Nonlinear Fitting" src="https://user-images.githubusercontent.com/65843134/151286963-fc74cf2e-dba2-42b7-aede-0ef108ef6da3.png">
</p>

the rest is coming soon!
