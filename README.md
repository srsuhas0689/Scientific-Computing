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
Here I fit randomly generated data to a nonlinear function with parameters `y = f (x; <b>`c`</b>)`. The least-squares fit is the objective function to be minimized. Although this problem is an optimization problem, it can be thought of as an overdetermined system of nonlinear equations for the parameter <b>`c`</b>. 
<p align="center">
<img width="287" alt="function" src="https://user-images.githubusercontent.com/65843134/151388259-8f87d014-d3c2-4586-a060-19c4c6f33b13.png">
</p>
where `f` is an exponentially-damped sinusoidal curve with the unknown parameter <b>`c`</b>. Here <b>`c`</b> is a vector of four scalars:
* Amplitude `c_1`
* Decay `c_2`
* Period `c_3`
* Phase `c_4`

## [Synthetic Data](https://github.com/erik-dali/Scientific-Computing/blob/bd333c7f5431c315695baf0e9f62da17661fea85/optimization.m)
I generate some synthetic data of 100 points randomly and uniformly distributed in `0 <= x <= 10`. Then I plot the true function and the peturbed data where `c = (1, 1/2, 2, 0)`.

<p align="center">
  <img width="659" alt="Nonlinear Fitting" src="https://user-images.githubusercontent.com/65843134/151286963-fc74cf2e-dba2-42b7-aede-0ef108ef6da3.png">
</p>

## [Gauss-Newton Method](https://github.com/erik-dali/Scientific-Computing/blob/bd333c7f5431c315695baf0e9f62da17661fea85/optimization.m)
I implement the Gauss-Newton method with the initial guess of `c0 = (0.95, 0.45, 1.95, 0.05)` and solve the normal equations:
<p align="center">
<img width="490" alt="GN" src="https://user-images.githubusercontent.com/65843134/151390198-4aefa0ef-9f94-4fbe-91e5-36e9c21497f1.png">
</p>
and set the stopping criteria to be 20 iterations. We find the partial derivatives of `f` by hand and compute the Jacobian matrix for each iteration. After 5 iterations the estimates converge to the actual results with full machine accuracy.

If we set the initial guess `c = (1, 1, 1, 1)` then the method diverges and after 8 iterations the program fails because the matrix is no longer positive-definite. Trying other initial guesses we find that the method diverges frequently when the initial guess are not close.
* Theoretically, convergence is only quadratic when our initial guess are close enough to the actual values.
* If we set ε = 0 then we get full machine accuracy after 6 iterations.
* Gauss-Newton method is derived from the Newton method for optimization and it can solve non-linear least squared problems. However, Gauss-Newton method has the advantage of not needing the second derivative which can be difficult to compute. (Source Wikipedia https://en.wikipedia.org/wiki/Gauss– Newton_algorithm)

## [Levenberg-Marquardt Algorithm](https://github.com/erik-dali/Scientific-Computing/blob/bd333c7f5431c315695baf0e9f62da17661fea85/optimization.m)
We set the initial guess to be `c = (1,1,1,1)` for which the Guass-Netwon method fails to converge. Experimentally, if we set `λ1 = 5.11e−3` then after 20 iterations the results converge to full machine accuracy.

<p align="center">
<img width="664" alt="Levenberg-Marquardt Algorithm" src="https://user-images.githubusercontent.com/65843134/151390703-bb0388d9-497f-4b3e-b9fe-88b72731293d.png">
</p>

For larger values of λ1 convergence becomes faster, however after a certain threshold, increasing λ1 does not improve the speed of convergence. It is more helpful to start with a better intial guess than to rely on increasing the magnitude of λ1.

