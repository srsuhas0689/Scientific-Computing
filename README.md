# Sampling from Probability Distributions in MatLab
This is a sample of my MatLab code related to scientific computing and probability. For the first update I'm defining a couple functions for visualizing and sampling random numbers from distributions.

## Histogram function `hst`
I have defined a function `hst` which creates a histogram using a sampler function and returns the parameters we are interested in. This function also plots the confidence intervals for our bins along with the pdf itself. The function is as follows:
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
I use the `hst` function to plot the histogram of f_E(x) = exp(−x) in the interval `[0, 10]` using 50 bins and
2,500 samples. The green curve and right y-axis show the graph of f_E.
<p align="center">
<img width="732" alt="Histogram of the exponential function." src="https://user-images.githubusercontent.com/65843134/150830562-93975576-13ea-4a4c-8ca2-efab777261b9.png">
</p>
