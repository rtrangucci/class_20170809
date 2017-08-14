## Course Description

This repo contains all of the course materials for
the 1-day Stan course on 08/09/2017. 

## Instructor

I'm Rob Trangucci, and I've been a Stan developer for the 
past 3 years. I've been funded by YouGov to do general 
research and software development for the Stan library.
I focus on hierarchical modeling, Gaussian processes, 
and writing the matrix functions for the [Stan math library](https://github.com/stan-dev/math).

Email: <robert.trangucci@gmail.com>
GitHub: <https://github.com/rtrangucci>

# Prerequisite software

## Necessary software

Come to class with [R](https://cran.r-project.org)  installed on your machines
and the following packages installed:

* `rstan`
* `rstanarm`
* `bayesplot`
* `dplyr`
* `ggplot2`

You can also just run this script to install the packages that you don't have.

```
if (!require("rstanarm")) install.packages("rstanarm")
if (!require("bayesplot")) install.packages("bayesplot")
```

These packages will install `ggplot2` and `dplyr` if they're not
already installed. As for installing `rstan`, Mac and Linux users
with recent OS versions and Xcode installed (for Mac users) should
have no problem using install.packages("rstan"). If you do run
into problems, follow [these instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Mac-or-Linux).

Windows users should follow [these instructions](https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Windows)
## Optional software

* `devtools`
* [RStudio](https://www.rstudio.com) 
  * This will make viewing the Rmd files much easier
* `shinystan`


# Class format

I'll walk through a presentation from 9:30 to 10:30 that'll introduce everyone
to Stan, and some bedrock principles needed to understand how Stan works. Then
we'll walk through four RMarkdown documents, which you'll have access to, and
we'll work through writing and fitting models to both real and generated data.
I'll highlight key takeaways for each example and how they apply to running
models over real data. When the class is over you'll have Stan and RStanArm
code that'll generalize across datasets, along with all the data prep code and
commentary in the RMarkdown docs.

# Schedule

 Time | Topic |
 -----| ------
 09:00 AM to 09:30 AM | Intro and coffee                                  
 09:30 AM to 11:30 AM | Hierarchical models in Stan                      
 11:30 AM to 12:00 PM | Coffee and questions                              
 12:00 PM to 01:30 PM | Logistic regression, postrat in Stan and rstanarm 
 01:30 PM to 02:00 PM | Lunch                                             
 02:30 PM to 04:30 PM | Multinomial logit in Stan                        
 04:30 PM to 05:30 PM | Review and questions                              

We'll hit the following key points in each of the sections

## Hierarchical models in stan

* One-way normal model
  * How to write a Stan model
  * What Stan does under the hood
  * Sampler output and diagnostics
  * What to do about divergences? (Non-centered parameterization)
  * How to modify the target density manually for speed (only if there is time)

* Linear regression with group-varying intercepts and coefficients
  * Priors for group-varying intercepts and coefficients
  * Priors for covariance matrices
  * LKJ prior for correlation matrices

## Logistic regression, poststrat in Stan and rstanarm

* Hierarchical logistic regression
  * Stan control parameters (adapt-delta, max-treedepth)
  * Poststrat in the generated quantities block
  * QR decomposition for nicer posteriors
  * Priors for hierarchical variance parameters

* Hierarchical logistic regression in rstanarm
  * Easier way to code hierarchical generalized linear models
  * Poststrat using rstanarm
  * Posterior predictive checks (PPC) in `bayesplot`
  * Structured priors for interactions

## Multinomial logistic regression

* Hierarchical multinomial logistic regression
  * Model description
  * Fake data generation for model checking
  * Equivalent models
  * Tricks for speed
  * Priors for hierarchical variance parameters
  * Fit to survey data

# Things to remember when optimizing a model

* Vectorize when possible
* Statistical efficiency and good model fit are the best antidotes to a slow
model (reducing treedepth for each iteration will reduce the time it takes to
       fit a model by an order of magnitude)
* Use native Stan functions rather than writing your own in the Stan language
  * Most of Stan's functions have been optimized a lot (see linear regression example)
* Can you use sampling independence to reduce the size of your dataset?
* The only expressions that matter are those that are in the kernel of the
posterior. Do not worry about constants!
* Use sufficient statistics to reduce the size of the expression tree whenever
possible (examples are IID normal observations at any level linear models)

# Useful links

* [Explanation of `state_pres_vote` data](http://www.slate.com/articles/technology/future_tense/2016/11/the_polls_of_the_future_will_be_reproducible_and_open_source.html)
* [Stan case study about why we fit hierarchical models (uses baseball data)](pool-binary-trials/pool-no-pool.Rmd)
  * [Knitted hmtl](pool-binary-trials/pool-no-pool.html)  
# Literature folder table of contents

* [Intuition behind different covariance priors](literature/Visualization-Covariance-Matrices.pdf)
  * Information about ways to parameterize a covariance matrix prior and how to evaluate the properties of different priors. Lighter on math, lots of nice graphs to give a better idea as to what is happening at high dimensions.
* [Conceptual introduction to Hamiltonian Monte Carlo](literature/Betancourt-Conceptual-introduction-to-Hamiltonian-Monte-Carlo.pdf)
  * Wonderful intro to why HMC is successful where other samplers fail and how and when HMC fails. Also lighter on math, lots of helpful graphics.
* [Challenges sampling from Hierarchical models](literature/Betancourt_Girolami_Hierarchical_Models.pdf)
  * In-depth treatment of non-centered parameterization for hierarchical models. Lighter on math.
