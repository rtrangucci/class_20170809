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

## Optional software

* `devtools`
* [RStudio](https://www.rstudio.com) 
  * This will make viewing the Rmd files much easier
* `shinystan`

[Instructions on installing rstan on Mac Linux and Windows](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)

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

## Logistic regression, poststrat in Stan and rstanarm

* Hierarchical logistic regression
  * Priors for hierarchical variance parameters
  * Poststrat in the generated quantities block
  * QR decomposition for nicer posteriors

* Hierarchical logistic regression in rstanarm
  * Easier way to code hierarchical generalized linear models
  * Poststrat using rstanarm
  * Structured priors for interactions

## Multinomial logistic regression

* Hierarchical multinomial logistic regression
  * Model description
  * Fake data generation for model checking
  * Equivalent models
  * Tricks for speed
  * Priors for hierarchical variance parameters
  * Fit to survey data
  * Stan control parameters (adapt-delta, max-treedepth)

# Things to remember when optimizing a model

* Vectorize when possible
* Statistical efficiency and good model fit are the best
antidotes to a slow model (reducing treedepth for each iteration
                           will reduce the time it takes to fit a model by an
                           order of magnitude)
* Use native Stan functions rather than writing your own in the Stan language
  * Most of Stan's functions have been optimized a lot (see linear regression example)
* Can you use sampling independence to reduce the size of your dataset?
* The only expressions that matter are those that are in the kernel of the
posterior. Do not worry about constants!
* Use sufficient statistics to reduce the size of the expression
tree whenever possible
