cvAUC.plus
==========

This small package contains a few simple functions that enhance the [cvAUC](https://github.com/ledell/cvAUC) package created by Erin LeDell. Notably, it allows for confidence intervals and hypothesis tests about differences between cross-validated AUC for two different prediction algorithms.

Installation of `cvAUC.plus`
----------------------------

The package can be installed directly from GitHub.

``` r
devtools::install_github("benkeser/cvAUC.plus")
library(cvAUC.plus)
```

Using wrap\_cvAUC
-----------------

The first demonstration shows how to define a learner that can be used with `wrap_cvAUC`, which will compute the cross-validated AUC for that learner. The learner that is passed to `wrap_AUC` is a function of a specific form. If the user is familiar with the format of wrappers passed to the `SuperLearner` function in the package of the same name, then defining these wrappers should be fast and easy. A proper wrapper for `wrap_cvAUC` should be a function that takes as input `Y`, `X`, and `newX`. The function estimates a prediction function using `Y` and `X` and returns predictions on `newX`. The output of the function should be a list with two entries: `fit` the prediction model fit (can be `NULL` if you'd like) and `pred` the predictions on `newX`. Here are two simple examples.

``` r
# a simple main terms GLM
myglm <- function(Y, X, newX){
    fm <- glm(Y~., data = X, family = binomial())
    pred <- predict(fm, newdata = newX, type = "response")
    return(list(fit = fm, pred = pred))
}
# a random forest
library(randomForest)
myrf <- function(Y, X, newX){
    require(randomForest)
    fm <- randomForest(x = X, y = factor(Y), xtest = newX)
    pred <- fm$test$votes[,2]
    return(list(fit=fm, pred = pred))
}
```

Now we can use `wrap_cvAUC` to estimate the cross-validated AUC for these algorithms.

``` r
# simulate data
n <- 500
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rbinom(n,1,0.5))
Y <- with(X, rbinom(n, 1, plogis(x1*x3 + x2^2/x1 + x3)))

# get CV-AUC of main terms GLM
auc_glm <- wrap_cvAUC(Y = Y, X = X, learner = "myglm", seed = 123)
auc_glm
```

    ## 10-fold CV-AUC for myglm      : 0.76
    ## 95% confidence interval        : 0.71 , 0.80
    ## One-sided p-value CV-AUC > 0.5 : < 0.01

``` r
# get CV-AUC of main terms 
auc_rf <- wrap_cvAUC(Y = Y, X = X, learner = "myrf", seed = 123)
auc_rf
```

    ## 10-fold CV-AUC for myrf      : 0.82
    ## 95% confidence interval        : 0.79 , 0.86
    ## One-sided p-value CV-AUC > 0.5 : < 0.01

Using diff\_cvAUC
-----------------

The main addition over the `cvAUC` package is the ability to test for differences in cross-validated AUC between two different model fits. This is achieved via the `diff_cvAUC` function.

``` r
# compare random forest to GLM
diff_auc <- diff_cvAUC(fit1 = auc_rf, fit2 = auc_glm)
diff_auc
```

    ## Diff 10-fold CV-AUC (myrf - myglm) : 0.07
    ## 95% confidence interval               : 0.02 , 0.12
    ## Two-sided p-value diff not 0          : = 0.01

Options for splitting
---------------------

The `wrap_cvAUC` provides some functionality to control the sample splitting, based on the code developed in the `SuperLearner` package by Eric Polley. You can check the documentation for `SuperLearner.CV.control` to see all of this functionality, as internally this is the function that is called by `wrap_AUC`.

Parallelization
---------------

The `wrap_cvAUC` function allows for parallelization via `foreach` if `parallel = TRUE`. This will parallelize the fitting of the `learner` algorithm over folds. The function will internally start (and stop) a cluster using `detectCores()` to determine how many cores to use.

``` r
# simulate data
n <- 2000
X <- data.frame(x1 = rnorm(n), x2 = rnorm(n), x3 = rbinom(n,1,0.5))
Y <- with(X, rbinom(n, 1, plogis(x1*x3 + x2^2/x1 + x3)))

# non-parallel
system.time(tmp <- wrap_cvAUC(Y = Y, X = X, learner = "myrf", seed = 123))
```

    ##    user  system elapsed 
    ##   2.375   0.029   2.477

``` r
# parallel
system.time(tmp <- wrap_cvAUC(Y = Y, X = X, learner = "myrf", seed = 123, parallel = TRUE))
```

    ##    user  system elapsed 
    ##   0.193   0.016   3.071
