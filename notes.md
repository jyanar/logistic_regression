## Logistic regression analysis of pbups data

We would like to model the binary choices of the poisson clicks task through the use of a logistic regression model, where each trial consists of an $M$-dimensional input vector $\bar{x}$.  

$\bar{x}$ can either model the click difference at each time bin during the trial or the number of clicks experienced on the left and right sides (in which case it will be of dimensionality $2M$).

We would like to use a logistic regression model because this allows us to map a continuous set of numbers (the click difference, for example) to a binary choice. The use of  the logistic sigmoid function,

$\sigma(a) = \frac{1}{1 + e^{-a}}$

Allows us to map a set of continuous values from $\R$ as we have here to the range $[0,1]$.

The inverse of the logistic sigmoid is given by

$a = ln(\frac{\sigma}{1 - \sigma})$ , and it represents the ratio of probabilities $ln[p(C_1|x)/p(C_2|x)]$ for the two classes, also known as the *log odds.*

### Fitting the logistic regression model

We will use maximum likelihood to determine the parameters of the logistic regression model. To do this, we simply take the derivative of the logistic sigmoid function:

$\frac{d\sigma}{dt} = \sigma(1 - \sigma)$

For a given data set {$\phi_n, t_n$}, 




