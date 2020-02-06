# Logistic regression

Julia code for fitting various logistic regression models to the rat behavioral
data from the Poisson clicks task. Allows us to assess L/R bias, as well as
temporal discounting. Uses the GLM.jl package.

In short, fits the following four models:

### 1D click difference
![equation](https://latex.codecogs.com/svg.latex?log%28%5Cfrac%7Bp_R%7D%7B1%20-%20p_R%7D%29%20%3D%20w%20%5CDelta_%7Bclicks%7D%20&plus;%20%5Cbeta)

### 2D, total number of right and left clicks
![equation](https://latex.codecogs.com/svg.latex?log%28%5Cfrac%7Bp_R%7D%7B1%20-%20p_R%7D%29%20%3D%20w_L%20N_L%20&plus;%20w_R%20N_R%20&plus;%20%5Cbeta)

### Click difference across time
![equation](https://latex.codecogs.com/svg.latex?log%28%5Cfrac%7Bp_R%7D%7B1%20-%20p_R%7D%29%20%3D%20w_1%20%5CDelta_1%20&plus;%20w_2%20%5CDelta_2%20&plus;%20...%20&plus;%20w_N%20%5CDelta_N%20&plus;%20%5Cbeta)

### Number of right and left clicks across time
![equation](https://latex.codecogs.com/svg.latex?log%28%5Cfrac%7Bp_R%7D%7B1%20-%20p_R%7D%29%20%3D%20w_%7BR1%7D%20N_%7BR1%7D%20&plus;%20w_%7BR2%7D%20N_%7BR2%7D%20&plus;%20...%20&plus;%20w_%7BL1%7D%20N_%7BL1%7D%20&plus;%20...%20%5Cbeta)

Producing, for each rat, a figure that looks like the following:

![image](figs/K317_frequency_40Hz_example_fig.png)

## To do

- [ ] Major refactor to make the code easier to use
- [ ] Use stratified (across gammas) kfold cv for the estimation of params and their error bars
    - [ ] Make sure that the X matrices contain gammas
- [ ] Check that other generalization metrics also show same trends as AUC(ROC)
- [ ] Visualize how separable the `w * X + beta` distributions are prior to feeding into the logit model
- [ ] Metarat analysis: could fit to one giant collated matrix or average logit weights
- [ ] Explore different bin sizes - are 20Hz rats fit worse bc of low click count per bin?

