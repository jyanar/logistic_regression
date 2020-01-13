# Logistic regression

We want to model the behavioral data using a logistic regression analysis,
where each side is weighed separately in order to assess how rats pay attention
to evidence on either the left or right side.

That is, we want to fit the following model to each rat:

```
log( p_R / (1 - p_R) ) = \Sum_{t} w_t^L L_t + w_t^R R_t
```


