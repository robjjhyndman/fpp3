# 5.1 A tidy forecasting workflow
The process of producing forecasts for time series data can be broken down into a 
few steps.
            
#             |->   Specify   ->|            
# Tidy -> Visualize         Estimate  ->  Forecast
#             |<-   Evaluate  <-|
              
# 5.1.1 Data preparation (tidy)
The first step in forecasting is to prepare data in the correct format. This process 
may involve loading in data, identifying missing values, filtering the time series, and 
other pre-processing tasks.

Many models have different data requirements; some require the series to be in time 
order, others require no missing values. Checking your data is an essential step to 
understanding its features and should always be done before models are estimated.

# 5.1.2 Plot the data (visualise)
Visualisation is an essential step in understanding the data. Looking at your data 
allows you to identify common patterns, and subsequently specify an appropriate model.

# 5.1.3 Define a model (specify)
There are many different time series models that can be used for forecasting, and 
much of this book is dedicated to describing various models. Specifying an appropriate 
model for the data is essential for producing appropriate forecasts.

# 5.1.4 Train the model (estimate)
Once an appropriate model is specified, we next train the model on some data.

# 5.1.5 Check model performance (evaluate)
Once a model has been fitted, it is important to check how well it has performed on 
the data.

# 5.1.6 Produce forecasts (forecast)
With an appropriate model specified, estimated and checked, it is time to produce 
the forecasts.




# 5.2 Some simple forecasting methods

# 5.2.1 Mean method
Here, the forecasts of all future values are equal to the average (or “mean”) of 
the historical data. If we let the historical data be denoted by y1,…,yT, then we 
can write the forecasts as

# ^y_T+h|T = ¯y = (y1 + ⋯ + yT)/T

The notation ^y_T+h|T is a short-hand for the estimate of y_T+h based on the data  
y1,…,yT.

# 5.2.2 Naïve method
For naïve forecasts, we simply set all forecasts to be the value of the last 
observation. That is,

# ^y_T+h|T = yT

This method works remarkably well for many economic and financial time series.

# 5.2.3 Seasonal naïve method
A similar method is useful for highly seasonal data. In this case, we set each 
forecast to be equal to the last observed value from the same season (e.g., the 
same month of the previous year). Formally, the forecast for time T+h is written as

# ^y_T+h|T = y_T+h−m(k+1)

where m = the seasonal period, and k is the integer part of (h−1)/m (i.e., the 
number of complete years in the forecast period prior to time T+h). This looks 
more complicated than it really is. For example, with monthly data, the forecast 
for all future February values is equal to the last observed February value. With 
quarterly data, the forecast of all future Q2 values is equal to the last observed 
Q2 value (where Q2 means the second quarter). Similar rules apply for other months 
and quarters, and for other seasonal periods.

# 5.2.4 Drift method
A variation on the naïve method is to allow the forecasts to increase or decrease 
over time, where the amount of change over time (called the drift) is set to be 
the average change seen in the historical data. Thus the forecast for time  
T+h is given by

# ^y_T+h|T = yT + (h/(T-1)) * sigma(yt - y_t-1, t=2->T) = yT + h((yT - y1)/(T-1))

This is equivalent to drawing a line between the first and last observations, and 
extrapolating it into the future.




# 5.3 Fitted values and residuals

# 5.3.1 Fitted values
Each observation in a time series can be forecast using all previous observations. 
We call these fitted values and they are denoted by ^y_t|t−1, meaning the forecast of  
yt based on observations y1,…,y_t−1. We use these so often, we sometimes drop part 
of the subscript and just write ^yt instead of ^y_t|t−1. Fitted values almost always 
involve one-step forecasts.

Actually, fitted values are often not true forecasts because any parameters involved 
in the forecasting method are estimated using all available observations in the time 
series, including future observations. For example, if we use the mean method, the 
fitted values are given by

# ^yt = ^c

where ^c is the average computed over all available observations, including those 
at times after t. Similarly, for the drift method, the drift parameter is estimated 
using all available observations. In this case, the fitted values are given by

# ^yt = y_t−1 + ^c

where ^c = (yT − y1)/(T−1). In both cases, there is a parameter to be estimated 
from the data. The “hat” above the c reminds us that this is an estimate. When 
the estimate of c involves observations after time t, the fitted values are not 
true forecasts. On the other hand, naïve or seasonal naïve forecasts do not involve 
any parameters, and so fitted values are true forecasts in such cases.

# 5.3.2 Residuals
The “residuals” in a time series model are what is left over after fitting a model. 
The residuals are equal to the difference between the observations and the 
corresponding fitted values:
  
# e_t = yt − ^yt

If a transformation has been used in the model, then it is often useful to look 
at residuals on the transformed scale. We call these “innovation residuals”. 
For example, suppose we modelled the logarithms of the data, wt = log(yt). Then 
the innovation residuals are given by w_t − ^w_t whereas the regular residuals 
are given by yt − ^yt. If no transformation has been used then the innovation 
residuals are identical to the regular residuals, and in such cases we will simply 
call them “residuals”.

Residuals are useful in checking whether a model has adequately captured the information 
in the data. For this purpose, we use innovation residuals.

If patterns are observable in the innovation residuals, the model can probably 
be improved.




# 5.4 Residual diagnostics
A good forecasting method will yield innovation residuals with the following properties:
  
1. The innovation residuals are uncorrelated. If there are correlations between 
innovation residuals, then there is information left in the residuals which should 
be used in computing forecasts.
2. The innovation residuals have zero mean. If they have a mean other than zero, 
then the forecasts are biased.

Any forecasting method that does not satisfy these properties can be improved. 
However, that does not mean that forecasting methods that satisfy these properties 
cannot be improved. It is possible to have several different forecasting methods 
for the same data set, all of which satisfy these properties. Checking these properties 
is important in order to see whether a method is using all of the available information, 
but it is not a good way to select a forecasting method.

If either of these properties is not satisfied, then the forecasting method can 
be modified to give better forecasts. Adjusting for bias is easy: if the residuals 
have mean m, then simply add m to all forecasts and the bias problem is solved.

In addition to these essential properties, it is useful (but not necessary) for 
the residuals to also have the following two properties.

3. The innovation residuals have constant variance. This is known as “homoscedasticity”.
4. The innovation residuals are normally distributed.

These two properties make the calculation of prediction intervals easier. However, 
a forecasting method that does not satisfy these properties cannot necessarily be 
improved. Sometimes applying a Box-Cox transformation may assist with these properties, 
but otherwise there is usually little that you can do to ensure that your innovation 
residuals have constant variance and a normal distribution. Instead, an alternative 
approach to obtaining prediction intervals is necessary.

# 5.4.1 Portmanteau tests for autocorrelation
In addition to looking at the ACF plot, we can also do a more formal test for 
autocorrelation by considering a whole set of r_k values as a group, rather than 
treating each one separately.

Recall that r_k is the autocorrelation for lag k. When we look at the ACF plot to 
see whether each spike is within the required limits, we are implicitly carrying 
out multiple hypothesis tests, each one with a small probability of giving a false 
positive. When enough of these tests are done, it is likely that at least one will 
give a false positive, and so we may conclude that the residuals have some remaining 
autocorrelation, when in fact they do not.

In order to overcome this problem, we test whether the first ℓ autocorrelations are 
significantly different from what would be expected from a white noise process. A 
test for a group of autocorrelations is called a portmanteau test, from a French 
word describing a suitcase or coat rack carrying several items of clothing.

One such test is the Box-Pierce test, based on the following statistic

# Q = T*sigma((r_k)^2, k=1->ℓ)

where ℓ is the maximum lag being considered and T is the number of observations. If each  
r_k is close to zero, then Q will be small. If some r_k values are large (positive 
or negative), then Q will be large. We suggest using ℓ=10 for non-seasonal data and  
ℓ=2m for seasonal data, where m is the period of seasonality. However, the test is 
not good when ℓ is large, so if these values are larger than T/5, then use ℓ=T/5

A related (and more accurate) test is the Ljung-Box test, based on

# Q* = T(T+2)*sigma((T-k)^(-1) * (r_k)^2, k=1->ℓ)

Again, large values of Q* suggest that the autocorrelations do not come from a white 
noise series.

How large is too large? If the autocorrelations did come from a white noise series, then both  
Q and Q* would have a χ^2 distribution with ℓ degrees of freedom.




# 5.5 Distributional forecasts and prediction intervals

# 5.5.1 Prediction intervals
A prediction interval gives an interval within which we expect yt to lie with a 
specified probability. For example, assuming that distribution of future observations 
is normal, a 95% prediction interval for the h-step forecast is

# ^y_T+h|T ± 1.96 * ^σ_h

where ^σ_h is an estimate of the standard deviation of the h-step forecast distribution.

More generally, a prediction interval can be written as

# ^y_T+h|T ± c * ^σ_h

where the multiplier c depends on the coverage probability.

The value of prediction intervals is that they express the uncertainty in the forecasts. 
If we only produce point forecasts, there is no way of telling how accurate the forecasts 
are. However, if we also produce prediction intervals, then it is clear how much uncertainty 
is associated with each forecast. For this reason, point forecasts can be of almost 
no value without the accompanying prediction intervals.

# 5.5.2 One-step prediction intervals
When forecasting one step ahead, the standard deviation of the forecast distribution 
can be estimated using the standard deviation of the residuals given by

# ^σ = sqrt((1/(T-K-M)) * sigma((e_t)^2, t=1->T)) (5.1)

where K is the number of parameters estimated in the forecasting method, and M
is the number of missing values in the residuals. (For example, M=1 for a naive 
forecast, because we can’t forecast the first observation.)

# 5.5.3 Multi-step prediction intervals
A common feature of prediction intervals is that they usually increase in length 
as the forecast horizon increases. The further ahead we forecast, the more uncertainty 
is associated with the forecast, and thus the wider the prediction intervals. That is,  
σ_h usually increases with h (although there are some non-linear forecasting methods 
which do not have this property).

To produce a prediction interval, it is necessary to have an estimate of σ_h. As 
already noted, for one-step forecasts (h=1), Equation (5.1) provides a good estimate 
of the forecast standard deviation σ1. For multi-step forecasts, a more complicated 
method of calculation is required. These calculations assume that the residuals 
are uncorrelated.

# 5.5.4 Benchmark methods
For the four benchmark methods, it is possible to mathematically derive the forecast 
standard deviation under the assumption of uncorrelated residuals. If ^σ_h denotes 
the standard deviation of the h-step forecast distribution, and ^σ is the residual 
standard deviation given by (5.1), then we can use the expressions shown in 
Table 5.2. Note that when h=1 and T is large, these all give the same approximate 
value ^σ.


Table 5.2: Multi-step forecast standard deviation for the four benchmark methods, where  
σ is the residual standard deviation, m is the seasonal period, and k
is the integer part of (h−1)/m (i.e., the number of complete years in the forecast 
period prior to time T+h).

# Benchmark method                h-step forecast standard deviation
# Mean                            ^σ_h = ^σ * sqrt(1 + 1/T)
# Naïve                           ^σ_h = ^σ * sqrt(h)
# Seasonal naïve                  ^σ_h = ^σ * sqrt(k+1)
# Drift                           ^σ_h = ^σ * sqrt(h(1 + h/(T-1)))

# 5.5.5 Prediction intervals from bootstrapped residuals
When a normal distribution for the residuals is an unreasonable assumption, one 
alternative is to use bootstrapping, which only assumes that the residuals are 
uncorrelated with constant variance. We will illustrate the procedure using a 
naïve forecasting method.

A one-step forecast error is defined as e_t = yt − ^y_t|t−1. For a naïve forecasting 
method, ^y_t|t−1 = y_t−1, so we can rewrite this as

# yt = y_t−1 + e_t

Assuming future errors will be similar to past errors, when t > T we can replace  
e_t by sampling from the collection of errors we have seen in the past (i.e., the 
residuals). So we can simulate the next observation of a time series using

# (y*)_T+1 = yT + (e*)_T+1

where (e*)_T+1 is a randomly sampled error from the past, and (y*)_T+1 is the possible 
future value that would arise if that particular error value occurred. We use a * 
to indicate that this is not the observed y_T+1 value, but one possible future that 
could occur. Adding the new simulated observation to our data set, we can repeat 
the process to obtain

# (y*)_T+2 = (y*)_T+1 + (e*)_T+2

where (e*)_T+2 is another draw from the collection of residuals. Continuing in this way, 
we can simulate an entire set of future values for our time series.




# 5.6 Forecasting using transformations
When forecasting from a model with transformations, we first produce forecasts of 
the transformed data. Then, we need to reverse the transformation (or back-transform) 
to obtain forecasts on the original scale. For Box-Cox transformations given by (3.1), 
the reverse transformation is given by

# yt = { exp(wt) if λ=0;
#      { sign(λ*w_t + 1) * |λ*w_t + 1|^(1/λ) otherwise

# 5.6.1 Prediction intervals with transformations
If a transformation has been used, then the prediction interval is first computed 
on the transformed scale, and the end points are back-transformed to give a prediction 
interval on the original scale. This approach preserves the probability coverage of 
the prediction interval, although it will no longer be symmetric around the point forecast.




# 5.7 Forecasting with decomposition
Assuming an additive decomposition, the decomposed time series can be written as

# yt = ^S_t + ^A_t

where ^A_t = ^T_t + ^R_t is the seasonally adjusted component. Or, if a multiplicative 
decomposition has been used, we can write

# yt = ^S_t * ^A_t

where ^A_t = ^T_t * ^R_t.

To forecast a decomposed time series, we forecast the seasonal component, ^S_t, 
and the seasonally adjusted component ^A_t, separately. It is usually assumed that 
the seasonal component is unchanging, or changing extremely slowly, so it is forecast 
by simply taking the last year of the estimated component. In other words, a seasonal 
naïve method is used for the seasonal component.

To forecast the seasonally adjusted component, any non-seasonal forecasting method 
may be used. For example, the drift method, or Holt’s method, or a non-seasonal 
ARIMA model, may be used.




# 5.8 Evaluating point forecast accuracy

# 5.8.1 Training and test sets
It is important to evaluate forecast accuracy using genuine forecasts. Consequently, 
the size of the residuals is not a reliable indication of how large true forecast 
errors are likely to be. The accuracy of forecasts can only be determined by considering 
how well a model performs on new data that were not used when fitting the model.

When choosing models, it is common practice to separate the available data into 
two portions, training and test data, where the training data is used to estimate 
any parameters of a forecasting method and the test data is used to evaluate its 
accuracy. Because the test data is not used in determining the forecasts, it should 
provide a reliable indication of how well the model is likely to forecast on new 
data.

The size of the test set is typically about 20% of the total sample, although this 
value depends on how long the sample is and how far ahead you want to forecast. The 
test set should ideally be at least as large as the maximum forecast horizon required. 
The following points should be noted.

- A model which fits the training data well will not necessarily forecast well.
- A perfect fit can always be obtained by using a model with enough parameters.
- Over-fitting a model to data is just as bad as failing to identify a systematic 
pattern in the data.

# 5.8.2 Forecast errors
A forecast “error” is the difference between an observed value and its forecast. 
Here “error” does not mean a mistake, it means the unpredictable part of an 
observation. It can be written as

# e_T+h = y_T+h − ^y_T+h|T

where the training data is given by {y1,…,yT} and the test data is given by 
{y_T+1, y_T+2,…}.

Note that forecast errors are different from residuals in two ways. First, residuals 
are calculated on the training set while forecast errors are calculated on the test 
set. Second, residuals are based on one-step forecasts while forecast errors can 
involve multi-step forecasts.

We can measure forecast accuracy by summarising the forecast errors in different ways.

# 5.8.3 Scale-dependent errors
The forecast errors are on the same scale as the data. Accuracy measures that are 
based only on e_t are therefore scale-dependent and cannot be used to make comparisons 
between series that involve different units.

The two most commonly used scale-dependent measures are based on the absolute errors 
or squared errors:
  
# Mean absolute error: MAE = mean(|e_t|),
# Root mean squared error: RMSE = sqrt(mean((e_t)^2))
  
When comparing forecast methods applied to a single time series, or to several 
time series with the same units, the MAE is popular as it is easy to both understand 
and compute. A forecast method that minimises the MAE will lead to forecasts of the 
median, while minimising the RMSE will lead to forecasts of the mean. Consequently, 
the RMSE is also widely used, despite being more difficult to interpret.

# 5.8.4 Percentage errors
The percentage error is given by p_t = 100e_t/yt. Percentage errors have the advantage 
of being unit-free, and so are frequently used to compare forecast performances between 
data sets. The most commonly used measure is:

# Mean absolute percentage error: MAPE = mean(|p_t|).
  
Measures based on percentage errors have the disadvantage of being infinite or undefined if  
yt = 0 for any t in the period of interest, and having extreme values if any yt 
is close to zero. Another problem with percentage errors that is often overlooked 
is that they assume the unit of measurement has a meaningful zero. For example, 
a percentage error makes no sense when measuring the accuracy of temperature forecasts 
on either the Fahrenheit or Celsius scales, because temperature has an arbitrary zero point.




# 5.10 Time series cross-validation
A more sophisticated version of training/test sets is time series cross-validation. 
In this procedure, there are a series of test sets, each consisting of a single observation. 
The corresponding training set consists only of observations that occurred prior to 
the observation that forms the test set. Thus, no future observations can be used in 
constructing the forecast. Since it is not possible to obtain a reliable forecast 
based on a small training set, the earliest observations are not considered as test 
sets.

The forecast accuracy is computed by averaging over the test sets. This procedure 
is sometimes known as “evaluation on a rolling forecasting origin” because the 
“origin” at which the forecast is based rolls forward in time.

With time series forecasting, one-step forecasts may not be as relevant as multi-step 
forecasts. In this case, the cross-validation procedure based on a rolling forecasting 
origin can be modified to allow multi-step errors to be used.