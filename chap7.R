# # 7.1 The linear model

us_change_long <- pivot_longer(us_change,
                               c(Consumption, Income, Production, Savings, Unemployment),
                               names_to="Series")
autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
# # 7.1.1 Simple linear regression
# In the simplest case, the regression model allows for a linear relationship
# between the forecast variable y and a single predictor variable x:
# # yt = β0 + β1*xt + εt

us_change_long <- pivot_longer(us_change,
                               c(Consumption, Income, Production, Savings, Unemployment),
                               names_to="Series")
autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
# An artificial example of data from such a model is shown in Figure 7.1. The coefficients  
# β0 and β1 denote the intercept and the slope of the line respectively. The intercept  
# β0 represents the predicted value of y when x=0. The slope β1 represents the 
# average predicted change in y resulting from a one unit increase in x.
# # draw the graph

us_change_long <- pivot_longer(us_change,
                               c(Consumption, Income, Production, Savings, Unemployment),
                               names_to="Series")
autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
# Notice that the observations do not lie on the straight line but are scattered around 
# it. We can think of each observation yt as consisting of the systematic or 
# explained part of the model, β0 + β1*xt, and the random “error”, εt. The “error” term 
# does not imply a mistake, but a deviation from the underlying straight line model. It 
# captures anything that may affect yt other than xt.

us_change_long <- pivot_longer(us_change,
                               c(Consumption, Income, Production, Savings, Unemployment),
                               names_to="Series")
autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
# # 7.1.2 Multiple linear regression
# When there are two or more predictor variables, the model is called a multiple 
# regression model. The general form of a multiple regression model is 
# # yt = β0 + β1*x_1,t + β2*x_2,t + ⋯ + βk*x_k,t + εt (7.1)

us_change_long <- pivot_longer(us_change,
                               c(Consumption, Income, Production, Savings, Unemployment),
                               names_to="Series")
autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where y is the variable to be forecast and x1,…,xk are the k predictor variables. 
# Each of the predictor variables must be numerical. The coefficients β1,…,βk measure the 
# effect of each predictor after taking into account the effects of all the other 
# predictors in the model. Thus, the coefficients measure the marginal effects of 
# the predictor variables.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.1.3 Assumptions
# When we use a linear regression model, we are implicitly making some assumptions about 
# the variables in Equation (7.1).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# First, we assume that the model is a reasonable approximation to reality; that is, the 
# relationship between the forecast variable and the predictor variables satisfies this 
# linear equation.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Second, we make the following assumptions about the errors (ε1,…,εT):

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# - they have mean zero; otherwise the forecasts will be systematically biased.
# - they are not autocorrelated; otherwise the forecasts will be inefficient, as there is 
# more information in the data that can be exploited.
# - they are unrelated to the predictor variables; otherwise there would be more information 
# that should be included in the systematic part of the model.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# It is also useful to have the errors being normally distributed with a constant variance σ^2
# in order to easily produce prediction intervals.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Another important assumption in the linear regression model is that each predictor x
# is not a random variable. If we were performing a controlled experiment in a laboratory, we 
# could control the values of each x (so they would not be random) and observe the resulting 
# values of y. With observational data (including most data in business and economics), it is
# not possible to control the value of x, we simply observe it. Hence we make this an assumption.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.2 Least squares estimation
# In practice, of course, we have a collection of observations but we do not know the values 
# of the coefficients β0, β1, …, βk. These need to be estimated from the data.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The least squares principle provides a way of choosing the coefficients effectively by 
# minimising the sum of the squared errors. That is, we choose the values of β0, β1, …, βk
# that minimise

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # sigma(εt^2, t=1->T) = sigma((yt - β0 - β1*x_1,t - β2*x_2,t - ... - βk*x_k,t)^2, t=1->T)
# 
# This is called least squares estimation because it gives the least value for the sum of 
# squared errors. Finding the best estimates of the coefficients is often called “fitting” 
# the model to the data, or sometimes “learning” or “training” the model.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.2.1 Fitted values
# Predictions of y can be obtained by using the estimated coefficients in the regression 
# equation and setting the error term to zero. In general we write,
# 
# # ^yt = ^β0 + ^β1*x_1,t + ^β2*x_2,t + ⋯ + ^βk*x_k,t

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Plugging in the values of x_1,t,…,x_k,t for t=1,…,T returns predictions of yt within 
# the training set, referred to as fitted values. Note that these are predictions of the 
# data used to estimate the model, not genuine forecasts of future values of y

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.2.2 Goodness-of-fit
# A common way to summarise how well a linear regression model fits the data is via the 
# coefficient of determination, or R^2. This can be calculated as the square of the 
# correlation between the observed y values and the predicted ^y values. Alternatively, 
# it can also be calculated as,

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # R^2 = sigma((^yt - yhat)^2) / sigma((yt - yhat)^2)
# 
# where the summations are over all observations. Thus, it reflects the proportion of 
# variation in the forecast variable that is accounted for (or explained) by the regression model.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# In simple linear regression, the value of R^2 is also equal to the square of the correlation 
# between y and x (provided an intercept has been included).
# 
# If the predictions are close to the actual values, we would expect R^2 to be close to 1. 
# On the other hand, if the predictions are unrelated to the actual values, then R^2=0 (again, 
# assuming there is an intercept). In all cases, R^2 lies between 0 and 1.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The R^2 value is used frequently, though often incorrectly, in forecasting. The value of  
# R^2 will never decrease when adding an extra predictor to the model and this can lead 
# to over-fitting. There are no set rules for what is a good R^2  value, and typical values of  
# R^2 depend on the type of data used. Validating a model’s forecasting performance on the 
# test data is much better than measuring the R^2 value on the training data.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.2.3 Standard error of the regression
# Another measure of how well the model has fitted the data is the standard deviation of 
# the residuals, which is often known as the “residual standard error”. This is shown in the 
# above output with the value 0.31.
# 
# # ^σ_e = sqrt((1 / T-k-1)*sigma(et^2, t=1->T))

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where k is the number of predictors in the model. Notice that we divide by T−k−1
# because we have estimated k+1 parameters (the intercept and a coefficient for each predictor 
# variable) in computing the residuals.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The standard error is related to the size of the average error that the model produces. We 
# can compare this error to the sample mean of y or with the standard deviation of y
# to gain some perspective on the accuracy of the model.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.3 Evaluating the regression model
# The differences between the observed y values and the corresponding fitted ^y values are the 
# training-set errors or “residuals” defined as,
# 
# # et = yt − ^yt
# #    = yt − ^β0 − ^β1*x_1,t − ^β2*x_2,t − ⋯ − ^βk*x_k,t

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# for t=1,…,T. Each residual is the unpredictable component of the associated observation.
# 
# The residuals have some useful properties including the following two:
# 
# # sigma(et, t=1->T) = 0 and sigma(x_k,t*et, t=1->T) = 0 for all k
#   
# As a result of these properties, it is clear that the average of the residuals is zero, 
# and that the correlation between the residuals and the observations for the predictor variable 
# is also zero. (This is not necessarily true when the intercept is omitted from the model.)

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# After selecting the regression variables and fitting a regression model, it is necessary to 
# plot the residuals to check that the assumptions of the model have been satisfied. There are 
# a series of plots that should be produced in order to check different aspects of the fitted 
# model and the underlying assumptions. We will now discuss each of them in turn.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.3.1 ACF plot of residuals
# With time series data, it is highly likely that the value of a variable observed in the current 
# time period will be similar to its value in the previous period, or even the period before that,
# and so on. Therefore when fitting a regression model to time series data, it is common to find 
# autocorrelation in the residuals. In this case, the estimated model violates the assumption of 
# no autocorrelation in the errors, and our forecasts may be inefficient — there is some information
# left over which should be accounted for in the model in order to obtain better forecasts. The 
# forecasts from a model with autocorrelated errors are still unbiased, and so they are not “wrong”, 
# but they will usually have larger prediction intervals than they need to. Therefore we should 
# always look at an ACF plot of the residuals.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.3.2 Outliers and influential observations
# Observations that take extreme values compared to the majority of the data are called 
# outliers. Observations that have a large influence on the estimated coefficients of a regression 
# model are called influential observations. Usually, influential observations are also outliers 
# that are extreme in the x direction.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# One source of outliers is incorrect data entry. Simple descriptive statistics of your data can 
# identify minima and maxima that are not sensible. If such an observation is identified, and it 
# has been recorded incorrectly, it should be corrected or removed from the sample immediately.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Outliers also occur when some observations are simply different. In this case it may not be wise 
# for these observations to be removed. If an observation has been identified as a likely outlier, 
# it is important to study it and analyse the possible reasons behind it. The decision to remove or 
# retain an observation can be a challenging one (especially when outliers are influential observations). 
# It is wise to report results both with and without the removal of such observations.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.3.3 Spurious regression
# More often than not, time series data are “non-stationary”; that is, the values of the time series 
# do not fluctuate around a constant mean or with a constant variance.
# 
# Regressing non-stationary time series can lead to spurious regressions. High  
# R^2 and high residual autocorrelation can be signs of spurious regression.
# 
# Cases of spurious regression might appear to give reasonable short-term forecasts, but they will 
# generally not continue to work into the future.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.4 Some useful predictors
# 
# # 7.4.1 Trend
# It is common for time series data to be trending. A linear trend can be modelled by simply using  
# x_1,t = t as a predictor,
# 
# # yt = β0 + β1*t + εt
# where t=1,…,T.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.4.2 Dummy variables
# So far, we have assumed that each predictor takes numerical values. But what about when 
# a predictor is a categorical variable taking only two values (e.g., “yes” and “no”)? Such 
# a variable might arise, for example, when forecasting daily sales and you want to take 
# account of whether the day is a public holiday or not. So the predictor takes value “yes” 
# on a public holiday, and “no” otherwise.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# This situation can still be handled within the framework of multiple regression models by 
# creating a “dummy variable” which takes value 1 corresponding to “yes” and 0 corresponding 
# to “no”. A dummy variable is also known as an “indicator variable”.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# A dummy variable can also be used to account for an outlier in the data. Rather than omit 
# the outlier, a dummy variable removes its effect. In this case, the dummy variable takes 
# value 1 for that observation and 0 everywhere else. An example is the case where a special 
# event has occurred. For example when forecasting tourist arrivals to Brazil, we will need to 
# account for the effect of the Rio de Janeiro summer Olympics in 2016.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# If there are more than two categories, then the variable can be coded using several dummy 
# variables (one fewer than the total number of categories). 

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.4.3 Seasonal dummy variables
# Suppose that we are forecasting daily data and we want to account for the day of the week as 
# a predictor. Then the following dummy variables can be created.
# #       d_1,t   d_2,t   d_3,t   d_4,t   d_5,t   d_6,t
# # Mon     1       0       0       0       0       0
# # Tue     0       1       0       0       0       0
# # Wed     0       0       1       0       0       0
# # Thu     0       0       0       1       0       0
# # Fri     0       0       0       0       1       0
# # Sat     0       0       0       0       0       1
# # Sun     0       0       0       0       0       0
# # Mon     1       0       0       0       0       0
# # ...    ...     ...     ...     ...     ...     ...

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Notice that only six dummy variables are needed to code seven categories. That is because the 
# seventh category (in this case Sunday) is captured by the intercept, and is specified when 
# the dummy variables are all set to zero.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Many beginners will try to add a seventh dummy variable for the seventh category. This is 
# known as the “dummy variable trap”, because it will cause the regression to fail. There will 
# be one too many parameters to estimate when an intercept is also included. The general rule 
# is to use one fewer dummy variables than categories. So for quarterly data, use three dummy 
# variables; for monthly data, use 11 dummy variables; and for daily data, use six dummy variables, 
# and so on.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The interpretation of each of the coefficients associated with the dummy variables is that it 
# is a measure of the effect of that category relative to the omitted category. In the above example, 
# the coefficient of d_1,t associated with Monday will measure the effect of Monday on the 
# forecast variable compared to the effect of Sunday. An example of interpreting estimated dummy 
# variable coefficients capturing the quarterly seasonality of Australian beer production follows.





autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.5 Selecting predictors
# 
# # 7.5.1 Adjusted R^2
# Computer output for a regression will always give the R^2 value, discussed in 
# Section 7.2. However, it is not a good measure of the predictive ability of a model. 
# It measures how well the model fits the historical data, but not how well the model will forecast future data.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# In addition, R^2 does not allow for “degrees of freedom”. Adding any variable tends 
# to increase the value of R^2, even if that variable is irrelevant. For these reasons, forecasters 
# should not use R^2 to determine whether a model will give good predictions, as it will lead 
# to overfitting.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# An equivalent idea is to select the model which gives the minimum sum of squared errors 
# (SSE), given by
# 
# # SSE = sigma(et^2, t=1->T)
# 
# Minimising the SSE is equivalent to maximising R^2 and will always choose the model with 
# the most variables, and so is not a valid way of selecting predictors.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# An alternative which is designed to overcome these problems is the adjusted R^2
# (also called “R-bar-squared”):
#   
# # ¯R^2 = 1 − (1 − R^2) * (T−1)/(T−k−1)

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where T is the number of observations and k is the number of predictors. This is an improvement on  
# R^2, as it will no longer increase with each added predictor. Using this measure, the 
# best model will be the one with the largest value of ¯R^2. Maximising ¯R^2 is equivalent to 
# minimising the standard error ^σ_e given in Equation (7.3).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Maximising ¯R^2 works quite well as a method of selecting predictors, although it does tend 
# to err on the side of selecting too many predictors.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.5.2 Cross-validation
# Time series cross-validation was introduced in Section 5.10 as a general tool for determining 
# the predictive ability of a model. For regression models, it is also possible to use 
# classical leave-one-out cross-validation to select predictors (Bergmeir et al., 2018). 
# This is faster and makes more efficient use of the data. The procedure uses the following steps:

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# 1. Remove observation t from the data set, and fit the model using the remaining data. Then 
# compute the error ((e*)_t = yt − ^yt) for the omitted observation. (This is not the same 
# as the residual because the t-th observation was not used in estimating the value of ^yt.)
# 2. Repeat step 1 for t=1,…,T.
# 3. Compute the MSE from (e*)_1,…,(e*)_T. We shall call this the CV.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Although this looks like a time-consuming procedure, there are fast methods of calculating 
# CV, so that it takes no longer than fitting one model to the full data set. The equation for 
# computing CV efficiently is given in Section 7.9. Under this criterion, the best model is the 
# one with the smallest value of CV.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.5.3 Akaike’s Information Criterion
# A closely-related method is Akaike’s Information Criterion, which we define as
# 
# # AIC = T*log(SSE/T) + 2(k+2)
# 
# where T is the number of observations used for estimation and k is the number of predictors 
# in the model. Different computer packages use slightly different definitions for the AIC, 
# although they should all lead to the same model being selected. The k+2 part of the equation 
# occurs because there are k+2 parameters in the model: the k coefficients for the predictors, 
# the intercept and the variance of the residuals. The idea here is to penalise the fit of the 
# model (SSE) with the number of parameters that need to be estimated.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The model with the minimum value of the AIC is often the best model for forecasting. For large 
# values of T, minimising the AIC is equivalent to minimising the CV value.
# 
# # 7.5.4 Corrected Akaike’s Information Criterion
# For small values of T, the AIC tends to select too many predictors, and so a bias-corrected 
# version of the AIC has been developed,

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # AIC_c = AIC + (2(k+2)(k+3))/(T-k-3)
# 
# As with the AIC, the AICc should be minimised.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.5.5 Schwarz’s Bayesian Information Criterion
# A related measure is Schwarz’s Bayesian Information Criterion (usually abbreviated to BIC, SBIC or SC):
#   
# # BIC = T*log(SSE/T) + (k+2)*log(T)
#   
# As with the AIC, minimising the BIC is intended to give the best model. The model chosen by 
# the BIC is either the same as that chosen by the AIC, or one with fewer terms. This is because 
# the BIC penalises the number of parameters more heavily than the AIC. For large values of T, 
# minimising BIC is similar to leave-v-out cross-validation when v=T[1−1/(log(T)−1)].

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.5.6 Which measure should we use?
# While ¯R^2 is widely used, and has been around longer than the other measures, its 
# tendency to select too many predictor variables makes it less suitable for forecasting.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Many statisticians like to use the BIC because it has the feature that if there is a true 
# underlying model, the BIC will select that model given enough data. However, in reality, 
# there is rarely, if ever, a true underlying model, and even if there was a true underlying 
# model, selecting that model will not necessarily give the best forecasts (because the 
# parameter estimates may not be accurate).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Consequently, we recommend that one of the AICc, AIC, or CV statistics be used, each of which 
# has forecasting as their objective. If the value of T is large enough, they will all lead to 
# the same model. In most of the examples in this book, we use the AIC_c value to select the 
# forecasting model.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.6 Forecasting with regression
# 
# # 7.6.1 Building a predictive regression model
# The great advantage of regression models is that they can be used to capture important 
# relationships between the forecast variable of interest and the predictor variables. However, 
# for ex ante forecasts, these models require future values of each predictor, which can be 
# challenging. If forecasting each predictor is too difficult, we may use scenario-based forecasting 
# instead, where we assume specific future values for all predictors.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# An alternative formulation is to use as predictors their lagged values. Assuming that we are 
# interested in generating a h-step ahead forecast we write
# 
# # y_t+h = β0 + β1*x_1,t +⋯+ βk*x_k,t + ε_t+h

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# for h=1,2…. The predictor set is formed by values of the xs that are observed h time periods 
# prior to observing y. Therefore when the estimated model is projected into the future, i.e., 
# beyond the end of the sample T, all predictor values are available.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Including lagged values of the predictors does not only make the model operational for easily 
# generating forecasts, it also makes it intuitively appealing. For example, the effect of a policy 
# change with the aim of increasing production may not have an instantaneous effect on consumption 
# expenditure. It is most likely that this will happen with a lagging effect.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.6.2 Prediction intervals
# As this involves some advanced matrix algebra we present here the case for calculating prediction intervals 
# for a simple regression, where a forecast can be generated using the equation,
# 
# # ^y = ^β0 + ^β1*x

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Assuming that the regression errors are normally distributed, an approximate 95% prediction 
# interval associated with this forecast is given by
# 
# # ^y ± 1.96*(^σ)_e*sqrt(1 + 1/T + (x-¯x)^2/((T-1)*(s_x)^2))  (7.4)

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where T is the total number of observations, ¯x is the mean of the observed x values, s_x
# is the standard deviation of the observed x values and ^σ_e is the standard error of the 
# regression given by Equation (7.3). Similarly, an 80% prediction interval can be obtained by 
# replacing 1.96 by 1.28. Other prediction intervals can be obtained by replacing the 1.96 with 
# the appropriate value given in Table 5.1. If the fable package is used to obtain prediction 
# intervals, more exact calculations are obtained (especially for small values of T) than what 
# is given by Equation (7.4).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Equation (7.4) shows that the prediction interval is wider when x is far from ¯x. That is, we 
# are more certain about our forecasts when considering values of the predictor variable close 
# to its sample mean.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.7 Nonlinear regression
# Although the linear relationship assumed so far in this chapter is often adequate, there are 
# many cases in which a nonlinear functional form is more suitable. To keep things simple in this 
# section we assume that we only have one predictor x.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# The simplest way of modelling a nonlinear relationship is to transform the forecast variable y
# and/or the predictor variable x before estimating a regression model. While this provides 
# a non-linear functional form, the model is still linear in the parameters. The most commonly
# used transformation is the (natural) logarithm (see Section 3.1).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# A log-log functional form is specified as
# 
# # logy = β0 + β1*logx + ε
# 
# In this model, the slope β1 can be interpreted as an elasticity: β1 is the average 
# percentage change in y resulting from a 1% increase in x. Other useful forms can also 
# be specified. The log-linear form is specified by only transforming the forecast 
# variable and the linear-log form is obtained by transforming the predictor.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Recall that in order to perform a logarithmic transformation to a variable, all of 
# its observed values must be greater than zero. In the case that variable x
# contains zeros, we use the transformation log(x+1); i.e., we add one to the value 
# of the variable and then take logarithms. This has a similar effect to taking 
# logarithms but avoids the problem of zeros. It also has the neat side-effect of zeros 
# on the original scale remaining zeros on the transformed scale.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# There are cases for which simply transforming the data will not be adequate and a 
# more general specification may be required. Then the model we use is
# 
# # y = f(x) + ε

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where f is a nonlinear function. In standard (linear) regression, f(x) = β0 + β1*x. In 
# the specification of nonlinear regression that follows, we allow f to be a more 
# flexible nonlinear function of x, compared to simply a logarithmic or other transformation.
# 
# One of the simplest specifications is to make f piecewise linear. That is, we 
# introduce points where the slope of f can change. These points are called knots. 
# This can be achieved by letting x1 = x and introducing variable x2 such that

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # x2 = (x−c)_+ = { 0 if x < c
# #                { x-c if x ≥ c
# 
# The notation (x−c)_+ means the value x−c if it is positive and 0 otherwise. This 
# forces the slope to bend at point c. Additional bends can be included in the relationship 
# by adding further variables of the above form.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Piecewise linear relationships constructed in this way are a special case of regression splines. 
# In general, a linear regression spline is obtained using
# 
# # x1 = x, x2 = (x − c1)_+, …, xk = (x − c_k−1)_+

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# where c1,…,c_k−1 are the knots (the points at which the line can bend). Selecting the number 
# of knots (k−1) and where they should be positioned can be difficult and somewhat arbitrary. 
# Some automatic knot selection algorithms are available, but are not widely used.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.7.1 Forecasting with a nonlinear trend
# In Section 7.4 fitting a linear trend to a time series by setting x=t was introduced. The 
# simplest way of fitting a nonlinear trend is using quadratic or higher order trends 
# obtained by specifying
# 
# # x_1,t = t, x_2,t = t^2, …

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# However, it is not recommended that quadratic or higher order trends be used in forecasting. 
# When they are extrapolated, the resulting forecasts are often unrealistic.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# A better approach is to use the piecewise specification introduced above and fit a piecewise 
# linear trend which bends at some point in time. We can think of this as a nonlinear trend 
# constructed of linear pieces. If the trend bends at time τ, then it can be specified by 
# simply replacing x=t and c=τ above such that we include the predictors,
# 
# # x_1,t = t
# # x_2,t = (t − τ)_+ = { 0 if t < τ
# #                   = { t-τ if t ≥ τ

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# in the model. If the associated coefficients of x_1,t and x_2,t are β1 and β2, then  
# β1 gives the slope of the trend before time τ, while the slope of the line after time  
# τ is given by β1 + β2. Additional bends can be included in the relationship by adding 
# further variables of the form (t−τ)_+ where τ is the “knot” or point in time at which 
# the line should bend.




autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.8 Correlation, causation and forecasting
# 
# # 7.8.1 Correlation is not causation
# It is important not to confuse correlation with causation, or causation with forecasting. A variable  
# x may be useful for forecasting a variable y, but that does not mean x is causing y. It is possible 
# that x is causing y, but it may be that y is causing x, or that the relationship between them is more 
# complicated than simple causality.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# For example, it is possible to model the number of drownings at a beach resort each month 
# with the number of ice-creams sold in the same period. The model can give reasonable 
# forecasts, not because ice-creams cause drownings, but because people eat more ice-creams
# on hot days when they are also more likely to go swimming. So the two variables (ice-cream
# sales and drownings) are correlated, but one is not causing the other. They are both caused 
# by a third variable (temperature). This is an example of “confounding” — where an omitted
# variable causes changes in both the response variable and at least one predictor variable.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# We describe a variable that is not included in our forecasting model as a confounder when
# it influences both the response variable and at least one predictor variable. Confounding 
# makes it difficult to determine what variables are causing changes in other variables, but 
# it does not necessarily make forecasting more difficult.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Similarly, it is possible to forecast if it will rain in the afternoon by observing the 
# number of cyclists on the road in the morning. When there are fewer cyclists than usual, 
# it is more likely to rain later in the day. The model can give reasonable forecasts, not
# because cyclists prevent rain, but because people are more likely to cycle when the published 
# weather forecast is for a dry day. In this case, there is a causal relationship, but 
# in the opposite direction to our forecasting model. The number of cyclists falls because
# there is rain forecast. That is, y (rainfall) is affecting x (cyclists).

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# It is important to understand that correlations are useful for forecasting, even when 
# there is no causal relationship between the two variables, or when the causality runs 
# in the opposite direction to the model, or when there is confounding.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# However, often a better model is possible if a causal mechanism can be determined. A 
# better model for drownings will probably include temperatures and visitor numbers and 
# exclude ice-cream sales. A good forecasting model for rainfall will not include cyclists, 
# but it will include atmospheric observations from the previous few days.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.8.2 Forecasting with correlated predictors
# When two or more predictors are highly correlated it is always challenging to accurately
# separate their individual effects. Suppose we are forecasting monthly sales of a company 
# for 2012, using data from 2000–2011. In January 2008, a new competitor came into the market 
# and started taking some market share. At the same time, the economy began to decline. In 
# your forecasting model, you include both competitor activity (measured using advertising 
# time on a local television station) and the health of the economy (measured using GDP). 
# It will not be possible to separate the effects of these two predictors because they 
# are highly correlated.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Having correlated predictors is not really a problem for forecasting, as we can still 
# compute forecasts without needing to separate out the effects of the predictors. However, 
# it becomes a problem with scenario forecasting as the scenarios should take account of 
# the relationships between predictors. It is also a problem if some historical analysis 
# of the contributions of various predictors is required.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# # 7.8.3 Multicollinearity and forecasting
# A closely related issue is multicollinearity, which occurs when similar information is 
# provided by two or more of the predictor variables in a multiple regression.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# It can occur when two predictors are highly correlated with each other (that is, they
# have a correlation coefficient close to +1 or -1). In this case, knowing the value of
# one of the variables tells you a lot about the value of the other variable. Hence, they
# are providing similar information. For example, foot size can be used to predict height,
# but including the size of both left and right feet in the same model is not going to 
# make the forecasts any better, although it won’t make them worse either.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Multicollinearity can also occur when a linear combination of predictors is highly 
# correlated with another linear combination of predictors. In this case, knowing the 
# value of the first group of predictors tells you a lot about the value of the second 
# group of predictors. Hence, they are providing similar information.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# An example of this problem is the dummy variable trap discussed in Section 7.4. Suppose 
# you have quarterly data and use four dummy variables, d1, d2, d3 and d4. Then 
# d4 = 1 − d1 − d2 − d3, so there is perfect correlation between d4 and d1 + d2 + d3.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# In the case of perfect correlation (i.e., a correlation of +1 or -1, such as in the
# dummy variable trap), it is not possible to estimate the regression model.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# If there is high correlation (close to but not equal to +1 or -1), then the estimation 
# of the regression coefficients is computationally difficult. In fact, some software 
# (notably Microsoft Excel) may give highly inaccurate estimates of the coefficients. 
# Most reputable statistical software will use algorithms to limit the effect of 
# multicollinearity on the coefficient estimates, but you do need to be careful. The major 
# software packages such as R, SPSS, SAS and Stata all use estimation algorithms to avoid
# the problem as much as possible.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# When multicollinearity is present, the uncertainty associated with individual regression
# coefficients will be large. This is because they are difficult to estimate. Consequently,
# statistical tests (e.g., t-tests) on regression coefficients are unreliable. (In 
# forecasting we are rarely interested in such tests.) Also, it will not be possible to 
# make accurate statements about the contribution of each separate predictor to the forecast.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Forecasts will be unreliable if the values of the future predictors are outside the 
# range of the historical values of the predictors. For example, suppose you have fitted 
# a regression model with predictors x1 and x2 which are highly correlated with each other, 
# and suppose that the values of x1 in the training data ranged between 0 and 100. Then 
# forecasts based on x1 > 100 or x1 < 0 will be unreliable. It is always a little dangerous 
# when future values of the predictors lie much outside the historical range, but it is 
# especially problematic when multicollinearity is present.

autoplot(filter(us_change_long, Series=="Consumption" | Series=="Income"), value) +
  labs(y = "% change")
ggplot(us_change, aes(x = Quarter)) +
  geom_line(aes(y = Consumption, colour = "Consumption")) +
  geom_line(aes(y = Income, colour = "Income")) +
  labs(x = "Quarter", y = "% change") +
  guides(colour = guide_legend(title = NULL))
# Note that if you are using good statistical software, if you are not interested in the
# specific contributions of each predictor, and if the future values of your predictor 
# variables are within their historical ranges, there is nothing to worry about — multicollinearity 
# is not a problem except when there is perfect correlation.
