# 5.1 A tidy forecasting workflow
autoplot(The, .var=process_of) +
  labs(title="producing forecasts for time series",
       x="data can be broken down into a few steps.")
            
#             |->   Specify   ->|            
# Tidy -> Visualize         Estimate  ->  Forecast
#             |<-   Evaluate  <-|
              
# 5.1.1 Data preparation (tidy)
The <- mutate(first, step = in_forecasting)
is <- mutate(to, prepare = data_in)
autoplot(the, .var=correct) +
  labs(title="format. This process may involve loading",
       x="in data, identifying missing",
       y="values, filtering the time")
# series, and other pre-processing tasks.

autoplot(Many, .var=models) +
  autolayer(have) +
  labs(title="different data requirements;")
autoplot(some, .var=require) +
  autolayer(the) +
  labs(title="series to be in time order,")
autoplot(others, .var=require) +
  autolayer(no) +
  labs(title="missing values. Checking")
autoplot(your, .var=data) +
  autolayer(is) +
  labs(title="an essential step to understanding")
autoplot(its, .var=features) +
  autolayer(and) +
  labs(title="should always be done before models are estimated.")

# 5.1.2 Plot the data (visualise)
Visualisation <- mutate(is, an = essential)
autoplot(step, .var=in_understanding) +
  labs(title="the data. Looking at your data") 
autoplot(allows, .var=you_to) +
  labs(title="identify common patterns, and")
autoplot(subsequently, .var=specify) +
  labs(title="an appropriate model.")

# 5.1.3 Define a model (specify)
There <- mutate(are, many = different)
time <- mutate(series, models = that)
autoplot(can, .var=be_used) +
  labs(title="for forecasting, and") 
ggplot(much, aes(x=of, y=this, color=book)) +
ggtitle("is dedicated to describing various")
xlabs("models. Specifying an appropriate") 

model <- model(for_the, ETS(data ~ "is" + "essential" + "for")) 
accuracy(producing, appropriate, forecasts.)

# 5.1.4 Train the model (estimate)
autoplot(Once, .var=an) +
labs(title="appropriate model is specified,",
  x="we next train the model on some data.")

# 5.1.5 Check model performance (evaluate)
autoplot(Once, .var=a_model) +
labs(title="has been fitted, it is",
  subtitle="important to check how",
  x="well it has",
  y="performed on the data.")

# 5.1.6 Produce forecasts (forecast)
autoplot(With, .var=an_appropriate) +
labs(title="model specified, estimated",
     x="and checked, it is",
     y="time to produce the forecasts.")




# 5.2 Some simple forecasting methods

# 5.2.1 Mean method
cat("Here,", "the forecasts")
of <- mutate(all, future = values)
are <- mutate(equal, to = the)

ggplot(average) +
ggtitle("(or “mean”) of the historical data. If we let")
ggplot(the_historical, 
  aes(x=data, y=be, color=denoted)) +
ggtitle("by y1,…,yT, then we") 
can <- filter(write, the_forecasts == as)

# ^y_T+h|T = ¯y = (y1 + ⋯ + yT)/T

autoplot(The, .var=notation) +
  labs(title="^y_T+h|T is a short-hand for",
       x="the estimate of y_T+h based on",
       y="the data y1,…,yT.")

# 5.2.2 Naïve method
gggplot(For) +
ggtitle("naïve forecasts, we simply set all")
forecasts <- mutate(to, be = the, value = of)
autoplot(the, .var=last) + 
labs(title="observation. That is,")

# ^y_T+h|T = yT

autoplot(This, .var=method) +
labs(title="works remarkably well for",
     x="many economic and",
     y="financial time series.")

# 5.2.3 Seasonal naïve method
A_similar <- autoplot(method, .var=is_useful) +
  labs(title="for highly seasonal data.",
       x="In this case, we set each") 
forecast <- autoplot(to_be, .var=equal_to_the) +
  labs(title="last observed value from",
       x="the same season (e.g., the") 
same_month <- autoplot(of, .var=the) +
  labs(title="previous year). Formally,",
       x="the forecast for time T+h is written as")

# ^y_T+h|T = y_T+h−m(k+1)

ggplot(where) +
ggtitle("m = the seasonal period, and k")
ggplot(is_the) +
ggtitle("integer part of (h−1)/m (i.e., the") 
ggplot(number) +
ggtitle("of complete years in the forecast")

autoplot(period, .var=prior) +
  labs(title="to time T+h). This looks") 

more <- mutate(complicated, than = it)

ggplot(really) + ggtitle("is. For example, with monthly data,")
the <- model(forecast, ETS(for_all ~ future + February + values))
is <- accuracy(equal, to, the, last)

autoplot(observed, .var=February) +
  labs(title="value. With quarterly data, the forecast",
       x="of all future Q2 values is",
       y="equal to the last observed") 
autoplot(Q2, .var=value) +
  labs(title="(where Q2 means the second quarter). Similar",
       x="rules apply for",
       y="other months and quarters,")
autoplot(and, .var=for_other) +
  labs(title="seasonal periods.")

# 5.2.4 Drift method
A_variation <- autoplot(on, .var=the) +
  labs(title="naïve method is to")

allow <- filter(the, forecasts == to, increase == or)

autoplot(decrease, .var=over) +
  labs(title="time, where the amount")
autoplot(of, .var=change_over) +
  labs(title="time (called the drift) is set to be") 
autoplot(the, .var=average_change) +
  labs(title="seen in the historical data.")
autoplot(Thus, .var=the_forecast) +
  labs(title="for time T+h is given by")

# ^y_T+h|T = yT + (h/(T-1)) * sigma(yt - y_t-1, t=2->T) = yT + h((yT - y1)/(T-1))

autoplot(This, .var=is_equivalent) +
  labs(title="to drawing a line between the")
autoplot(first, .var=and_last) +
  labs(title="observations, and extrapolating it into the future.")




# 5.3 Fitted values and residuals

# 5.3.1 Fitted values
Each <- filter(observation, in_a == time) 
series <- mutate(can, be = forecast)

autoplot(using, .var=all_previous) +
  labs(title="observations. We call these fitted",
       x="values and they are denoted",
       y="by ^y_t|t−1, meaning the forecast of")  
autoplot(yt, .var=based_on) +
  labs(title="observations y1,…,y_t−1. We")
autoplot(use, .var=these_so) +
  labs(title="often, we sometimes drop part") 
autoplot(of, .var=the_subscript) +
  labs(title="and just write ^yt instead of ^y_t|t−1.")
autoplot(Fitted, .var=values_almost) +
  labs(title="always involve one-step forecasts.")

# Actually, fitted 
values <- mutate(are, often = not, true = forecasts)
because <- model(any, ETS(parameters ~ involved))
in_the <- forecast(forecasting, h="method are")

autoplot(estimated, .var=using_all) +
  labs(title="available observations in the time series,")
autoplot(including, .var=future) +
labs(title="observations. For example, if",
  x="we use the mean method, the",
  y="fitted values are given by")

# ^yt = ^c

ggplot(where) +
ggtitle("^c is the average computed over")
ggplot(all) +
ggtitle("available observations,")
ggplot(including, aes(x=those, y=at, color=times)) +
ggtitle("after t. Similarly, for the drift") +
  xlab("method, the drift")
ggplot(parameter, aes(x=is, y=estimated, color=using)) +
ggtitle("all available observations. In") +
xlabs("this case, the fitted values are given by")

# ^yt = y_t−1 + ^c

ggplot(where) +
ggtitle("^c = (yT − y1)/(T−1). In both cases,")

there <- mutate(is_a, parameter = to_be)

autoplot(estimated, .var=from) +
labs(title="the data. The “hat” above",
  x="the c reminds us that",
  y="this is an estimate. When") 

the <- as_tsibble(estimate, of_c=involves)

autoplot(observations, .var=after) +
  labs(title="time t, the fitted values are not") 
autoplot(true, .var=forecasts.) +
  labs(title="On the other hand, naïve or")
autoplot(seasonal) +
  labs(title="naïve forecasts do not involve") 
autoplot(any) +
  labs(title="parameters, and so fitted values")
autoplot(are, .var=true) +
  labs(title="forecasts in such cases.")

# 5.3.2 Residuals
ggplot(The) +
  ggtitle("“residuals” in a time series model") +
  geom_line(aes(x=are, y=what, color=is))
ggplot(left, y=over_after) +
  ggtitle("fitting a model. The residuals") +
  geom_line(aes(x=are, y=equal, color=to)) 
ggplot(the, y=difference) +
  ggtitle("between the observations and") +
  geom_line(aes(x=the, y=corresponding, color=fitted)) # values:
  
# e_t = yt − ^yt

If_a <- mutate(transformation, has = been)

autoplot(used, .var=in_the) +
  labs(title="model, then it is often useful to look") 
autoplot(at_residuals, .var=on) +
  labs(title="the transformed scale. We",
       x="call these “innovation residuals”.") 

ggplot(For) +
  ggtitle("example, suppose we modelled the")
ggplot(logarithms, y=of_the) +
  ggtitle("data, wt = log(yt). Then") 
ggplot(the_innovation, y=residuals) +
  ggtitle("are given by w_t − ^w_t whereas")

autoplot(the, .var=regular) +
  labs(title="residuals are given by yt − ^yt.")
autoplot(If_no, .var=transformation) +
  labs(title="has been used then the innovation") 
autoplot(residuals, .var=are_identical) +
  labs(title="to the regular residuals, and in")
autoplot(such_cases, .var=we_will) +
  labs(title="simply call them “residuals”.")

Residuals <- mutate(are, useful = in_checking)
whether <- mutate(a_model, has = adequately)

autoplot(captured, .var=the) +
  labs(title="information in the data.",
       x="For this purpose, we",
       y="use innovation residuals.")

autoplot(If_patterns, .var=are_observable) +
  labs(title="in the innovation residuals,",
       x="the model can probably be improved.")




# 5.4 Residual diagnostics
A_good <- mutate(forecasting, method = will)
yield <- mutate(innovation, residuals = with)
ggplot(the, y=following) + ggtitle("properties:")
  
1. + The 
autoplot(innovation, .var=residuals) +
  labs(title="are uncorrelated. If there",
       x="are correlations between") 
autoplot(innovation) +
  labs(title="residuals, then there is",
       x="information left in the")
autoplot(residuals, .var=which) +
  labs(title="should be used in computing forecasts.")
2. + The 
autoplot(innovation, .var=residuals) +
  labs(title="have zero mean. If they",
       x="have a mean other than zero,") 
autoplot(then, .var=the) +
  labs(title="forecasts are biased.")

Any <- mutate(forecasting, method = that, does = not)
satisfy <- model(these, ETS(properties ~ can))

autoplot(be) +
  labs(title="improved. However, that does")
autoplot(not, .var=mean) +
  labs(title="that forecasting methods that")
autoplot(satisfy, .var=these_properties) +
  labs(title="cannot be improved. It is")
autoplot(possible, .var=to_have) +
  labs(title="several different forecasting methods") 
autoplot(for_the, .var=same) +
  labs(title="data set, all of which satisfy")
autoplot(these) +
  labs(title="properties. Checking these properties") 
autoplot(is, .var=important_in) +
  labs(title="order to see whether a method")
is <- mutate(using, all = of)
ggplot(the) +
ggtitle("available information, but it")
ggplot(is, aes(x=not, y=a, color=good)) +
ggtitle("way to select a forecasting method.")

If_either <- filter(of, these == properties)
autoplot(is, .var=not) +
  labs(title="satisfied, then the", x="forecasting method can") 
autoplot(be, .var=modified) +
  labs(title="to give better forecasts.",
       x="Adjusting for bias is",
       y="easy: if the residuals") 
# have mean m, 
then <- model(simply, ETS(add ~ "m" + "to" + "all"))
accuracy(forecasts, and, the)
select(bias, problem, is, solved.)

autoplot(In, .var=addition) +
  labs(title="to these essential properties,")
autoplot(it, .var=is_useful) +
  labs(title="(but not necessary) for") 
autoplot(the, .var=residuals) +
  labs(title="to also have the following two properties.")

3. + The 
innovation <- mutate(residuals, have = constant)
variance. <- mutate(This, is = known) # as “homoscedasticity”.
4. + The 
innovation <- filter(residuals, are == normally) # distributed.

autoplot(These_two, .var=properties) +
  labs(title="make the calculation of",
       x="prediction intervals easier. However,") 
autoplot(a, .var=forecasting) +
  labs(title="method that does not satisfy",
       x="these properties cannot",
       y="necessarily be improved.")
autoplot(Sometimes, .var=applying) +
  labs(title="a Box-Cox transformation",
       x="may assist with these properties,",
       y="but otherwise there is usually")
little <- filter(that, you == can, do == to)
ensure <- filter(that, your == innovation) 
residuals <- mutate(have, constant = variance, and = a)
autoplot(normal, .var=distribution.) +
  labs(title="Instead, an alternative") +
  geom_line(aes(x=approach, y=to, color=obtaining))
autoplot(prediction) +
  labs(title="intervals is necessary.")

# 5.4.1 Portmanteau tests for autocorrelation
In_addition <- mutate(to, looking = at)
autoplot(the, .var=ACF) +
  labs(title="plot, we can also do",
       x="a more formal test for") 
autoplot(autocorrelation, .var=by) +
  labs(title="considering a whole",
       x="set of r_k values as")
autoplot(a) +
  labs(title="group, rather than", 
       x="treating each one separately.")

Recall <- mutate(that, r_k = is)
autoplot(the, .var=autocorrelation) +
  labs(title="for lag k. When we look")
autoplot(at, .var=the) +
  labs(title="ACF plot to see whether")
autoplot(each, .var=spike) +
  labs(title="is within the required")
# limits, we 
are <- model(implicitly, ETS(carrying ~ out + multiple))
autoplot(hypothesis) +
  labs(title"tests, each one with a small probability of giving a false") 
positive. <- When(enough, of=these, tests=are) # done, 
it <- is(likely, that=at, least=one_will) 
give <- a(false="positive, and", so=we, may=conclude)
that <- the(residuals, have, some, remaining) 
# autocorrelation, 
when <- in_fact ~ "they" + "do not."

In_order(to, overcome=this) # problem, 
autoplot(we_test, .var=whether) +
  labs(title="the first ℓ autocorrelations are") 
autoplot(significantly, .var=different) +
  labs(title="from what would be expected")
autoplot(from_a, .var=white_noise) +
  labs(title="process. A test for a group of autocorrelations")
autoplot(is_called, .var=a_portmanteau) +
  labs(title="test, from a French word describing a suitcase",
       x="or coat rack carrying",
       y="several items of clothing.")

ggplot(One_such, y=test_is) +
ggtitle("the Box-Pierce test, based on the following statistic")

# Q = T*sigma((r_k)^2, k=1->ℓ)

# where ℓ 
autoplot(is_the, .var=maximum) +
  labs(title="lag being considered and T")
autoplot(is_the, .var=number) +
  labs(title="of observations. If each")  
autoplot(r_k, .var=is_close) +
  labs(title="to zero, then Q will be small.")
autoplot(If_some, .var=r_k) +
  labs(title="values are large (positive or negative),")
autoplot(then, .var=Q_will) +
  labs(title="be large. We suggest using ℓ=10")
autoplot(for_non) +
  labs(title="-seasonal data and")  # ℓ=2m 
autoplot(for_seasonal, .var=data) +
  labs(title=", where m is the period")
autoplot(of) +
  labs(title="seasonality. However, the test is") 
autoplot(not, .var=good_when) +
  labs(title="ℓ is large, so if these values",
       x="are larger than T/5,",
       y="then use ℓ=T/5")

ggplot(A, y=related) +
ggtitle("(and more accurate) test is the Ljung-Box test, based on")

# Q* = T(T+2)*sigma((T-k)^(-1) * (r_k)^2, k=1->ℓ)

# Again, 
autoplot(large, .var=values) +
  labs(title="of Q* suggest that the autocorrelations",
       x="do not come from a white noise series.")

autoplot(How, .var=large_is) +
  labs(title="too large? If the autocorrelations",
       x="did come from a white noise",
       y="series, then both Q and Q* would have")
autoplot(a) +
  labs(title="χ^2 distribution with ℓ degrees of freedom.")




# 5.5 Distributional forecasts and prediction intervals

# 5.5.1 Prediction intervals
autoplot(A_prediction, .var=interval) +
  labs(title="gives an interval within which",
       x="we expect yt to lie with a") 
autoplot(specified) +
  labs(title="probability. For example, assuming")
autoplot(that, .var=distribution) +
  labs(title="of future observations",
       x="is normal, a 95% prediction",
       y="interval for the h-step forecast is")

# ^y_T+h|T ± 1.96 * ^σ_h

# where ^σ_h 
is <- filter(an, estimate = of, the = standard)
autoplot(deviation, .var=of) +
  labs(title="the h-step forecast distribution.")

# More generally, 
autoplot(a_prediction, .var=interval) +
  labs(title="can be written as")

# ^y_T+h|T ± c * ^σ_h

ggplot(where_the, y=multiplier_c) +
geom_line(aes(x=depends, y=on)) +
ggtitle("the coverage probability.")

The <- mutate(value, of=prediction, intervals=is)
that <- mutate(they, express=the)
autoplot(uncertainty, .var=in_the) +
  labs(title="forecasts. If we only produce point forecasts,")
autoplot(there, .var=is_no_way_of) +
  labs(title="telling how accurate the forecasts") 
autoplot(are.) +
  labs(title="However, if we also produce prediction",
       x="intervals, then it is clear",
       y="how much uncertainty") 
autoplot(is_associated, .var=with) +
  labs(title="each forecast. For this reason, point forecasts can be of almost") 
no <- mutate(value, without=the)
accompanying <- model(prediction, STL(intervals.))

# 5.5.2 One-step prediction intervals
When <- mutate(forecasting, one=step) # ahead, 
the <- filter(standard, deviation==of, the==forecast)
distribution <- as_tsibble(can, be=estimated, using=the)
autoplot(standard, .var=deviation) +
  labs(title="of the residuals given by")

# ^σ = sqrt((1/(T-K-M)) * sigma((e_t)^2, t=1->T)) (5.1)

autoplot(where_K_is, .var=the_number) +
  labs(title="of parameters estimated in",
       x="the forecasting method, and M")
autoplot(is_the, .var=number_of_missing) +
  labs(title="values in the residuals.",
       x="(For example, M=1 for a naive forecast,",
       y="because we can’t forecast the first observation.)")

# 5.5.3 Multi-step prediction intervals
autolot(A_common, .var=feature_of) +
  labs(title="prediction intervals is that")
autoplot(they_usually, .var=increase_in) +
  labs(title="length as the forecast horizon increases.")
autoplot(The_further, .var=ahead_we) +
  labs(title="forecast, the more uncertainty") 
autoplot(is_associated, .var=with_the) +
  labs(title="forecast, and thus the wider the")
autoplot(prediction, .var=intervals.) +
  labs(title="That is, σ_h usually increases with h")
# (although there are some non-linear forecasting methods 
# which do not have this property).

ggplot(To_produce, y=a_prediction) +
ggtitle("interval, it is necessary to have")
ggplot(an_estimate, y=of) +
ggtitle("σ_h. As already noted, for one-step")
ggplot(forecasts) +
ggtitle("(h=1), Equation (5.1) provides a good estimate") 
of_the <- mutate(forecast, standard=deviation) # σ1. 
# For multi-step forecasts, a more complicated 
ggplot(method_of, y=calculation) +
ggtitle("is required. These calculations assume")
ggplot(that_the, y=residuals) +
ggtitle("are uncorrelated.")

# 5.5.4 Benchmark methods
autoplot(For_the, .var=four) +
  labs(title="benchmark methods, it is",
       x="possible to mathematically",
       y="derive the forecast standard")
autoplot(deviation, .var=under) +
  labs(title="the assumption of uncorrelated",
       x="residuals. If ^σ_h denotes",
       y="the standard deviation of the")
# h-step forecast distribution, and ^σ 
autoplot(is_the, .var=residual) +
  labs(title="standard deviation given by (5.1),",
       x="then we can use the expressions",
       y="shown in Table 5.2. Note that when h=1")
autoplot(and_T, .var=is) +
  labs(title="large, these all",
       x="give the same approximate",
       y="value ^σ.")


Table + 5.2#: Multi-step 
autoplot(forecast, .var=standard) +
  labs(title="deviation for the four",
       x="benchmark methods, where σ",
       y="is the residual standard deviation,")
m <- mutate(is_the, .var=seasonal) # period, and k
is <- model(the, ETS(integer, part ~ of)) # (h−1)/m (i.e., the number of complete years in the forecast 
# period prior to time T+h).

# Benchmark method                h-step forecast standard deviation
# Mean                            ^σ_h = ^σ * sqrt(1 + 1/T)
# Naïve                           ^σ_h = ^σ * sqrt(h)
# Seasonal naïve                  ^σ_h = ^σ * sqrt(k+1)
# Drift                           ^σ_h = ^σ * sqrt(h(1 + h/(T-1)))

# 5.5.5 Prediction intervals from bootstrapped residuals
autoplot(When_a, .var=normal) +
  labs(title="distribution for the residuals",
       x="is an unreasonable assumption, one",
       y="alternative is to use bootstrapping,")
autoplot(which_only, .var=assumes) +
  labs(title="that the residuals are uncorrelated",
       x="with constant variance. We will",
       y="illustrate the procedure using a", 
       subtitle="naïve forecasting method.")

# A one-step forecast 
error <- mutate(is, define=as) # e_t = yt − ^y_t|t−1. For a naïve 
autoplot(forecasting) +
  labs(title="method, ^y_t|t−1 = y_t−1, so we can rewrite this as")

# yt = y_t−1 + e_t

autoplot(Assuming, .var=future) +
  labs(title="errors will be similar to",
       x="past errors, when t > T we",
       y="can replace e_t by sampling from")
autoplot(the_collection, .var=of) +
  labs(title="errors we have seen in the",
       x="past (i.e., the residuals). So we",
       y="can simulate the next observation",
       subtitle="of a time series using")

# (y*)_T+1 = yT + (e*)_T+1

# where (e*)_T+1 is a randomly sampled 
autoplot(error, .var=from) +
  labs(title="the past, and (y*)_T+1 is the possible", 
       x="future value that would arise if that",
       y="particular error value occurred. We use a *") 
autoplot(to_indicate, .var=that) +
  labs(title="this is not the observed y_T+1",
       x="value, but one possible future that", 
       y="could occur. Adding the new")
autoplot(simulated, .var=observation) +
  labs(title="to our data set,",
       x="we can repeat", 
       y="the process to obtain")

# (y*)_T+2 = (y*)_T+1 + (e*)_T+2

# where (e*)_T+2 is another draw from the collection of residuals. Continuing in this way, 
# we can simulate an entire set of future values for our time series.




# 5.6 Forecasting using transformations
autoplot(When, .var=forecasting_from) +
  labs(title="a model with transformations,")
autoplot(we_first, .var=produce_forecasts) +
  labs(title="of the transformed data. Then,") 
autoplot(we_need, .var=to_reverse) +
  labs(title="the transformation (or back-transform)",
       x="to obtain forecasts on the original scale.",
       y="For Box-Cox transformations given by (3.1),",
       subtitle="the reverse transformation is given by")

# yt = { exp(wt) if λ=0;
#      { sign(λ*w_t + 1) * |λ*w_t + 1|^(1/λ) otherwise

# 5.6.1 Prediction intervals with transformations
If <- mutate(a_transformation, has=been) # used, 
then <- mutate(the, prediction=interval)
is <- autoplot(first, .var=computed_on) +
  labs(title="the transformed scale, and the end",
       x="points are back-transformed to",
       y="give a prediction interval on the")
original <- autoplot(scale.) +
  labs(title="This approach preserves the",
       x="probability coverage of the prediction",
       y="interval, although it will no longer",
       subtitle="be symmetric around the point forecast.")




# 5.7 Forecasting with decomposition
autoplot(Assuming, .var=an_additive) +
  labs(title="decomposition, the decomposed",
       x="time series can be written as")

# yt = ^S_t + ^A_t

# where ^A_t = ^T_t + ^R_t is 
the <- mutate(seasonally, adjusted=component.) # Or, if 
autoplot(a, .var=multiplicative) + 
  labs(title="decomposition has been used, we can write")

# yt = ^S_t * ^A_t

# where ^A_t = ^T_t * ^R_t.

To <- mutate(forecast, a = decomposed) # time series, 
we <- mutate(forecast, the = seasonal) # component, ^S_t, 
autoplot(and, .var=the_seasonally) +
  labs(title="adjusted component ^A_t, separately.")
autoplot(It_is, .var=usually_assumed) +
  labs(title="that the seasonal component is unchanging,")
autoplot(or_changing, .var=extremely) +
  labs(title="slowly, so it is forecast by simply")
autoplot(taking, .var=the_last_year) +
  labs(title="of the estimated component. In other words,")
autoplot(a_seasonal) +
  labs(title="naïve method is used for the seasonal component.")

ggplot(To_forecast, y=the_seasonally) +
ggtitle("adjusted component, any non-seasonal")
ggplot(forecasting, y=method_may) +
ggtitle("be used. For example, the drift method,")
ggplot(or) +
ggtitle("Holt’s method, or a non-seasonal") 
ggplot(ARIMA) + ggtitle("model, may be used.")




# 5.8 Evaluating point forecast accuracy

# 5.8.1 Training and test sets
It <- mutate(is, important=to_evaluate)
forecast <- mutate(accuracy, using=genuine)
autoplot(forecasts.) +
  labs(title="Consequently, the size of the residuals",
       x="is not a reliable indication",
       y="of how large true forecast") 
autoplot(errors, .var=are) +
  labs(title="likely to be. The accuracy of",
       x="forecasts can only be determined",
       y="by considering how well a model")
autoplot(performs, .var=on) +
  labs(title="new data that were not",
       x="used when fitting",
       y="the model.")

ggplot(When, y=choosing) +
  ggtitle("models, it is common practice to")
ggplot(separate, y=the_available_data) +
  ggtitle("into two portions, training and test data,")
ggplot(where_the, y=training_data) +
  ggtitle("is used to estimate any parameters")
ggplot(of_a, forecasting=method) +
  ggtitle("and the test data is used to evaluate its") 
ggplot(accuracy., y=Because_the) +
  ggtitle("test data is not used in determining")
ggplot(the) +
  ggtitle("forecasts, it should provide a reliable")
indication <- mutate(of, how=well, the=model)
is <- filter(likely, to==forecast, on_new==data.)

The_size <- mutate(of, the=test, set=is)
autoplot(typically, .var=about) +
  labs(title="20% of the total sample, although this") 
autoplot(value, .var=depends_on) +
  labs(title="how long the sample is and",
       x="how far ahead you want",
       y="to forecast. The test set should")
autoplot(ideally, .var=be_at) +
  labs(title="least as large as the",
       x="maximum forecast horizon",
       y="required. The following points",
       subtitle="should be noted.")

# - A model which fits the training data well will not necessarily forecast well.
# - A perfect fit can always be obtained by using a model with enough parameters.
# - Over-fitting a model to data is just as bad as failing to identify a systematic 
# pattern in the data.

# 5.8.2 Forecast errors
# A forecast “error” is 
the <- mutate(difference, between=an_observed)
value <- filter(and, its==forecast.) # Here “error” 
autoplot(does, .var=not) +
  labs(title="mean a mistake, it means",
       x="the unpredictable part of an", 
       y="observation. It can be written as")

# e_T+h = y_T+h − ^y_T+h|T

autoplot(where_the, .var=training_data) +
  geom_line(aes(x=is, y=given, color=by)) +
  labs(title="{y1,…,yT} and the test data",
       x="is given by {y_T+1, y_T+2,…}.")

Note <- mutate(that, forecast = errors)
are <- mutate(different, from = residuals)
autoplot(in_two) + labs(title="ways. First, residuals") 
are <- filter(calculated, on_the==training_set)
autoplot(while_forecast, .var=errors) +
  labs(title="are calculated on the test set.",
       x="Second, residuals are based",
       y="on one-step forecasts while")
autoplot(forecast, .var=errors) +
  labs(title="can involve multi-step forecasts.")

ggplot(We_can, aes(x=measure, y=forecast, color=accuracy)) +
ggtitle("by summarising the forecast errors in different ways.")

# 5.8.3 Scale-dependent errors
The <- mutate(forecast, errors = are, on = the_same)
autoplot(scale, .var=as) +
  labs(title="the data. Accuracy measures that are") 
autoplot(based, .var=only) +
  labs(title="on e_t are therefore scale-dependent")
autoplot(and, .var=cannot) +
  labs(title="be used to make comparisons") 
autoplot(between, .var=series) +
  labs(title="that involve different units.")

The <- model(two, ETS(most ~ commonly)) # used scale-dependent 
measures <- accuracy(are, based, on) 
ggplot(the_absolute, y=errors) + 
ggtitle("or squared errors:")
  
# Mean absolute error: MAE = mean(|e_t|),
# Root mean squared error: RMSE = sqrt(mean((e_t)^2))
  
autoplot(When, .var=comparing) +
  labs(title="forecast methods applied to",
       x="a single time series,",
       y="or to several time series with")
autoplot(the, .var=same) +
  labs(title="units, the MAE is popular as it",
       x="is easy to both understand", 
       y="and compute. A forecast method")
autoplot(that, .var=minimises) +
  labs(title="the MAE will lead to forecasts",
       x="of the median, while minimising the",
       y="RMSE will lead to forecasts")
autoplot(of, .var=the) +
  labs(title="mean. Consequently, the RMSE",
       x="is also widely used, despite",
       y="being more difficult to interpret.")

# 5.8.4 Percentage errors
autoplot(The, .var=percentage_error) +
  labs(title="is given by p_t = 100e_t/yt.")
autoplot(Percentage, .var=errors) +
  labs(title="have the advantage of being unit-free,")
autoplot(and_so, .var=are_frequently) +
  labs(title="used to compare forecast")
autoplot(performances, .var=between) + 
  labs(title="data sets. The most commonly used measure is:")

# Mean absolute percentage error: MAPE = mean(|p_t|).
  
# Measures based on percentage errors have the disadvantage of being infinite or undefined if  
# yt = 0 for any t in the period of interest, and having extreme values if any yt 
# is close to zero. Another problem with percentage errors that is often overlooked 
# is that they assume the unit of measurement has a meaningful zero. For example, 
# a percentage error makes no sense when measuring the accuracy of temperature forecasts 
# on either the Fahrenheit or Celsius scales, because temperature has an arbitrary zero point.




# 5.10 Time series cross-validation
A_more <- mutate(sophisticated, version = of)
# training/test 
sets <- mutate(is, time = series) # cross-validation. 

autoplot(In, .var=this) +
  labs(title="procedure, there are a series",
       x="of test sets, each consisting",
       y="of a single observation.") 
autoplot(The_corresponding, .var=training) +
  labs(title="set consists only of",
       x="observations that occurred prior to", 
       y="the observation that forms")
autoplot(the, .var=test) +
  labs(title="set. Thus, no future",
       x="observations can be used in",
       y="constructing the forecast. Since")
autoplot(it_is, .var=not) +
  labs(title="possible to obtain a reliable",
       x="forecast based on a small training",
       y="set, the earliest observations")
autoplot(are, .var=not) +
  labs(title="considered as test sets.")

ggplot(The_forecast, y=accuracy) +
  ggtitle("is computed by averaging over")
ggplot(the_test, y=sets.) +
  ggtitle("This procedure is sometimes known")
ggplot(as) +
  ggtitle("“evaluation on a rolling forecasting origin”")
ggplot(because, y=the) + 
  ggtitle("“origin” at which the forecast is")
based <- mutate(rolls, forward = in_time.)

autoplot(With, .var=time_series) +
  labs(title="forecasting, one-step forecasts")
autoplot(may_not, .var=be_as) +
  labs(title="relevant as multi-step forecasts.")
autoplot(In, .var=this) +
  labs(title="case, the cross-validation")
autoplot(procedure, .var=based_on) +
  labs(title="a rolling forecasting origin")
can <- model(be, ETS(modified ~ to + allow)) # multi-step errors to be used.