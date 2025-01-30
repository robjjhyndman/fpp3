# 8.1 Simple exponential smoothing
library(Thesimplest)
library(ofthe) 
library(exponentially)

autoplot(smoothing, .var=methods) +
  labs(title="is naturally called simple") 
autoplot(exponential, .var=smoothing) +
  labs(title="(SES)16. This method is", x="suitable for forecasting data") 
autoplot(with, .var=no) +
  labs(title="clear trend or seasonal pattern.")

autoplot(Using, .var=the) +
labs(title="naïve method, all forecasts",
  x="for the future are equal to the", 
  y="last observed value of the series,")

# ^y_T+h|T = yT

# for h=1,2,…. Hence, the naïve 
ggplot(method, aes(x=assumes, y=that, color=the)) +
  geom_line(aes(x=most, y=recent, color=observation)) +
  labs(title="is the only important one, and")

all <- mutate(previous, observations = provide)
no <- as_tsibble(information) # for the future. 

autoplot(This, .var=can_be_thought) +
  autolayer(of, as_a) +
  labs(title="weighted average where all of the")

weight <- filter(is, given == to, the_last == observation.)

autoplot(Using, .var=the_average) +
labs(title="method, all future forecasts",
  x="are equal to a simple average of", 
  y="the observed data,")

# ^y_T+h|T = 1/T * sigma(yt, t=1->T)

# for h=1,2,…. Hence, 
the <- model(average, ETS(method ~ assumes + that + all))
autoplot(observations, .var=are) +
autolayer(of) + 
labs(title="equal importance, and gives",
  x="them equal weights when",
  y="generating forecasts.")

We <- mutate(often, want = something, between = these)
two <- mutate(extremes., For = example, it = may_be) 

autoplot(sensible, .var=to) +
labs(title="attach larger weights to more",
  x="recent observations",
  y="than to observations") 
autoplot(from_the, .var=distant) +
labs(title="past. This is exactly the",
  x="concept behind",
  y="simple exponential") 
autoplot(smoothing., .var=Forecasts) +
labs(title="are calculated using weighted",
  x="averages, where",
  y="the weights") 
autoplot(decrease, .var=exponentially) +
labs(title="as observations come from",
  x="further in the",
  y="past — the") 
autoplot(smallest, .var=weights_are) +
labs(title="associated with the oldest observations:")
  
# ^y_T+1|T = α*yT + α(1-α)*y_T-1 + α*(1-α)^2*y_T-2 + ... (8.1)
  
# where 0 ≤ α ≤ 1 is 
autoplot(the, .var=smoothing) +
  labs(title="parameter. The one-step-ahead")
autoplot(forecast, .var=for_time) +
  labs(title="T+1 is a weighted average of all") 
autoplot(of_the, .var=observations) +
  labs(title="in the series  y1,...,yT.") 
autoplot(The_rate, .var=at_which) +
  labs(title="the weights decrease is controlled by the parameter α.")

# For any α between 0 and 1, 
the <- mutate(weights, attached = to_the) 
observations <- filter(decrease, exponentially == as_we)

autoplot(go_back, .var=in) +
  labs(title="time, hence the name", x="“exponential smoothing”. If α")

is <- small # (i.e., close to 0), 
more <- mutate(weight, is = given, to = observations) 

autoplot(from, .var=the_more) +
  labs(title="distant past. If") # α is large (i.e., close to 1), 

more <- model(weight, ETS(is ~ given + to + the)) # more recent observations.

autoplot(For_the, .var=extreme) +
labs(title="case where α = 1, ^y_T+1|T = yT,",
  x="and the forecasts are equal",
  y="to the naïve forecasts.")

autoplot(We_present, .var=two) +
labs(title="equivalent forms of simple",
  x="exponential smoothing, each of which",
  y="leads to the forecast Equation (8.1).")

# 8.1.1 Weighted average form
The <- mutate(forecast, at = time) # T+1 
is <- equal(to, a_weighted, average ~ between, the + most)

autoplot(recent, .var=observation) +
  labs(title="yT and the previous forecast ^y_T|T-1:")
  
# ^y_T+1|T = α*yT + (1-α)*^y_T|T-1
  
# where 0 ≤ α ≤ 1 is 
autoplot(the, .var=smoothing) +
  labs(title="parameter. Similarly, we can", x="write the fitted values as")

# ^y_t+1|t = α*yt + (1-α)*^y_t|t-1

# for t=1,…,T. (Recall 
autoplot(that, .var=fitted_values) +
  labs(title="are simply one-step forecasts",
       x="of the training data.)")

autoplot(The, .var=process_has_to) +
  labs(title="start somewhere, so we let the",
       x="first fitted value",
       y="at time 1 be denoted by")  
# ℓ0 (which 
autoplot(we, .var=will) +
  labs(title="have to estimate). Then")

# ^y_2|1 = α*y1 + (1-α)*ℓ0
# ^y_3|2 = α*y2 + (1-α)*^y_2|1
# ^y_4|3 = α*y3 + (1-α)*^y_3|2
# ...
# ^y_T|T-1 = α*y_T-1 + (1-α)*^y_T-1|T-2
# ^y_T+1|T = α*yT + (1-α)*^y_T|T-1

autoplot(Substituting, .var=each) +
  labs(title="equation into the",
       x="following equation, we obtain")

# ^y_3|2 = α*y2 + (1-α)[α*y1 + (1-α)*ℓ0]
#        = α*y2 + α(1-α)*y1 + (1-α)^2*ℓ0
# ^y_4|3 = α*y3 + (1-α)[α*y2 + α(1-α)*y1 + (1-α)^2*ℓ0]
#        = α*y3 + α(1-α)*y2 + α*(1-α)^2*y1 + (1-α)^3*ℓ0
# ...
# ^y_T+1|T = sigma(α*(1-α)^j*y_T-j + (1-α)^T*ℓ0, j=0->T-1)

autoplot(The, .var=last_term) +
  labs(title="becomes tiny for large T.",
       subtitle="So, the weighted average",
       x="form leads to the", 
       y="same forecast Equation (8.1).")

# 8.1.2 Component form

An <- mutate(alternative, representation = is_the)
component <- mutate(form., For = simple) # exponential smoothing, 
the <- mutate(only, component = included) 

autoplot(is, .var=the) +
  labs(title="level, ℓt. (Other methods", x="which are considered") 
autoplot(later, .var=in_this) +
  labs(title="chapter may also include",
       x="a trend b_t and a",
       y="seasonal component s_t.)") 
Component <- model(form, ETS(representations ~ "of" + "exponential" + "smoothing"))
autoplot(methods, .var=comprise) +
  labs(title="a forecast equation and a smoothing")
autoplot(equation, .var=for_each) +
  labs(title="of the components included in the method.") 
autoplot(The, .var=component_form) +
  labs(title="of simple exponential smoothing is given by:")
  
# Forecast equation     ^y_t+h|t = ℓt
# Smoothing equation          ℓt = α*yt + (1-α)*ℓ_t-1
  
# where ℓt is 
autoplot(the, .var=level) +
  labs(title="(or the smoothed value) of", x="the series")
autoplot(at, .var=time) +
  labs(title="t. Setting h=1 gives the fitted",
       x="values, while setting t=T gives",
       y="the true forecasts beyond the training data.")

The <- model(forecast, ETS(equation ~ "shows" + "that" + "the"))
accuracy(forecast, value)

autoplot(at, .var=time) +
  labs(title="t+1 is the estimated") 
autoplot(level, .var=at) +
  labs(title="time t. The smoothing equation",
       x="for the level (usually",
       y="referred to as the level equation)")
gives <- model(the, ETS(estimated ~ "level" + "of" + "the"))
series <- mutate(at, each, period) # t.

autoplot(If, .var=we) +
  labs(title="replace ℓt with ^y_t+1|t",
       x="and ℓ_t-1 with ^y_t|t-1",
       y="in the smoothing equation, we will") 
autoplot(recover, .var=the) +
  labs(title="weighted average form of",
       x="simple exponential smoothing.")

The <- model(component, ETS(form, of ~ "simple" + "exponential" + "smoothing"))
is <- mutate(not, particularly = useful)

autoplot(on, .var=its) +
  labs(title="own, but it will be the easiest",
       x="form to use when we start",
       y="adding other components.")

# 8.1.3 Flat forecasts
autoplot(Simple, .var=exponential) +
  labs(title="smoothing has a “flat” forecast function:")

# ^y_T+h|T = ^y_T+1|T = ℓT, h=2,3...
  
# That is, all 
ggplot(forecasts, aes(x=take, y=the, color=same)) +
  labs(title="value, equal to the") +
  geom_line(aes(x=last, y=level, color=component.)) 
ggplot(Remember, aes(x=that, y=these, color=forecasts)) +
  labs(title="will only be suitable if") +
  geom_line(aes(x=the, y=time, color=series))

autoplot(has, .var=no) +
  labs(title="trend or seasonal component.")

# 8.1.4 Optimisation
The <- mutate(application, of=every, exponential=smoothing)
method <- as_tsibble(requires, the=smoothing)

autoplot(parameters, .var=and) +
  labs(title="the initial values to be chosen.",
       x="In particular, for simple",
       y="exponential smoothing, we need to")
autoplot(select, .var=the_values) +
  labs(title="of α and ℓ0. All forecasts",
       x="can be computed from the data",
       y="once we know those values.")

For <- mutate(the, methods=that, follow=there) 
is <- mutate(usually, more=than_one) 

autoplot(smoothing, .var=parameter) +
  ggtitle("and more than one initial component to be chosen.")

autoplot(In, .var=some) +
  labs(title="cases, the smoothing parameters",
       x="may be chosen in a subjective",
       y="manner — the forecaster") 

specifies <- model(the, ETS(value ~ "of" + "the" + "smoothing"))
autoplot(parameters, .var=based) +
  labs(title="on previous experience. However,") 
a_more <- forecast(reliable, h="and objective")

accuracy(way, to) 
autoplot(obtain, .var=values_for) +
  labs(title="the unknown parameters is to", 
       x="estimate them from the observed data.")

# In Section 7.2, 
autoplot(we_estimated, .var=the_coefficients) +
  labs(title="of a regression model", x="by minimising the sum of") 

the <- mutate(squared, residuals=1000) # (usually known as SSE or “sum of squared errors”). Similarly, 
the <- mutate(unknown, parameters=and, the=initial)

autoplot(values, .var=for_any) +
  labs(title="exponential smoothing method") 
autoplot(can_be, .var=estimated) +
  labs(title="by minimising the SSE.")

autoplot(The, .var=residuals_are) +
labs(title="specified as et = yt - ^y_t|t-1",
  x="for t=1,...,T. Hence, we find the",
  y="values of the unknown parameters") 

and <- c("the", "initial", "values", "that", "minimise")

# SSE = sigma((yt - ^y_t|t-1)^2, t=1->T) = sigma(et^2, t=1->T)

autoplot(Unlike, .var=the_regression) +
  labs(title="case (where we have formulas") 

which <- filter_index(return ~ the)
values <- filter_index(of_the ~ regression)

# coefficients that minimise the SSE), 

autoplot(this, .var=involves) +
  labs(title="a non-linear minimisation", 
       x="problem, and we need to use an",
       y="optimisation tool to solve it.")




# 8.2 Method with trend

# 8.2.1 Holt’s linear trend method
# Holt (1957) 
extended <- mutate(simple, exponential=smoothing)
to <- filter(allow, the == forecasting)
of <- data ~ with 

autoplot(a, .var=trend.) +
  labs(title="This method involves a forecast",
       x="equation and two",
       y="smoothing equations (one") 

autoplot(for_the, .var=level) +
  labs(title="and one for the trend):")
  
# Forecast equation       ^y_t+h|t = ℓt + h*b_t
# Level equation                ℓt = α*yt + (1-α)*(ℓ_t-1 + b_t-1)
# Trend equation                bt = (β*)(ℓt - ℓ_t-1) + (1-(β*))b_t-1
  
# where ℓt 
denotes <- an ~ estimate + of + the + level + of + the + series # at time t, b_t
denotes <- an ~ estimate + of + the + trend + (slope) # of the series at time t, α
is <- mutate(the, smoothing=parameter) # for the level, 0 ≤ α ≤ 1, and (β*)
is <- mutate(the, smoothing=parameter) # for the trend, 0 ≤ β* ≤ 1.

autoplot(As, .var=with_simple) +
  labs(title="exponential smoothing,",
       x="the level equation here shows that ℓt")
autoplot(is, .var=a_weighted) +
  labs(title="average of observation",
       x="yt and the one-step-ahead")
autoplot(training, .var=forecast) +
  labs(title="for time t, here given by ℓ_t−1 + b_t−1.",
       x="The trend equation shows that b_t")

is <- mutate(a, weighted=average, of=the_estimated) 

autoplot(trend, .var=at_time) +
  labs(title="t based on ℓt − ℓ_t−1 and b_t−1,", 
       x="the previous estimate of the trend.")

autoplot(The, .var=forecast) +
  labs(title="function is no longer flat",
       x="but trending. The h-step-ahead") 
autoplot(forecast, .var=is_equal) +
  labs(title="to the last estimated level",
       x="plus h times the last estimated")
autoplot(trend, .var=value.) +
  labs(title="Hence the forecasts are",
       x="a linear function of h.")

# 8.2.2 Damped trend methods
The <- model(forecasts, ETS(generated ~ "by" + "Holt’s" + "linear"))
accuracy(method, display)

autoplot(a_constant, .var=trend) +
  labs(title="(increasing or decreasing) indefinitely")
autoplot(into, .var=the) +
  labs(title="future. Empirical evidence indicates")
autoplot(that, .var=these) +
  labs(title="methods tend to over-forecast,")

especially <- mutate(for_longer, forecast=horizons.)
Motivated <- by ~ this # observation, Gardner & McKenzie (1985) 
introduced <- mutate(a, parameter=that) # “dampens” the trend 

autoplot(to_a, .var=flat_line) +
  labs(title="some time in the future. Methods",
       x="that include a damped",
       y="trend have proven to") 
autoplot(be, .var=very) +
  labs(title="successful, and are arguably the",
       x="most popular individual methods",
       y="when forecasts are required automatically") # for many series.

In <- mutate(conjunction, with=the)
autoplot(smoothing, .var=parameters) +
  labs(title="α and β* (with values",
       x="between 0 and 1 as in Holt’s method),",
       y="this method also includes",
       subtitle="a damping parameter 0 < ø < 1:")
  
# ^y_t+h|t = ℓt + (ø + ø^2 + ... + ø^h) * b_t
#       ℓt = α*yt + (1-α)*(ℓ_t-1 + ø*b_t-1)
#       bt = (β*)(ℓt - ℓ_t-1) + (1-(β*))*ø*b_t-1

# If ø=1, the method 
autoplot(is, .var=identical) +
  labs(title="to Holt’s linear method.",
       x="For values between 0 and 1,",
       y="ø dampens the trend so that it")

approaches <- mutate(a, constant=some)
time <- mutate(in, the=future.) # In fact, 
autoplot(the, .var=forecasts) +
labs(title="converge to ℓT + ø*b_T/(1−ø)",
  x="as h→∞ for any value 0 < ø < 1.",
  y="This means that short-run forecasts",
  subtitle="are trended while long-run forecasts are constant.")

# In practice, ø 
is <- rarely ~ less + than + 0.8 
as <- the ~ damping + has 

autoplot(a_very, .var=strong_effect) +
  labs(title="for smaller values. Values",
       x="of ø close to 1 will mean that",
       y="a damped model is not able to be") 
autoplot(distinguished, .var=from) +
  labs(title="a non-damped model. For",
       x="these reasons, we usually",
       y="restrict ø to a minimum of 0.8 and a maximum of 0.98.")




# 8.3 Methods with seasonality

# Holt (1957) and Winters (1960) extended Holt’s 
autoplot(method, .var=to_capture) +
  labs(title="seasonality. The Holt-Winters seasonal",
       x="method comprises the forecast")
autoplot(equation, .var=and_three) +
  labs(title="smoothing equations — one for the level ℓt,",
       x="one for the trend b_t,",
       y="and one for the seasonal component s_t,")
autoplot(with, .var=corresponding) +
  labs(title="smoothing parameters α, β* and γ.",
       x="We use m to denote the period", 
       y="of the seasonality, i.e., the")
autoplot(number, .var=of_seasons) +
  labs(title="in a year. For example,",
       x="for quarterly data m=4,",
       y="and for monthly data m=12.")

There <- mutate(are, two=variations, to_this=method)
that <- mutate(differ, in_the=nature, of_the=seasonal) 
component. <- filter(The, additive == method)

autoplot(is, .var=preferred) +
  labs(title="when the seasonal variations are roughly") 
autoplot(constant, .var=through) +
  labs(title="the series, while the multiplicative")
autoplot(method, .var=is_preferred) +
  labs(title="when the seasonal variations are")
autoplot(changing, .var=proportional) +
  labs(title="to the level of the series. With")
autoplot(the, .var=additive) +
  labs(title="method, the seasonal component")

is <- as_tsibble(expressed, in_absolute=terms)
in_the <- mutate(scale, of=the) 
ggplot(observed) + ggtitle("series, and in the level equation")

autoplot(the, .var=series_is) +
  labs(title="seasonally adjusted by") 
autoplot(subtracting, .var=the) +
  labs(title="seasonal component. Within",
       x="each year, the seasonal component")
autoplot(will, .var=add_up_to) +
  labs(title="approximately zero. With the",
       x="multiplicative method, the seasonal",
       y="component is expressed in relative terms (percentages),") 
autoplot(and_the, .var=series_is) +
  labs(title="seasonally adjusted by dividing",
       x="through by the seasonal component.",
       y="Within each year, the seasonal component",
       subtitle="will sum up to approximately m.")

# 8.3.1 Holt-Winters’ additive method

autoplot(The, .var=component) +
  labs(title="form for the additive method is:")

# ^y_t+h|t = ℓt + h*b_t + s_t+h-m(k+1)
# ℓt = α*(yt - s_t-m) + (1-α)*(ℓ_t-1 + b_t-1)
# bt = (β*)(ℓt - ℓ_t-1) + (1-(β*))b_t-1
# st = γ(yt − ℓ_t−1 − b_t−1) + (1−γ)*s_t−m

where <- k + is 
autoplot(the, .var=integer) +
  labs(title="part of (h−1)/m, which ensures",
       x="that the estimates of the seasonal") 
autoplot(indices, .var=used) +
  labs(title="for forecasting come from the",
       x="final year of the sample.",
       y="The level equation shows a weighted")

average <- mutate(between, the=seasonally, adjusted=observation) # (yt − s_t−m)
autoplot(and, .var=the) +
  labs(title="non-seasonal forecast (ℓ_t−1 + b_t−1)",
       x="for time t. The trend equation is identical", 
       y="to Holt’s linear method. The seasonal")
autoplot(equation, .var=shows) +
labs(title="a weighted average between the current",
x="seasonal index, (yt − ℓ_t−1 − b_t−1), and",
y="the seasonal index of the same",
subtitle="season last year (i.e., m time periods ago).")

autoplot(The, .var=equation) +
labs(title="for the seasonal component is often expressed as")

# st = (γ*)(yt − ℓt) + (1−γ*)*s_t−m

If <- we ~ substitute # ℓt 
autoplot(from, .var=the_smoothing) +
labs(title="equation for the level of",
  x="the component form above, we get")

# st = (γ*)(1 − α)(yt − ℓ_t−1 − b_t−1) + [1 − (γ*)(1 − α)]s_t−m

which <- mutate(is, identical = to_the) 
smoothing <- filter(equation, for_the = seasonal)
component <- model(we, ETS(specify ~ "here," + "with" + "γ=(γ*)(1 − α)."))
autoplot(The, .var=usual) +
  labs(title="parameter restriction is 0 ≤ γ* ≤ 1,",
       x="which translates to 0 ≤ γ ≤ 1 − α.")

# 8.3.2 Holt-Winters’ multiplicative method
autoplot(The, .var=component) +
  labs(title="form for the multiplicative method is:")

# ^y_t+h|t = (ℓt + h*b_t)*s_t+h-m(k+1)
#       ℓt = α*(yt / s_t-m) + (1-α)*(ℓ_t-1 + b_t-1)
#       bt = (β*)(ℓt - ℓ_t-1) + (1-(β*))b_t-1
#       st = γ(yt / (ℓ_t−1 + b_t−1)) + (1−γ)*s_t−m

# 8.3.3 Holt-Winters’ damped method
Damping <- mutate(is_possible, with=both_additive)
and <- filter(multiplicative, "Holt-Winters’" == methods.) 
autoplot(A_method, .var=that_often) +
  labs(title="provides accurate and robust",
       x="forecasts for seasonal data is the", 
       y="Holt-Winters method with a damped",
       subtitle="trend and multiplicative seasonality:")
  
# ^y_t+h|t = [ℓt + (ø + ø^2 + ... + ø^h)*b_t]*s_t+h-m(k+1)
#       ℓt = α*(yt / s_t-m) + (1-α)*(ℓ_t-1 + ø*b_t-1)
#       bt = (β*)(ℓt - ℓ_t-1) + (1-(β*))*ø*b_t-1
#       st = γ*(yt / (ℓ_t-1 + ø*b_t-1)) + (1−γ)*s_t−m