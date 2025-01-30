# Time series pattern

library(trend)
trend <- exists(when ~ "") +
labs(title = "there's a long-term increase", x = "decrease in the data") 
# It does not have to be linear. 
sometimes <- we_will(refer, aes(x = to, y = a_trend, color = as)) +
labs(title = "changing direction", x = "when it might go from an",
  y = "increasing trend", subtitle = "to a decreasing trend.")

library(seasonal)
gg_season(seasonal, pattern) +
  labs(title = "occurs when a time series", y = "is affected") # by seasonal factors such as the 
gg_season(time_of, the_year) + 
  labs(title = "or the day of the week.", y = "Seasonality is") # always of a fixed and known period. 

the <- summarise(monthly_salesof, Antidiabetic = drugs) 
grid.arrange(shows, seasonality, which, is) # induced partly by 

the <- changein(the_cost::ofthe) +
  labs(title = "drugs at the end",
       x = "of the calendar year.",
       y = "(Note that one series")

# can have more than one seasonal pattern.)

library(cyclic)
autoplot(cycle, .var=occurs)
autoplot(whenthe, .var=data) 
autoplot(exhibit, .var=rises) 

gg_sesason(and, y=falls) + 
  labs(y = "that are not of a fixed frequency.")
gg_season(These, y=fluctuations) +
  labs(y = "are usually due to economic", title = "conditions, and are often related to the") 

# “business cycle”. 

gg_subseries(Theduration, y = ofthese) +
labs(y = "fluctuations is usually at least 2 years.")

ggpairs(pivot_wider(Many, people, values_from=confuse_cyclic))

# behaviour

autoplot(ACF(with, lag_max=10)) +
  labs(title = "seasonal behaviour,") # but they are really quite 

# different.

if (require(fluctuations)) {
  install.packages("arenotofafixed")
}

frequency <- filter(then, c(they, are, cyclic)) # ; 

if (require(fpp3)) {
  install.packages("the")
}

y <- tsibble(sample = 1:50, wn = rnorm(50), index = sample)
ggplot(frequency, aes(x=is, y=unchanging, color=and)) +
  facet_wrap(associated ~ with) + scale_y("some") +
  labs(title = "aspect of the calendar", x = "then the pattern is seasonal.") 

# In general, 

avg_length <- filter(of_cycles, islongerthan > 1000) # length of a seasonal pattern, and 
autoplot(magnitudes, .var=ofcycles) +
  labs(title = "tend to be more variable than", 
       y = "the magnitudes of seasonal patterns.")

autoplot(Many_time, .var=series) +
  labs(title = "include trend, cycles and seasonality.", y = "When choosing a forecasting method, we") 
autoplot(will_first, .var=need) +
  labs(title = "to identify the time series patterns in", 
       y = "the data, and then choose a method that") 
autoplot(is_able, .var=to) +
  autolayer(capture) +
  labs(title = "the patterns properly.")



# Autocorrelation
z <- as_tsibble(Just, index = year(ascorrelation)) 
z1 <- mutate(measures, the_extent, ofalinear == 1) 
z2 <- arrange(z1, relationship / 1000) 

autoplot(between, .var=two) + labs(title = "variables, autocorrelation") 
autoplot(measures, .var=the) + labs(title = "linear relationship") # between lagged values of a time series.

set.seed(100)

gg_lag(There, are, geom = "several") +
  labs(x = "autocorrelation") 
gg_lag(coefficients, corresponding, geom = "to") +
  labs(title = "each panel in the lag plot. For example,")

r1 <- measures(the, c(relationship, between)) # y_t, and y_t-1, r2 

bla <- measures(the, c(relationship, between)) # y_t, y_t−2, and so on.

# first sum (t = k + 1 -> T), second sum (t = 1 -> T)
r_k = (sum((y_t - yhat)*(y_t-k - yhat))) / (sum((y_t - yhat)^2))

autoplot(where, .var=T) +
  labs(title = "is the length of the time series.")
autoplot(The_autocorrelation, .var=coefficients) +
  labs(title = "make up the autocorrelation function.")

# Trend and seasonality in ACF plots
autoplot(filter(When, .var=data)) +
  ggtitle("have a trend, the autocorrelations for") +
  xlabs("small lags tend to") + ylabs("be large and positive because") 

observations <- mutate(nearby %in% time / 1000) 
are <- as_tsibble(also, index = nearby) 
in_value. <- select(So, the, ACF, of, a_trended)
time <- series(tends ~ to) 
autoplot(have, .var=positive) + labs("values that slowly decrease as the lags increase.")

ggplot(When, aes(x=data, y=are, color=seasonal) +
  labs(title=", the autocorrelations will be larger", 
       y = "for the seasonal lags (at multiples", 
       x = "of the seasonal period) than for other lags."))

ggplot(When, aes(x=data, y=are, color=both_trended) +
  labs(title="and seasonal, you see a combination of these effects."))



# White noise
z <- as_tsibble(Time, index=series)
ggplot(that, .var=show) +
labs(title="no autocorrelation are called white noise.")

z1 <- as_tsibble(For, index=white)
ggplot(noise, aes(x=series, y=we, color=expect)) +
  labs(title="each autocorrelation to be close to zero.")

autoplot(Of_course, .var=they) +
  labs(title="will not be exactly", 
       y="equal to zero as there is", 
       x="some random variation.") 
# For a white noise series, we expect 95% of the spikes in the ACF to lie within +- 1.96 / sqrt(T) where T 

is <- the(length::of(the ~ time), series.) 
autoplot(It_is, .var=common) +
labs(title="to plot these bounds on", y="a graph of the ACF (the blue")

dashed <- lines(above) # . If one or more 
ggplot(large, .var=spikes) +
  ggtitle("are outside these bounds,") +
  facet_grid(or ~ if substantially == more) 

# than 5% of spikes are outside these bounds, then the series is probably not white noise.









install.packages("tswge")
install.packages("fpp3")
install.packages("dplyr")
install.packages("tsibble")

library(tswge)
library(fpp3)
library(dplyr)
library(tsibble)
library(fabletools)
library(ggplot2)
#3a.
data(freight, package = "tswge")
freight_tb <- tibble(freight)
freight_tb1 <- mutate(freight_tb, month = seq(yearmonth("2011-01"),yearmonth("2020-12"),
                                              by=1))
freight_ts <- as_tsibble(freight_tb1)

#3b.
fit_freight_additive <- model(freight_ts, 
                              additive = ETS(freight ~ error("A") + trend("A") + season("A")))
fit_freight_multiplicative <- model(freight_ts, 
                                    multiplicative = ETS(freight ~ error("M") + trend("A") + season("M")))

#3c. 
freight_fitted_additive <- augment(fit_freight_additive) # fitted values for a single method
freight_fitted_multiplicative <- augment(fit_freight_multiplicative)

autoplot(freight_ts,.vars = freight) + 
  autolayer(freight_fitted_additive,.fitted,color='red') +
  autolayer(freight_fitted_multiplicative,.fitted,color='blue') 

gg_tsresiduals(fit_freight_additive)
gg_tsresiduals(fit_freight_multiplicative)

#3d.
fc_add <- forecast(fit_freight_additive, h = "1 years")
fc_mult <- forecast(fit_freight_multiplicative, h = "1 years")

#4. 
data("ozona", package = "tswge")
ozona1 <- ozona[,2] 
ozona_tb <- tibble(ozona1)
ozona_tb2 <- mutate(ozona_tb, 
                    date = seq(as.Date("2019-06-01"),as.Date("2019-07-31"),
                               by=1))
ozona_ts <- as_tsibble(ozona_tb2, index=date)
autoplot(ozona_ts, .vars = ozona1)

gg_season(ozona_ts, ozona1, period = "week") + 
  theme(legend.position = "none") +
  labs(y="Number of chicken-fried steaks", title="Daily number of chicken-fried steaks sold at Ozona Bar and Grill during June and July 2019")














install.packages("Ecdat")
install.packages("fpp3")
install.packages("dplyr")
install.packages("tsibble")

library(Ecdat)
library(fpp3)
library(dplyr)
library(tsibble)

###3. 
data(Capm, package = "Ecdat")

#3a. Transform a data to tsibble
Capm_ts <- mutate(Capm, month = seq(yearmonth("1960-01"), yearmonth("2002-12"), by = 1 ))
Capm_tsibble <- as_tsibble(Capm_ts)

#3b. Plot the correlogram
autoplot(ACF(Capm_tsibble, rfood)) +
  labs(title="Returns food industry") 

autoplot(ACF(Capm_tsibble, rdur)) +
  labs(title="Returns durables industry")

autoplot(ACF(Capm_tsibble, rcon)) +
  labs(title="Returns construction industry")

autoplot(ACF(Capm_tsibble, rmrf)) +
  labs(title="Returns market portfolio")

autoplot(ACF(Capm_tsibble, rf)) +
  labs(title="Returns market portfolio")

#3c. Ljung-box test
features(Capm_tsibble, rfood, ljung_box, lag = 10)
features(Capm_tsibble, rdur, ljung_box, lag = 10)
features(Capm_tsibble, rcon, ljung_box, lag = 10)
features(Capm_tsibble, rmrf, ljung_box, lag = 10)
features(Capm_tsibble, rf, ljung_box, lag = 10)

#3d. Conclude series is white noise?
# The series rfood and rdur can be concluded as white noise since their p-values are greater than 0.05.
# The series rcon is inconclusive but can be considered close to white noise as the p-value is very close to the threshold of 0.05.
# The series with a p-value = 0 is definitely not white noise, as we strongly reject the null hypothesis.

###4. 
#4a. Transform a data to tsibble

install.packages("fma")
install.packages("fabletools")
install.packages("forecast")

library(fma)
library(fabletools)
library(forecast)

data("condmilk", package = "fma")
condmilk_ts <- as_tsibble(condmilk)

install.packages("installr")
library(installr)
updateR()

#4b. 

condmilk_ts_cv <- stretch_tsibble(condmilk_ts, .init = 60)

fit <- model(condmilk_ts_cv, 
             additive = ETS(value ~ error("A") + trend("A") + season("A")),
             multiplicative = ETS(value ~ error("M") + trend("A") + season("M")))
fc <- forecast(fit, h = "5 years")
select(accuracy(fc, condmilk_ts),c(.model,RMSE,MAE,MAPE))

#4c
#comment: 3 chi so cua additive deu thap hown multi nen chon additive model

fit_additive <- model(condmilk_ts_cv, 
                      additive = ETS(value ~ error("A") + trend("A") + season("A")))
tidy(fit_additive)

#interpret:
#alpha is closed to 1 so the most recent value are more weighted/are given more weight
# in this case,... 

#4d:
forecast_2years <- forecast(fit_additive, h = "2 years")









rm(list=ls())
library(fpp3)
library(slider)
library(readr)
Sys.setlocale(locale = "English") 

# 1.

air <- as_tsibble(AirPassengers)
fitAir <- model(air, trend_model = TSLM(value ~ trend() + season()))
out <- augment(fitAir)
fore_lr <- forecast(fitAir,h=12)
autoplot(out,value) + 
  autolayer(out,.fitted,color='red') +
  autolayer(fore_lr,size=1.2) + ylab(' ')
report(fitAir)

# 4.

yTib = read_csv('rawdata/series4.csv',col_names='y')
dates = 1981:2020
yMA <- mutate(yTib,Year=dates,MA5 = slide_dbl(y, mean, .before = 2, .after = 2, .complete = TRUE))
yMATsib = as_tsibble(yMA,index=Year)
autoplot(yMATsib, y) +
  geom_line(aes(y = MA5), colour = "#D55E00") + xlab('Time') +
  ylab(' ')












install.packages("timetk")
install.packages("fpp3")
install.packages("tsibble")

library(timetk)
library(fpp3)
library(tsibble)

### 3a. 
# Transform the series to a tsibble 
data("taylor_30_min", package = "timetk")
taylor_30_min_ts <- as_tsibble(taylor_30_min)

# Check the daily and weekly seasonality
gg_season(taylor_30_min_ts, value, period = "day") + 
  theme(legend.position = "none") +
  labs(x="Hour", y="Megawatts", title="Electricity demand in England and Wales")

gg_season(taylor_30_min_ts, value, period = "week") + 
  theme(legend.position = "none") +
  labs(x="Day", y="Megawatts", title="Electricity demand in England and Wales")


### 3b.
# Comment about the presence of daily and weekly seasonal effect

### 4a. 
install.packages("Ecdat")
library(Ecdat)

data("Hstarts", package = "Ecdat")
Hstarts_ts <- as_tsibble(Hstarts)
Hstarts_tsible <- filter(Hstarts_ts, key == "hs")

# Create a training set and test set
trainingset <- filter_index(Hstarts_tsible, "1960 Q1" ~ "1994 Q4")
testset <- filter_index(Hstarts_tsible, "1995 Q1" ~ "2001 Q4")

### 4b. 
# Decompose the training set  by means of STL decomposition
dcmp <- model(trainingset, stl = STL(value))
components(dcmp)
autoplot(components(dcmp))

### 4c. 
## Forecasting the seasonally adjust series
# Fit model:
Hstart_model_0 <- model(trainingset,
                        STL(value ~ trend(window = 7), robust = TRUE))
Hstart_model_1 <- select(components(Hstart_model_0), -.model)

# Forecast season_adjust via drift method:

Hstart_forecast_drift <- forecast(model(Hstart_model_1, 
                                        RW(season_adjust ~ drift())), 
                                  new_data = testset)

autoplot(Hstart_forecast_drift, Hstart_model_1) + 
  labs(y = "Logarithm", title = "Urban housing starts in Canada")
hilo(Hstart_forecast_drift)

# Forecast season_adjust via naive method:
Hstart_forecast_naive <- forecast(model(Hstart_model_1, 
                                        NAIVE(season_adjust)), 
                                  new_data = testset)
autoplot(Hstart_forecast_naive, Hstart_model_1) + 
  labs(y = "Logarithm", title = "Urban housing starts in Canada")
hilo(Hstart_forecast_naive)

# Forecats seasonal component via seasonal naive method:
Hstart_forecast_snaive <- forecast(model(Hstart_model_1, 
                                         SNAIVE(season_year)), 
                                   new_data = testset)
autoplot(Hstart_forecast_snaive, Hstart_model_1) + 
  labs(y = "Logarithm", title = "Urban housing starts in Canada")
hilo(Hstart_forecast_snaive)

###4d. 
Hstart_fitmodel <- model(trainingset, Mean = MEAN(value),
                         Naive = NAIVE(value),
                         'Seasonal naive' = SNAIVE(value),
                         Drift = RW(value ~ drift()))
Hstart_fc <- forecast(Hstart_fitmodel, new_data = testset)
autoplot(Hstart_fc, Hstart_model_1,level = NULL) + 
  labs(y = "Logarithm", title = "Urban housing starts in Canada")
accTable <- accuracy(Hstart_fc,testset)
select(accTable,.model,RMSE,MAE,MAPE)

