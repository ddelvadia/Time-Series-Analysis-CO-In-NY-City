# Final Project
# Name: Dhaval Delvadia
# Date: 2019-03-06
# Version: 5
# Write what this dataset is about

# install packages
install.packages("tidyverse")

# Load libraries
library(tidyverse)
library(tseries)
library(lubridate)
library(fBasics)
library(fUnitRoots)
library(TSA)
library(forecast)
library(lmtest)
library(dplyr)


# set working directory
setwd("~/CSC 425/Final/pollution_us_2000_2016.csv")

# read complete dataset of the US pollution
df = read.csv("pollution_us_2000_2016.csv")

# show the first five points of the dataset dataframe
head(df)

# let's look at the structure of the dataset
str(df, vec.len=1)

# based on the structure, it looks like column X is just the index column
# Address, state, country, city, Date.Local, No2.Units, O3.Units, SO2.Units, Co.Units are all strings
# rest of the variables are either int or num

# let's convert the date field from string to actual date and then convert that date to time series
df$Date.Local = as.Date(df$Date.Local, format="%Y-%m-%d")
head(df$Date.Local)


# we are only interested in city of New York dataset to be specific.
# So we will make new dataframe for New York
ny=dplyr::filter(df, City=="New York")

#summary of the new york dataset Only
summary(ny)

# STEP 1
# Since I chose to do the Ozone Series from the NY dataset, we then focus on the mean values of each
# pollutants only.
# aggregate mean data of O3 based on week
ny_week=aggregate(ny$O3.Mean ~ week(ny$Date.Local) + year(ny$Date.Local), data = ny, FUN=mean)
#class(ny_week)


# STEP 2:
# Plot of the weekly time series
tsWeek=ts(ny_week$`ny$O3.Mean`, start=c(2000,1), frequency = 53)
plot.ts(tsWeek, main="Weekly NY Mean O3 value from 2000 to 2016", xlab='Time', ylab='O3 mean')


# STEP 3:
# let's see decompose graph of weekly and monthly O3.Mean for the year 2000 to 2016
# We see clearly there is a kind of exponential trend upward
# Seasonality distribution for frequency 53 was best. It was 0.02 half as much as data distribution
# 0.04
# Based on the weekly decompose graph, we can see exponential upward trend and may be seasonality
plot(decompose(tsWeek, type="additive")) 
plot(stl(tsWeek, s.window=53)) #decompose with moving average window of 53 for # of weeks in a year


# STEP 4:
# let's find out the distribution of weekly and monthly O3 mean values in NY
hist(tsWeek, prob=TRUE, main = "Histogram of the Weekly NY O3 Mean", xlab='Aggragated O3 Mean Per Week')
xfit<-seq(min(tsWeek), max(tsWeek))
yfit<-dnorm(xfit, mean=mean(tsWeek), sd=sd(tsWeek))
lines(xfit,yfit, col='blue', lwd=2)

# STEP 5:
# Normal plots for weekly mean Ozone for NY
plot(qqnorm(tsWeek), main = 'normal probability plot of weekly mean Ozone for NY',
     xlab='Mean Ozone',
     ylab='Mean Ozone')
qqline(tsWeek, col = 'red') 


# STEP 6:
# NORMALITY TESTS FOR WEEKLY MEAN of pollutant O3
# Perform Jarque-Bera normality test.
normalTest(tsWeek,method=c('jb')) 


## ny weekly return
## nyWkRtn = diff(ny_week$`ny$O3.Mean`)/ny_week$`ny$O3.Mean`[-nrow(ny_week)]
## plot(nyWkRtn, type="l")
## qplot(nyWkRtn)


# STEP 7
# ACF Weekly: 
# Analyze if the time series is serially correlated using the ACF plot and the Ljung Box test.
#acf(df$index, lag.max = 15, plot = F) # we do not plot the output
acf(tsWeek, lag.max = 50, main="ACF plot of weekly Mean O3 for NY")
#acf(ny_week$`ny$O3.Mean`, lag.max=50, main ='ACF plot of weekly Mean O3 for NY')


# STEP 8
# Let's analyze the PACF graphs
pacf(tsWeek, lag.max = 20, main='PACF plot of Weekly Mean O3 for NY')
#pacf(ny_week$`nyO3mean$O3.Mean`, lag.max = 50, main='PACF plot of Weekly O3 mean for NY')

# STEP 9
# compute the Ljung-Box Test for weekly mean O3 for lag 5 and 25. Both shows p-values less than
# 2.2e-16. This means we can reject null hypothesis that lags are independent and atleast lags are
# not independent.
Box.test(tsWeek, lag=5, type='Ljung') 
Box.test(tsWeek, lag=25, type='Ljung') 
#Box.test(ny_week$`ny$O3.Mean`, lag=5, type='Ljung') 
#Box.test(ny_week$`ny$O3.Mean`, lag=25, type='Ljung')


# STEP 10
# APPLY EACF MODEL SELECTION PROCEDURES
# use TSA package
# install.packages("TSA")
# library(TSA)
print('Weekly AR/MA')
#compute weekly eacf values for p<=10 and q<=10
eacf(tsWeek, 10,10)


# Based on the eacf of the ny_week output, we can say that the model can be made from
# ARMA with AR = 1 and MA = 2

#print(m1$eacf, digits=2)


# STEP 11
# FIT AR(P) MODEL (ORDER AUTOMATICALLY SELECTED)
# Based on the output of the ACF and PACF plot, we can see that ACF does not go
# to zero quickly while PACF does. Thus, it can be a AR model. We will try to use ar()
# to find the AR model.
m1 = ar(tsWeek, method="mle")
m1
# Based on the output looks like, the order was selected as 12 

#to print internal variables computed for AR model m1
names(m1)
# print aic values
m1$aic


#Another approach to fit AR model of given order
m2=arima(tsWeek, c(1,0,2), method="ML")  
coefficients(m2)
coeftest(m2)


# STEP 11:
# let's find the models for weekly using the ARIMA
# ny_weekly = mean O3 values per given week in NY
# ln_wk = log of mean O3 values in a given week in NY
plot(tsWeek, type='l') # plot just the weekly time series of mean O3 of NY 
# Now, create the log of mean O3 for each week
ln_wk = log(tsWeek)
plot(ln_wk, type='l') # still shows variation and trend slowly going up; let's check using decompose
plot(decompose(ln_wk)) # shows time trend going up so we will need to apply return

eacf(ln_wk) # still see the AR(0) row having insignificant values so we will have to take the difference


# And finally the first difference = the log returns
rtn_wk = diff(ln_wk)
plot(rtn_wk, type='l') # now looks stationary; however, has variations; will have to apply GARCH later
eacf(rtn_wk)


# STEP 12 
# COMPUTE SUMMARY STATISTICS
basicStats(tsWeek)
basicStats(ln_wk)
basicStats(rtn_wk)


# STEP 13
################### HW#3 Problem 2: ########################
# The weekly time series looks like a linear trend along with seasonality. 
# There may just be one return series. let's check with difference 

# 1. Run a Dickey-Fuller test, there are three versions of this test, 
# so each member should try a different version of this test. If you 
# have more than 3 participants, then you should try get related series 
# that you can test as well.

# Dickey Fuller Test on the raw weekly Ozone series for New York with various lag values
# However, we since the data shows seasonality and this is weekly series with data going
# more than one year cycle.  We will have to start atleast at lag 53 to see the correct
# output from the Dickey Fuller Test. 
adfTest(coredata(tsWeek), lags=53, type=c("nc")) # DFT with "nc" = no intercept nor time trend (53 weeks lags)
adfTest(coredata(tsWeek), lags=70, type=c("nc")) # DFT with "nc" = no intercept nor time trend (70 weeks lags)
adfTest(coredata(tsWeek), lags=90, type=c("nc")) # DFT with "nc" = no intercept nor time trend (90 weeks lags)

adfTest(coredata(tsWeek), lags=53, type=c("c")) # DFT with "c" = constant intercept but no time trend (53 weeks lags)
adfTest(coredata(tsWeek), lags=70, type=c("c")) # DFT with "c" = constant intercept but no time trend (70 weeks lags)
adfTest(coredata(tsWeek), lags=90, type=c("c")) # DFT with "c" = constant intercept but no time trend (90 weeks lags)

adfTest(coredata(tsWeek), lags=53, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)
adfTest(coredata(tsWeek), lags=70, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)
adfTest(coredata(tsWeek), lags=90, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)


plot(rtn_wk)

# Dickey Fuller Test on the returned weekly Ozone series for New York with various log values
adfTest(coredata(rtn_wk), lags=3, type=c("nc")) 
adfTest(coredata(rtn_wk), lags=6, type=c("nc")) 

adfTest(coredata(rtn_wk), lags=3, type=c("c")) 
adfTest(coredata(rtn_wk), lags=6, type=c("c"))

adfTest(coredata(rtn_wk), lags=3, type=c("ct")) 
adfTest(coredata(rtn_wk), lags=6, type=c("ct"))


# STEP 14
# let's check the order of this 
# Based on the outcome of the raw data for eacf on tsWeek, we can see that it is has seasonality
# as the colum 6 has signifacnt values
eacf(tsWeek)

# Based on the output of the eacf for the difference return
# we can see that we can pick AR=0 and MA=1
# we can check with auto.arima to see if we can same value
eacf(ln_wk)
eacf(rtn_wk)
acf(ln_wk)
pacf(ln_wk)

acf(rtn_wk)
pacf(rtn_wk)

# STEP 15
# FIT AN INITIAL ARIMA(1,0,1) MODEL based on the eacf(rtn_wk) return
#
# Fit an ARIMA model to a univariate time series.
# arima(x, order = c(p, 1, q)), since one level of differencing
# resulted in a stationary series
#lnWeek=log(tsWeek)
m3 = arima(ln_wk, order=c(0,1,1), include.mean=T)
m3
coeftest(m3)   


# STEP 16
### Based on the initial results from problem 1, run an initial model for one of the series 
# in your set. You can choose any order ARIMA (or GARCH or other models if you want to read 
# ahead ???) that you want, but each member should try a different model. 
#m3=auto.arima(tsWeek, trace=TRUE, ic="bic", stationary = FALSE, allowdrift = F)
m4=auto.arima(ln_wk, trace=TRUE, ic="bic", stationary = FALSE, allowdrift = TRUE) #Best model ARIMA(1,0,0)(0,1,0)
coeftest(m4)


# All models m2, m3, and m4 are different from one another.

# this graph seems to be more like a white noise which is what it should look like after taking the 
# first difference of the linear model;

# Now we should check to see if model m3 residual contain white noise using Ljung Box Test
Box.test(m3$residuals, lag=53, type='Ljung')
Box.test(m3$residuals, lag=60, type='Ljung')

# We will perform Jarque-Bera normality test on m1 residual.
normalTest(m3$residuals, method=c('jb'))  


# Now we should check to see if model m4 residual contain white noise using Ljung Box Test
Box.test(m4$residuals, lag=53, type='Ljung')
Box.test(m4$residuals, lag=60, type='Ljung')

# We will perform Jarque-Bera normality test on m1 residual.
normalTest(m4$residuals, method=c('jb'))  


plot(m2$residuals)
plot(m3$residuals)
plot(m4$residuals)

# we will now perform backtesting with model m1 (weekly) to see which model is better
source("backtest.R")
length(tsWeek)
test = round(length(tsWeek) * .9)
pm2 = backtest(m2, orig=test, rt=tsWeek, 1)
pm3 = backtest(m3, orig=test, rt=ln_wk, 1)
pm4 = backtest(m4, orig=test, rt=ln_wk, 1)
# Based on the results we see that Root Mean Square Error of out-of-sample forecast is 612. This means
# square root of squared errors of out of samples 17 next samples are 612. This seems a lot but the
# mean absolute percentage error is 3.4% which is not that much.  So we can say our model m3 is 96.6% 
# accurate.


# Now let's check weekly model
# Now we should check to see if model m2 residual contain white noise using Ljung Box Test
Box.test(m2$residuals, lag=6, type='Ljung')
Box.test(m2$residuals, lag=53, type='Ljung')

Box.test(m4$residuals, lag=6, type='Ljung')
Box.test(m4$residuals, lag=53, type='Ljung')

# We will perform Jarque-Bera normality test on m1 residual.
normalTest(m2$residuals, method=c('jb'))
normalTest(m4$residuals, method=c('jb'))




# seasonality on weekly series
# I do not see any seasonality from the decompose plot. However, the distribution of the seasonal
# is half the amound for the actual model. Eacf is also does not show any colums completely insignificatnt


# GARCH effects on weekly series







################################ NOW MONTHLY O3 MEAN ##############################################
###################################################################################################
# STEP 1: 
# aggregate mean data of O3 based on month
ny_month=aggregate(ny$O3.Mean ~ month(ny$Date.Local) + year(ny$Date.Local), data = ny, FUN=mean)

# STEP 2:
# monthly time series
tsMonth=ts(ny_month$`ny$O3.Mean`, start=c(2000,1), frequency = 12)
plot.ts(tsMonth, main="Monthly NY Mean O3 value from 2000 to 2016", xlab='Time', ylab='O3 mean')

# STEP 3:
# let's see decompose graph of weekly and monthly O3.Mean for the year 2000 to 2016
plot(decompose(tsMonth)) #Based on monthly decompose graph, we can see upward trend and seasonality

# STEP 4:
# let's find out the distribution of weekly and monthly O3 mean values in NY
hist(tsMonth, main = "Histogram of the Monthly NY O3 Mean", xlab='Aggragated O3 Mean Per Month') # this is skewed to the right

# STEP 5:
# Normal plots for monthly mean Ozone for NY
plot(qqnorm(tsMonth), main = 'normal probability plot of monthly mean Ozone for NY',
     xlab='Mean Ozone',
     ylab='Mean Ozone')
qqline(tsMonth, col = 'red') 


# STEP 6:
# NORMALITY TESTS FOR MONTHLY MEAN of pollutant O3
# Perform Jarque-Bera normality test.
normalTest(tsMonth,method=c('jb')) 


# STEP 7:
# Monthly:
# Analyze if the time series is serially correlated using the ACF plot and the Ljung Box test.
#acf(df$index, lag.max = 15, plot = F) # we do not plot the output
acf(tsMonth, lag.max = 50, main="ACF plot of Monthly Mean O3 for NY")
#acf(ny_month$`ny$O3.Mean`, lag.max=50, main ='ACF plot of Mean O3 for NY')


# STEP 8:
# Let's analyze the PACF graphs
pacf(tsMonth, lag.max = 50, main='PACF plot of Monthly Mean O3 for NY')
#pacf(ny_month$`ny$O3.Mean`, lag.max = 50, main='PACF plot of Monthly Mean O3 for NY')


# STEP 9:
# compute the Ljung-Box Test for monthly mean O3. Both shows p-values less than
# 2.2e-16. This means we can reject null hypothesis that lags are independent and atleast lags are
# not independent.
Box.test(tsMonth, lag=5, type='Ljung')
Box.test(tsMonth, lag=20, type='Ljung')


# STEP 10:
# APPLY EACF MODEL SELECTION PROCEDURES
# use TSA package
# install.packages("TSA")
# library(TSA)
print('Monthly AR/MA')
#compute monthly eacf values
mm1=eacf(tsMonth, 10,10)
names(mm1)
# Based on the eacf of the ny_month ouput, we can say that the ARMA model will have
# p=2 and q=1

#print(mm1$eacf, digits=2)


# STEP 11:
# let's find the models for monthly using the ARIMA
# ny_month = mean O3 values in a given month in NY 
# ln_mth = log of mean O3 values in a given month in NY

plot(tsMonth, type='l') # plot just the weekly time series of mean O3 of NY 

# Now, create the log of mean O3 for each week
ln_mth = log(tsMonth)
plot(ln_mth, type='l') # still shows variation and trend slowly going up; let's check using decompose
plot(decompose(ln_mth)) # shows time trend going up so we will need to apply return
# And finally the first difference = the log returns
rtn_mth = diff(ln_mth, 12)
plot(rtn_mth, type='l') # now looks stationary; however, has variations; will have to apply GARCH later


# STEP 12 
# COMPUTE SUMMARY STATISTICS
basicStats(tsMonth)
basicStats(ln_mth)
basicStats(rtn_mth)


# STEP 13
################### HW#3 Problem 2: ########################
# The weekly time series looks like a linear trend along with seasonality. 
# There may just be one return series. let's check with difference 

# 1. Run a Dickey-Fuller test, there are three versions of this test, 
# so each member should try a different version of this test. If you 
# have more than 3 participants, then you should try get related series 
# that you can test as well.

# Dickey Fuller Test on the raw weekly Ozone series for New York with various lag values
# However, we since the data shows seasonality and this is weekly series with data going
# more than one year cycle.  We will have to start atleast at lag 53 to see the correct
# Dickey Fuller Test on the raw monthly Ozone series for New York with various log values
adfTest(coredata(tsMonth), lags=12, type=c("nc")) # DFT with "nc" = no intercept nor time trend (53 weeks lags)
adfTest(coredata(tsMonth), lags=24, type=c("nc")) # DFT with "nc" = no intercept nor time trend (70 weeks lags)
adfTest(coredata(tsMonth), lags=36, type=c("nc")) # DFT with "nc" = no intercept nor time trend (90 weeks lags)

adfTest(coredata(tsMonth), lags=12, type=c("c")) # DFT with "c" = constant intercept but no time trend (53 weeks lags)
adfTest(coredata(tsMonth), lags=24, type=c("c")) # DFT with "c" = constant intercept but no time trend (70 weeks lags)
adfTest(coredata(tsMonth), lags=36, type=c("c")) # DFT with "c" = constant intercept but no time trend (90 weeks lags)

adfTest(coredata(tsMonth), lags=12, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)
adfTest(coredata(tsMonth), lags=24, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)
adfTest(coredata(tsMonth), lags=36, type=c("ct")) # DFT with "ct"=constant intercept and a time trend (53 weeks lags)

# let's try return/difference for weekly time series with weekly adjusted data by subtracting 53 weeks
#diff_tsMonth=diff(log(tsMonth), 12) # adjusting by subtracting 53 weeks cycle
#plot(diff_tsMonth)

# Dickey Fuller Test on the returned weekly Ozone series for New York with various log values
adfTest(coredata(rtn_mth), lags=3, type=c("nc")) 
adfTest(coredata(rtn_mth), lags=6, type=c("nc")) 

adfTest(coredata(rtn_mth), lags=3, type=c("c")) 
adfTest(coredata(rtn_mth), lags=6, type=c("c"))

adfTest(coredata(rtn_mth), lags=3, type=c("ct")) 
adfTest(coredata(rtn_mth), lags=6, type=c("ct")) 

# STEP 14
# let's check the order of this 
# Based on the outcome of the raw data for eacf on tsMonth, we can see that it is has seasonality
# as the colums 0 and 11 are signifacnt values
eacf(tsMonth)

# Based on the output of the eacf for the difference return
# we can see that we can pick AR=0 and MA=1
# we can check with auto.arima to see if we can same value
eacf(rtn_mth)



# STEP 15
# FIT AN INITIAL ARIMA(0,1,1) MODEL based on the eacf(rtn_wk) return in step 15
#
# Fit an ARIMA model to a univariate time series.
# arima(x, order = c(p, 1, q)), since one level of differencing
# resulted in a stationary series
mm2 = arima(ln_mth, order=c(0,1,1), include.mean=T)
mm2
coeftest(mm2) 


# STEP 16
### Based on the initial results from problem 1, run an initial model for one of the series 
# in your set. You can choose any order ARIMA (or GARCH or other models if you want to read 
# ahead ???) that you want, but each member should try a different model. 
#m3=auto.arima(tsWeek, trace=TRUE, ic="bic", stationary = FALSE, allowdrift = F)
mm3=auto.arima(ln_mth, trace=TRUE, ic="bic", stationary = FALSE, allowdrift = TRUE) #Best model ARIMA(1,0,0)(0,1,0)

# Both m1 and m2 are different not sure why

# this graph seems to be more like a white noise which is what it should look like after taking the 
# first difference of the linear model;

# Now we should check to see if model m1 residual contain white noise using Ljung Box Test
Box.test(mm3$residuals, lag=6, type='Ljung')
Box.test(mm3$residuals, lag=12, type='Ljung')

# We will perform Jarque-Bera normality test on m1 residual.
normalTest(mm3$residuals, method=c('jb'))  

# we will now perform backtesting with model m1 (weekly) to see which model is better
source("backtest.R")
length(tsWeek)
window = length(tsWeek[780:866]) * .9
pm1 = backtest(m3, tsWeek[780:866], window, 1)
# Based on the results we see that Root Mean Square Error of out-of-sample forecast is 612. This means
# square root of squared errors of out of samples 17 next samples are 612. This seems a lot but the
# mean absolute percentage error is 3.4% which is not that much.  So we can say our model m3 is 96.6% 
# accurate.

# Now let's check monthly model
# Now we should check to see if model m2 residual contain white noise using Ljung Box Test
Box.test(m2$residuals, lag=6, type='Ljung')
Box.test(m2$residuals, lag=12, type='Ljung')

# We will perform Jarque-Bera normality test on m1 residual.
normalTest(m2$residuals, method=c('jb'))  

# we will now perform backtesting with model m1 (weekly) to see which model is better
length(tsMonth)
window = length(tsMonth[780:866]) * .9
pm2 = backtest(m2, tsMonth[780:866], window, 1)