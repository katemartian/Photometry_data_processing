
# If you would like to use this code, please cite our Jove paper:
#  Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry to Record
#  Neural Activity in Freely Moving Animal. J. Vis. Exp. (152), e60278, 
#  doi:10.3791/60278 (2019).
#  https://www.jove.com/video/60278/multi-fiber-photometry-to-record-neural-activity-freely-moving


########### Your data ##########################################################

df <- read.csv('example.csv') # Change to your file and directory
head(df, n=5)

raw_reference <- df$MeanInt_410nm[2:nrow(df)]
raw_signal <- df$MeanInt_470nm[2:nrow(df)]

# Plot raw data
par(mfrow=c(2,1))
plot(raw_reference, type='l', col='purple')
plot(raw_signal, type='l', col='blue')


########## Use function ########################################################

# Run function get_zdFF.R
zdFF = get_zdFF(raw_reference,raw_signal)

# Plot z-score dF/F
par(mfrow=c(1,1))
plot(zdFF, type='l', col='black')

########## Analysis step by step ###############################################

### Smooth

# Run function movavg.R
smooth_win = 10
smooth_reference <- movavg(raw_reference, smooth_win)
smooth_signal <- movavg(raw_signal, smooth_win)

# Plot smoothed signal
par(mfrow=c(2,1))
plot(smooth_reference, type='l', col='purple')
plot(smooth_signal, type='l', col='blue')


### Find slope baseline

# Install package airPLS 
install.packages('devtools')
library(devtools)
httr::set_config( httr::config( ssl_verifypeer = 0L ) )
install_github("zmzhang/airPLS_R")

# Call airPLS library
library(airPLS)

base_r <- airPLS(smooth_reference,5e4,1,50)
base_s <- airPLS(smooth_signal,5e4,1,50)

# Plot data and baselines
par(mfrow=c(2,1))
plot(smooth_reference, type='l', col='purple')
lines(base_r)
plot(smooth_signal, type='l', col='blue')
lines(base_s)


### Remove baseline and the begining

remove = 200
n = length(raw_reference)
reference <- smooth_reference[remove:n] - base_r[remove:n]
signal <- smooth_signal[remove:n] - base_s[remove:n]

# Plot flatten data
par(mfrow=c(2,1))
plot(reference, type='l', col='purple')
plot(signal, type='l', col='blue')


### Standardize signals

z_reference <- (reference - median(reference)) / sd(reference)
z_signal <- (signal - median(signal)) / sd(signal)

# Plot standardized data
par(mfrow=c(2,1))
plot(z_reference, type='l', col='purple')
plot(z_signal, type='l', col='blue')


### Linear robust fit
require(MASS)
fit <- rlm(z_signal ~ z_reference)

# Plot fit
par(mfrow=c(1,1))
plot(reference, signal)
abline(fit, col="red")


### Align signals

z_reference_fit <- predict(fit)

# Plot aligned signals
plot(signal, type='l', col='blue')
lines(reference, type='l', col='purple')


### Calculate z-score dF/F

zdFF <- z_signal - z_reference_fit

# Plot z-score dF/F
plot(zdFF, type='l', col='black')


########## Contuct us ##########################################################

# If you have any questions, please contact us: ekaterina.martianova.1@ulaval.ca
