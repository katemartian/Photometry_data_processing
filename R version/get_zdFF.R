get_zdFF <- function(reference, signal, smooth_win=10, remove=200, lambda=5e4, itermax=50, differences=1)
{
  # Calculates z-score dF/F signal based on fiber photometry calcium-idependent 
  # and -dependent signals.
  #
  # This program is a translation in R of the Python source code of get_zdFF.py
  #
  # Input 
  #     reference: calcium-independent signal (usually 405-420 nm excitation)
  #     signal: calcium-dependent signal (usually 465-490 nm excitation 
  #             for green fluorescent proteins, or ~560 nm for red)
  #     smooth_win: window for moving average smooth 
  #     remove: the beginning of the traces with a steep slope one would like 
  #             to remove
  #  Inputs for airPLS:
  #   lambda: lambda is an adjustable parameter, it can be adjusted by user. 
  #             The larger lambda is, the smoother baseline will be 
  #     itermax: maximum iteration times
  #     differences
  #
  #  Output
  #     zdFF - z-score dF/F, 
  #    
  #  Examples:
  #     zdFF = get_zdFF(reference, signal);
  #     zdFF = get_zdFF(reference, signal, 10, 200, 5e4, 50, 1);
  #
  #  Reference:
  #    (1) Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry 
  #         to Record Neural Activity in Freely Moving Animal. J. Vis. Exp. 
  #         (152), e60278, doi:10.3791/60278 (2019)
  #         https://www.jove.com/video/60278/multi-fiber-photometry-to-record-neural-activity-freely-moving
  #
  #  March 2020 Ekaterina Martianova ekaterina.martianova.1@ulaval.ca
    
  
 # Smooth signals
  reference = movavg(reference, smooth_win)
  signal = movavg(signal, smooth_win)
  
 # Find slope using airPLS algorithm
  require(airPLS)
  base_r <- airPLS(reference, lambda, differences, itermax)
  base_s <- airPLS(signal, lambda, differences, itermax)
  
 # Remove slope and the begining of the recordings
  n = length(reference)
  reference = reference[remove:n] - base_r[remove:n]
  signal = signal[remove:n] - base_s[remove:n]
  
 # Standardize signals
  reference = (reference - median(reference)) / sd(reference)
  signal = (signal - median(signal)) / sd(signal)
  
 # Linear robust fit
  require(MASS)
  fit <- rlm(signal ~ reference)
  reference <- predict(fit)
  
 # Calculate z-score dF/F
  zdFF <- signal - reference
  
  return(zdFF)
  
}