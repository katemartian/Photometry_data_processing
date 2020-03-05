'''
get_zdFF.py calculates standardized dF/F signal based on calcium-idependent 
and calcium-dependent signals comomnly recorded using fiber photometry calcium imaging

Ocober 2019 Ekaterina Martianova ekaterina.martianova.1@ulaval.ca 

Reference:
  (1) Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry 
      to Record Neural Activity in Freely Moving Animal. J. Vis. Exp. 
      (152), e60278, doi:10.3791/60278 (2019)
      https://www.jove.com/video/60278/multi-fiber-photometry-to-record-neural-activity-freely-moving

'''

def get_zdFF(raw_reference,raw_signal,smooth_win,remove,lambd,porder,itermax): 
  '''
  Calculates z-score dF/F signal based on fiber photometry calcium-idependent 
  and calcium-dependent signals
  
  Input
      raw_reference - calcium-independent signal (usually 405-420 nm excitation)
      raw_signal - calcium-dependent signal (usually 465-490 nm excitation for 
                   green fluorescent proteins, or ~560 nm for red)
      smooth_win - window for moving average smooth 
      remove - the beginning of the traces with a big slope one would like to remove
      Inputs for airPLS:
      lambd - parameter that can be adjusted by user. The larger lambda is,  
              the smoother the resulting background, z
      porder - adaptive iteratively reweighted penalized least squares for baseline fitting
      itermax: maximum iteration times
  Output
      zdFF - z-score dF/F, numpy vector
  '''
  
  import numpy as np
  from sklearn.linear_model import Lasso

 # Smooth signal
  reference = np.array(raw_reference.rolling(window=smooth_win).mean()).reshape(len(raw_reference),1)
  signal = np.array(raw_signal.rolling(window=smooth_win).mean()).reshape(len(raw_signal),1)
  
 # Remove slope using airPLS algorithm
  r_base=airPLS(raw_reference.T,lambda_=lambd,porder=porder,itermax=itermax).reshape(len(raw_reference),1)
  s_base=airPLS(raw_signal.T,lambda_=lambd,porder=porder,itermax=itermax).reshape(len(raw_signal),1)  

 # Remove baseline and the begining of recording
  reference = (reference[remove:] - r_base[remove:])
  signal = (signal[remove:] - s_base[remove:])   

 # Standardize signals    
  z_reference = (reference - np.median(reference)) / np.std(reference)
  z_signal = (signal - np.median(signal)) / np.std(signal)
  
 # Align reference signal to calcium signal using non-negative robust linear regression
  lin = Lasso(alpha=0.0001,precompute=True,max_iter=1000,
              positive=True, random_state=9999, selection='random')
  lin.fit(z_reference, z_signal)
  z_reference_fitted = lin.predict(z_reference).reshape(len(z_reference),1) 

 # z dFF    
  zdFF = (z_signal - z_reference_fitted)
 
  return zdFFf