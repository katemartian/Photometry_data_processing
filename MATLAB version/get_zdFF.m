function zdFF = get_zdFF(reference, signal, smooth_win, remove, lambda, itermax, order, wep, p)
            
% Calculates z-score dF/F signal based on fiber photometry calcium-idependent 
% and calcium dependent signals.
%
% This program is a translation in MATLAB of the Python source code of
% get_zdFF.py
%
% Input 
%     reference: calcium-independent signal (usually 405-420 nm excitation)
%     signal: calcium-dependent signal (usually 465-490 nm excitation 
%              for green fluorescent proteins, or ~560 nm for red)
%     smooth_win: window for moving average smooth 
%     remove: the beginning of the traces with a steep slope one would like to remove
%   Inputs for airPLS:
%     lambda: lambda is an adjustable parameter, it can be adjusted by user. 
%             The larger lambda is, the smoother baseline will be 
%     itermax: maximum iteration times
%     order: an integer indicating the order of the difference of penalties
%     wep: weight exception proportion at both the start and end
%     p: asymmetry parameter for the start and end
%
%  Output
%     zdFF - z-score dF/F, vector
%    
%  Examples:
%     zdFF = get_zdFF(reference, signal);
%     zdFF = get_zdFF(reference, signal, 10, 200, 5e9, 50, 2, 0.5, 0.5);
%
%  Reference:
%     (1) Martianova, E., Aronson, S., Proulx, C.D. Multi-Fiber Photometry 
%         to Record Neural Activity in Freely Moving Animal. J. Vis. Exp. 
%         (152), e60278, doi:10.3791/60278 (2019)
%         https://www.jove.com/video/60278/multi-fiber-photometry-to-record-neural-activity-freely-moving
%
%  March 2020 Ekaterina Martianova ekaterina.martianova.1@ulaval.ca


  % Preset some parameters
    if nargin<9
        p = 0.5;        
        if nargin<8
            wep = 0.5;            
            if nargin<7
                order = 2;                
                if nargin<6
                    itermax = 50;                    
                    if nargin<5
                        lambda=5e9;                        
                        if nargin<4
                            remove=200;                            
                            if nargin<3
                                smooth_win=10;                                
                            end
                        end
                    end
                end
            end
        end
    end
    
  % Smooth signals using moving average
    reference = movmean(reference,smooth_win);
    signal = movmean(signal,smooth_win);
    
  % Find and remove slope in signals using airPLS
    [reference, ~]= airPLS(reference,lambda,order,wep,p,itermax);
    [signal, ~]= airPLS(signal,lambda,order,wep,p,itermax);
    
  % Remove the begining of recordings with fast drop
    reference = reference(remove:end);
    signal = signal(remove:end);
    
  % Standardize signals
    reference = (reference - median(reference)) / std(reference);
    signal = (signal - median(signal)) / std(signal);
    
  % Robust non negative linear regression
    fitdata = fit(reference',signal',fittype('poly1'),'Robust','on');
  % Align reference trace to signal using the fit
    reference = fitdata(reference)';
    
  % z-score dFF
    zdFF = signal - reference;   
end