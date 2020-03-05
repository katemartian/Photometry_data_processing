function zdFF = get_zdFF(reference, signal, ...
                smooth_win, remove, lambda, itermax, order, wep, p)

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
                        lambda=10e9;                        
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