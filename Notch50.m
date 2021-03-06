function y = Notch50(x)
%NOTCH50 Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 8.5 and the DSP System Toolbox 9.0.
% Generated on: 18-May-2016 14:45:55
clear Hd;
persistent Hd;

if isempty(Hd)
    
    % IIR Notching filter designed using the IIRNOTCH function.
    
    % All frequency values are in Hz.
    Fs = 100000;  % Sampling Frequency
    
    Fnotch = 50;  % Notch Frequency
    BW     = 1;   % Bandwidth
    Apass  = 1;   % Bandwidth Attenuation
    
    [b, a] = iirnotch(Fnotch/(Fs/2), BW/(Fs/2), Apass);Hd = dsp.IIRFilter( ...
        'Structure', 'Direct form II', ...
        'Numerator', b, ...
        'Denominator', a);
end

y = step(Hd,x);


% [EOF]
