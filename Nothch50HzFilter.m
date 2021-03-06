function y = Nothch50HzFilter(x)
%50HZ FILTER Filters input x and returns output y.

% MATLAB Code
% Generated by MATLAB(R) 8.5 and the DSP System Toolbox 9.0.
% Generated on: 16-May-2016 22:52:07

%#codegen

% To generate C/C++ code from this function use the codegen command.
% Type 'help codegen' for more information.

persistent Hd;

if isempty(Hd)
    
    % The following code was used to design the filter coefficients:
    %
    % N      = 10;      % Order
    % Fpass1 = 45;      % First Passband Frequency
    % Fpass2 = 55;      % Second Passband Frequency
    % Apass  = 1;       % Passband Ripple (dB)
    % Fs     = 100000;  % Sampling Frequency
    %
    % h = fdesign.bandstop('n,fp1,fp2,ap', N, Fpass1, Fpass2, Apass, Fs);
    %
    % Hd = design(h, 'cheby1', ...
    %     'SystemObject', true);
    
    Hd = dsp.BiquadFilter( ...
        'Structure', 'Direct form II', ...
        'SOSMatrix', [1 -1.99999022909928 1 1 -1.99992548531524 ...
        0.999937431757629; 1 -1.99999022909928 1 1 -1.99994083466247 ...
        0.999948825691103; 1 -1.99999022909928 1 1 -1.99959561906516 ...
        0.999608622419489; 1 -1.99999022909928 1 1 -1.99969859532012 ...
        0.999705934613004; 1 -1.99999022909928 1 1 -1.99782218486054 ...
        0.997831945169336], ...
        'ScaleValues', [0.999971514828317; 0.999971514828317; ...
        0.999828538420112; 0.999828538420112; 0.998915972584668; 1]);
end

y = step(Hd,x);