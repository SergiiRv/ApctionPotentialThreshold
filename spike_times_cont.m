function [N ,out1] = spike_times_cont(trace,threshold1)

%   This function detects and locates the time points of action potentials in a trace of 
%   membrane potential as a function of time in a neuron. The trace should represent
%   a current clamp recording from a neuron.
%   Input: 
%   "trace" is the membrane voltage array of the neuron
%   "Theshold" is the value for the spike to cross to be detected.
%   Output:
%   The output array is the index location of spikes.
%
%   Rune W. Berg 2006
%   rune@berg-lab.net
%   www.berg-lab.net
%   Modified by Rune Berg May 2015

 gim=trace;

    clear('set_crossgi')
    set_crossgi=find(gim > threshold1);  % setting the threshold
    %% %SR
    %first = min(set_crossgi)
    %last = max(set_crossgi)
    %length_set_crossgi = length (set_crossgi)
    %% %SR
    %plot (gim); hold on
    %plot (set_crossgi(1750:1800)) 
    %hold on
    %plot (set_crossgi(1:10), gim(set_crossgi(1:10)))
    %plot (set_crossgi, gim(set_crossgi))
    %%
    clear('index_shift_neggi');clear('index_shift_pos');

if isempty(set_crossgi) < 1     % This to make sure there is a spike otherwise the code below gives problems. There is an empty else statement below.

    clear('set_cross_plusgi');clear('set_cross_minus')
    %SR potentially fix here...
    index_shift_posgi(1)=min(set_crossgi);
    index_shift_neggi(length(set_crossgi))=max(set_crossgi);
    %SR
    %%
% SR the condition was changed from 1 to 5 for better sustainability...
    for i=1:length(set_crossgi)-1
     if set_crossgi(i+1) > set_crossgi(i)+5 ; 
     index_shift_posgi(i+1)=i;
     index_shift_neggi(i)=i;
     end
    end

    %Identifying up and down slopes:

    set_cross_plusgi=  set_crossgi(find(index_shift_posgi)) ;  % find(x) returns nonzero arguments.
    set_cross_minusgi=  set_crossgi(find(index_shift_neggi));   % find(x) returns nonzero arguments.
    set_cross_minusgi(length(set_cross_plusgi))= set_crossgi(end);
    nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows.

    %% Getting the spike coords
    %spikemax = zeros(1,nspikes);
    for i=1: nspikes
            spikemax(i)=min(find(gim(set_cross_plusgi(i):set_cross_minusgi(i)) == max(gim(set_cross_plusgi(i):set_cross_minusgi(i))))) +set_cross_plusgi(i)-1;
    end

else
    spikemax=[];
    display('no spikes in trace')
end 
N=length(spikemax) ;
out1=spikemax;
end