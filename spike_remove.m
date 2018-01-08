function [N,out1] = spike_remove(trace1, prewin, postwin, threshold1)   
%   This function replaces all the spikes in the input trace
%   "trace1" in a window "prewin" before the spike and "postwin" after
%   the spike, with the values right before and right after, connected by a straight line. 
% 
%%  input:
%
%   trace1:  The intracellular recording with action potentials to be
%   replaced.
%   prewin, postwin: The windows before and after the peak of spike to
%   exclude in data points (not in time).
%   threshold1: The membrane potential value for the action potentials to
%   cross in order to be identified. This is typically 0, but if the spikes 
%   are smaller than zero then this has to be adjusted.
%
%%  Output:   
%   N:  Number of spikes detected and replaced
%   out1: The output trace with the spikes replaced by straight lines.
%
%   Sample:
%   [N,out1] = spike_replace(trace1, 40, 40, 0)
%
%   Rune W. Berg 2015
%   University of Copenhagen
%   rune@berg-lab.net, runeb@sund.ku.dk
%   www.berg-lab.net
%   Write me for comments of questions!
%   
traceout=trace1;

clear('set_crossgi')
set_crossgi=find(trace1(1:end) > threshold1)  ;  % setting the threshold

if isempty(set_crossgi) < 1     % This to make sure there is a spike otherwise the code below gives problems. There is an empty else statement below.
    
    clear('index_shift_neggi');clear('index_shift_pos');
    clear('set_cross_plusgi');clear('set_cross_minus')

    index_shift_posgi(1)=min(set_crossgi);
    index_shift_neggi(length(set_crossgi))=max(set_crossgi);

   for i=1:length(set_crossgi)-1
     if set_crossgi(i+1) > set_crossgi(i)+1 ; 
         index_shift_posgi(i+1)=i;
        index_shift_neggi(i)=i;
     end
   end

%These are the coords in the nerve based smooothing:
    set_cross_plusgi=  set_crossgi(find(index_shift_posgi));   % find(x) returns nonzero arguments.
    set_cross_minusgi=  set_crossgi(find(index_shift_neggi));   % find(x) returns nonzero arguments.
    set_cross_minusgi(length(set_cross_plusgi))= set_crossgi(end);

    nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows.

    clear('y');clear('y2')
    base1gi=0:1:length(trace1)-1;    y(1:nspikes)=0;
    base2gi=0:1:length(set_crossgi)-1;     y2(1:length(set_crossgi))=0;
    base3gi=0:1:length(trace1)-1;    y(1:nspikes)=0;
    
%% Making sure that the beginning and end are inside the whole window
% avoiding the endpoint problems:
        subtrac=nspikes-min(find(set_cross_minusgi > length(trace1)-postwin ))+1 ;
        if isempty(subtrac) > 0 
            subtrac1=0;
        else
         subtrac1=subtrac
        end
    
    for i=1: nspikes - subtrac1
        spikemax(i)=min(find(trace1(set_cross_plusgi(i):set_cross_minusgi(i)) == max(trace1(set_cross_plusgi(i):set_cross_minusgi(i))))) +set_cross_plusgi(i)-1;
   
         newbase=0:(postwin+prewin);        
         if spikemax(i) - prewin > 1
            traceout(spikemax(i) - prewin:spikemax(i) + postwin)=[[trace1(spikemax(i) + postwin+1)-trace1(spikemax(i) - prewin)]/max(newbase)]*newbase +trace1(spikemax(i) - prewin);    
         else
             traceout=trace1;
         end
    end

else

end

if isempty(set_crossgi) >0
    nspikes = 0;
end

N=nspikes; out1 = traceout ;
end