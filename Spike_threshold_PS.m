
function [mean_thresh sd_thresh mean_thresh_hh sd_thresh_hh N thresvalue thres_coords_orginal thresvalue_hh thres_coords_orginal_hh AP_Data_Table AP_power_list]=Spike_threshold_PS(traceD,dt,win, win2, Fc,factordvdt)

%%  Spike threshold detection alogrithm - Method I
%   This function calculates the spike threshold of a membrane potential
%   recording based on the method by Sekerli et al IEEE Trans Biom Eng
%   51(9): 1665-1672, 2004. In this paper they suggest two methods of
%   estimating threshold. Here we use the first method, which is based on
%   revealing the threshold in a phase plot of V versus the derivative
%   dVdt. The largest positive slope is an indication of the threshold. 

%%  INPUT:
%   trace: Vm recording containing the spikes for which the threshold to
%   determine
%   dt: sampling interval in seconds
%   Fc: cut-off frequency for low pass filtering, default value is 5000 Hz
%   win: the half-window size in data points around the spike to show, default is
%   200 points.
%   win2: the window from the peak backwards in time in which the threshold
%   is found. Default value is 100 datapoints
%   factordvdt: is the factor times the max derivative in Vm, as a minimum to
%   include in the analysis. This is to limit the divergent data points.
%   Default is 0.03
%
%%  OUTPUT
%   mean_thresh: mean of threshold values
%   sd_thresh: standard deviation of threshold values
%   N: number of spikes
%   thresvalue: the threshold values
%   thres_coords_orginal: the x-cords of the filtered trace, which is
%   approximately the same as the original.

%   Sample :  [mean_thresh sd_thresh N thresvalue thres_coords_orginal]=Spike_threshold_PS(trace,5*10^-4,100, 100, 5000,0.03)

%   Rune W. Berg 2015
%   University of Copenhagen
%   rune@berg-lab.net, runeb@sund.ku.dk
%   www.berg-lab.net
%  
% the code is modified by S. Romanenko in 2016-17 (see details on sources in Comment.txt)
% The University of Western Australia; CRAWLEY WA 6009;
% +61 8 6488 7014; sergii.romanenko@uwa.edu.au; 

%% 

%factordvdt=0.03 ;
%win=200; win2=100;S

% Filtering the trace
%Fc =[5000]; 
%%
% SR removing the baseline:
clear('traceT')
clear('trace2')
clear('temp_X')
temp_X = traceD;
    % figure(200);
    % subplot(3, 1, 1);
    % plot(temp_X);
[baseline, baselinebig, noSpike] = baselinesmooth(traceD, dt, win, 1.2*win); %win2 after dt
    % subplot(3, 1, 2)
    % plot(baseline)
%traceD = temp_X - baseline;
traceD = smooth(temp_X - baseline,100); %traceD = smooth(temp_X - noSpike + baseline - baselinebig,100);
    % subplot(3, 1, 3)
    % plot(traceD)

%filtering
[b1,a1]=butter(2,2*Fc*dt, 'low'); % low pass butterworth filter
traceT = filter(b1,a1,traceD);  % Filtered trace to minimize noise
trace2 = Notch50(traceT);  % filtering of the 50 Hz interferience
trace2(1:100)= mean(traceD(1:100)); %rid of the initial point filtering distortions 
%%
%SR % deternmination of the spike detection operational threshold
dmax = max(trace2(:));
dmin = min(trace2(:));
height = abs(0.5*(dmin - dmax)) + dmin;

%%SR
%% 
N=0; i=0;
[N ,spikecords]=spike_times(trace2,height);   % N is the number of spikes in the trace. 
base=1:2*win+1 ; g=zeros(N,2*win+1); hh=zeros(N,2*win+1);
spiketrace = zeros(N,2*win+1); dvdt = zeros(N,2*win+1);

Am_list = zeros(2, N); Am50_list = zeros(2, N); RiseMax_list = zeros(2, N); DecayMax_list = zeros(2, N); 
ChargeIn_list = zeros(2, N); ChargeOut_list = zeros(2, N); MaxInCurr_list = zeros(2, N); MaxOutCurr_list = zeros(2, N); AP_power_list = zeros(1, N);
thres_coord = zeros(1, N); thres_coord_hh = zeros(1, N); thresvalue = zeros(1, N); thresvalue_hh = zeros(1, N);

for i=1:N
    spiketrace(i,:)=trace2(spikecords(i)-win:spikecords(i)+win);        
    
    [Am, Am_Coord] = max(spiketrace(i,:));  
    
    %spiketrace_h = spiketrace(i,1:win);
    dvdt(i,:)=[diff(spiketrace(i,:),1) 0];               
    dvdtT(i,:)= dvdt(i,1:win) ; %SR1
    %dvdtIntermed = low500pass(dvdtT(i,:)'); dvdt(i,:) = dvdtIntermed'; 
    setNZ = find((dvdtT(i,:)>(factordvdt*max(dvdtT(i,:))))&(dvdtT(i,:)>0)) ;
    
          %figure(201)
    %figure(201,'Name','dV/dt');
           %plot (dvdt(i,:),'b'); hold on; plot (dvdtT(i,:),'r'); hold off;
  
    %d2vdt2(i,:)=[diff(spiketrace(i,:),2) 0 0];          
    d2vdt2(i,:)= [diff(dvdt(i,:),1) 0];                       
    d3vdt3(i,:)= [diff(d2vdt2(i,:),1) 0];                
    g(i,setNZ)=d2vdt2(i,setNZ)./dvdt(i,setNZ);                                                      % Metric for which the maximum indicate threshold; method 1
    hh(i,setNZ)=(d3vdt3(i,setNZ).*dvdt(i,setNZ)-d2vdt2(i,setNZ).^2)./(dvdt(i,setNZ).^3);            % Metric for which the maximum indicate threshold; method 3
    %SR2 setNZwin2=find((setNZ<win) & (setNZ>(win-win2)))+win-win2 ;         % points in win2 before peak of spike
    thres_coord(i)=min(find(g(i,win-win2:win) == max(g(i,win-win2:win))))+win-win2-1 ;
    thres_coord_hh(i)=min(find(hh(i,win-win2:win) == max(hh(i,win-win2:win))))+win-win2-1 ;
    thresvalue(i)=spiketrace(i,thres_coord(i));
    thresvalue_hh(i)=spiketrace(i,thres_coord_hh(i));
    thres_coords_orginal(i)=spikecords(i)-win+thres_coord(i) ;
    thres_coords_orginal_hh(i)=spikecords(i)-win+thres_coord_hh(i) ;
    
    Am_list(1, i) = Am - thresvalue(i);
    Am_list(2, i) = spikecords(i)-win+Am_Coord;
    a = find(spiketrace(i,:) >= (Am_list(1, i)/2));
    Am50_list(1,i) = max(a) - min(a);  Am50_list(2,i) = Am_list(2, i);
    [maxInCurr, RMCoord] = max(dvdt(i,:)); 
    RiseMax_list(1, i) = spiketrace(i,RMCoord); RiseMax_list(2, i) = spikecords(i)-win+RMCoord; 
    MaxInCurr_list(1, i) = maxInCurr; MaxInCurr_list(2, i) = RiseMax_list(2, i);
    [maxOutCurr,DMCoord] = min(dvdt(i,:)); 
    DecayMax_list(1, i) = spiketrace(i,DMCoord); DecayMax_list(2, i) = spikecords(i)-win+DMCoord;  
    MaxOutCurr_list(1, i) = maxOutCurr; MaxOutCurr_list(2, i) = DecayMax_list(2, i);
    ChargeIn_list(1,i) = trapz(dvdt(find(dvdt(i,:)>0))); ChargeIn_list(2,i) = Am_list(2, i);
    ChargeOut_list(1,i) = trapz(dvdt(find(dvdt(i,:)<0))); ChargeOut_list(2,i) = Am_list(2, i);
    %AP_power_list(i) = trapz(spiketrace((dvdt(i,:)>0)), dvdt((dvdt(i,:)>0))) + trapz(spiketrace((dvdt(i,:)<0)), dvdt((dvdt(i,:)<0)));
    AP_power_list(i) = trapz(spiketrace(i,:), dvdt(i,:));
    
    %     figure(204)
    %     plot(spiketrace(i,:)); hold on;
    %     plot(base(setNZ),dvdt(i,setNZ),'k');
    %     plot(base(setNZ),d2vdt2(i,setNZ),'b');
    %     plot(base(setNZ),d3vdt3(i,setNZ),'m');
    %     %SR2 plot(base(setNZwin2),10*d2vdt2(i,setNZwin2)./dvdt(i,setNZwin2),'g')
    %     hold off ; ylim([-50 30])
    %     
    %     figure(202)
    %     plot(spiketrace(i,:));hold on
    %     plot(50*g(i,:),'r');
    %     plot(50*hh(i,:),'b');
    %     plot(thres_coord(i), spiketrace(i,thres_coord(i)) ,'or');
    %     plot(thres_coord_hh(i), spiketrace(i,thres_coord_hh(i)) ,'ob');
    %     plot(setNZ, 10, 'ob'); 
    %     ylim([-50 50])
    % plot(setNZwin2, 5, 'ok')
    hold off ;
    clear setNZ; 
    % clear setNZwin2
           %drawnow
end

mean_thresh=mean(thresvalue) ;   sd_thresh=std(thresvalue) ;
mean_thresh_hh=mean(thresvalue_hh) ;   sd_thresh_hh=std(thresvalue_hh) ;
%%
%     figure(203)
%     subplot(131); plot(spiketrace','color',[.7 .7 .7]); hold on
%     for i=1:N
%         plot(thres_coord(i),spiketrace(i,thres_coord(i)),'.r')
%         plot(thres_coord_hh(i),spiketrace(i,thres_coord_hh(i)),'.b')
%     end
%     plot(mean(spiketrace,1)','b');hold off; grid on
%     xlabel('datapoints')
%     subplot(132);plot(spiketrace',dvdt','color',[.7 .7 .7]);hold on
%     plot(mean(spiketrace,1)',mean(dvdt,1)','r');hold off; grid on
%     xlabel('V'); ylabel('dV/dt')
%     subplot(133);hist(thresvalue); xlabel('thresholds');
%     
%     figure(205)
%     plot(traceD);hold on
%     plot(thres_coords_orginal,trace2(thres_coords_orginal),'.r')
%     plot([1 length(traceD)], [mean_thresh mean_thresh],'k')
%     plot([1 length(traceD)], [mean_thresh-sd_thresh mean_thresh-sd_thresh],'--k')
%     plot([1 length(traceD)], [mean_thresh+sd_thresh mean_thresh+sd_thresh],'--k')
%     hold off

% figure(204)
%     for i=1:N
%         plot(spiketrace(i,:),dvdt(i,:),'color',[colorline_r(i), colorline_g(i), colorline_b(i)]); hold on
%     end
% xlabel('V, mV'); ylabel('dV/dt, mV/msec')
% hold off;
% figure(205)
% hist(thresvalue); xlabel('thresholds, mV');
%--------------------------------------------------------------------------
 
HeadHorz = {'Am' 'Am50' 'RiseMax' 'DecayMax' 'MaxInCurr' 'MaxOutCurr' 'ChargeIn' 'ChargeOut'};
AP_Data_Table = table(Am_list', Am50_list', RiseMax_list', DecayMax_list', MaxInCurr_list', MaxOutCurr_list', ChargeIn_list', ChargeOut_list', 'VariableNames', HeadHorz)

clear('traceD')
end

function [N ,out1] = spike_times(traceD,threshold1)

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

 gim=traceD;

    clear('set_crossgi')
    set_crossgi=find(gim > threshold1)  ;  % setting the threshold
    clear('index_shift_neggi');clear('index_shift_pos');

if isempty(set_crossgi) < 1     % This to make sure there is a spike otherwise the code below gives problems. There is an empty else statement below.

    clear('set_cross_plusgi');clear('set_cross_minus')
    index_shift_posgi(1)=min(set_crossgi);
    index_shift_neggi(length(set_crossgi))=max(set_crossgi);

    % SR the condition was changed from 1 to 2 for better sustainability...
    for i=1:length(set_crossgi)-1
     if set_crossgi(i+1) > set_crossgi(i)+2; 
     index_shift_posgi(i+1)=i;
     index_shift_neggi(i)=i;
     end
    end

    %Identifying up and down slopes:

    set_cross_plusgi=  set_crossgi(find(index_shift_posgi));   % find(x) returns nonzero arguments.
    set_cross_minusgi=  set_crossgi(find(index_shift_neggi));   % find(x) returns nonzero arguments.
    set_cross_minusgi(length(set_cross_plusgi))= set_crossgi(end);
    nspikes= length(set_cross_plusgi); % Number of pulses, i.e. number of windows.

    %% Getting the spike coords
    %spikemax = zeros(1,nspikes)
    for i=1: nspikes
            spikemax(i)=min(find(gim(set_cross_plusgi(i):set_cross_minusgi(i)) == max(gim(set_cross_plusgi(i):set_cross_minusgi(i))))) +set_cross_plusgi(i)-1;
    end

else
    spikemax=[];
    display('no spikes in trace')
end 
N=length(spikemax)
out1=spikemax;
end