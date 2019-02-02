function [dataoutx, dataoutxbig, roughtrace] = baselinesmooth(data, dt, prewin, postwin)
spikes = 0; xx = 0; cd = 0; yff = 0; 
clearvars coord a1 a2 a3 yff cd a N x spike_coord deltabs;

deltabs = abs(max(data)-min(data));


[spikes,xx] = spike_remove(data, round(0.5*prewin),round(0.5*postwin), (0.6*deltabs+min(data))); %0.995
 %roughtrace = nanfastsmooth(xx, 50000,3,0.8); %
 roughtrace = xx;
% 
% [b1,a1]=butter(2,2*5*dt, 'low'); % low pass butterworth filter
% x = filter(b1,a1,xx);  % Filtered trace to minimize noise
% mm = abs(0.9*mean(xx(1:(round(0.05*length(x))))));
% for i = (1:(round(0.05*length(x))))
%     if abs(x(i))< mm
%         x(i) = xx(i);
%     end
%     
% end
% x(1:(round(0.05*length(x)))) = smooth(x(1:(round(0.05*length(x)))),round(0.005*length(x)));

%% digitisation for peaks search simplification
[N ,spike_coord] = spike_times_cont(data, (0.6*deltabs+min(data))); 
a1 = [1 spike_coord];
a2 = [spike_coord length(data)];
cd = round((a2+a1)./2); 
mids = cd; 

yff = data(cd); 
linfit = fit(cd',yff,'linearinterp');
xbig= linfit([1:length(data)]);
clearvars cd cd2;

n=20; % number of bids per inter spike interval
coord = [];
for i=1:(length(a1))
    a3 = [a1(i):(abs((a2(i)-a1(i))/n)):a2(i)];
    a = coord;
    coord = [a(1:end-1) a3];
end

cd2 = round(coord);
yff = xx(cd2);
linfit = fit(cd2',yff,'linearinterp');
yff2=linfit([1:length(xx)]);

x = yff2;

%% final smoothing
fin = nanfastsmooth(x, 5000,3,0.8); %smooth(x,5000); %------------------------------- 

%figure(2); plot(fin); hold on; plot(data-xx+fin-xbig)
clear coord a1 a2 a3 yff cd a N x spike_coord deltabs;
dataoutx = fin; dataoutxbig =xbig; % nanfastsmooth(xbig, 50000,3,0.8); 
end

    