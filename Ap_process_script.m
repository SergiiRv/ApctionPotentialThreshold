
%%
clear
passname ='put here the forlder where all processed files are placed';
        
    headers_horz = {'1.abf', '2.abf', '3.abf', 'and so-on.abf',};
        
nnn = int32(length(headers_horz));
fnlist = {' '};
for j = 1:nnn
    fnlist(j) = strcat(passname,headers_horz(j));
end
       
%input
nn = int32(length(fnlist));
        
% start and end of the trace to process (unit: s) 
startlist = [ , , , ];
endlist =  [ , , , ];

       
%% output
mean_thresh_list = zeros(nn,1); mean_thresh_list_hh = zeros(nn,1);
sd_thresh_list = zeros(nn,1); sd_thresh_list_hh = zeros(nn,1);
N_list = [1:nn]'.*0;
thresvalue_list = zeros(500,nn); thresvalue_list_hh = zeros(500,nn);
thres_coords_orgina_list = zeros(500,nn); thres_coords_orgina_list_hh = zeros(500,nn);
temp_list = zeros(500,nn); temp_list_hh = zeros(500,nn); 
        
channel = {'IN 0'};
channel_temp = {'IN 6'};
AP_HeadVert = {}; APnewBig = []; AP_HeadVertBig =[];
AP_HeadHorz = {'File_name_processed' 'Am' 'Am_temp' 'Am50' 'Am50_temp' 'RiseMax' 'RiseMax_temp' 'DecayMax' 'DecayMax_temp' 'MaxInCurr' 'MaxInCurr_temp' 'MaxOutCurr' 'MaxOutCurr_temp' 'ChargeIn' 'ChargeIn_temp' 'ChargeOut' 'ChargeOut_temp'};
[~, secS] = size(AP_HeadHorz);
regCellBig = cell(nn, (secS-1)/2); CorelLineBig = cell(nn, (secS-1)/2); CorelLineBigTH =  cell(nn, 2); SCorelLineBig = cell(nn, 1);
temp_Grid = [25:0.001:50]; TH1_Grid_list = zeros(length(temp_Grid), nn); TH2_Grid_list = zeros(length(temp_Grid), nn); Param_Grid_3d = [];
temp_list_UP = zeros(500,nn); temp_list_DW= zeros(500,nn); temp_list_hh_UP = zeros(500,nn); temp_list_hh_DW= zeros(500,nn);
thresvalue_list_UP = zeros(500,nn); thresvalue_list_DW = zeros(500,nn); thresvalue_list_hh_UP = zeros(500,nn); thresvalue_list_hh_DW = zeros(500,nn);

%% 
%% Calculus
for  i =1:nn
   [data, dt, info] = abfload(char(fnlist(i)),'start',startlist(i),'stop',endlist(i),'channels',channel);
   %[data_temp] = abfload(char(fnlist(i)),'start',startlist(i),'stop',endlist(i),'channels',channel_temp);
   [mean_thresh sd_thresh mean_thresh_hh sd_thresh_hh Np thresvalue thres_coords_orginal thresvalue_hh thres_coords_orginal_hh AP_Data_Table_small]=Spike_threshold_PS(data,dt/1000000,45000/dt,20000/dt, 500, 0.01);
   
   mean_thresh_list(i)= mean_thresh;
   sd_thresh_list(i) = sd_thresh;
   mean_thresh_list_hh(i)= mean_thresh_hh;
   sd_thresh_list_hh(i) = sd_thresh_hh;
   N_list(i) = Np;
   thresvalue_list(1:length(thresvalue'),i) = thresvalue';
   thres_coords_orgina_list(1:length(thres_coords_orginal'),i) = thres_coords_orginal';
   thresvalue_list_hh(1:length(thresvalue_hh'),i) = thresvalue_hh';
   thres_coords_orgina_list_hh(1:length(thres_coords_orginal_hh'),i) = thres_coords_orginal_hh';

end

filename = 'C:\...\finalReport.xlsx';

writetable(AP_Data_Table_small,filename,'Sheet','Last file AP parameters','Range','B2');
xlswrite(filename,headers_horz,'Th_values_hh','A1');
xlswrite(filename,thresvalue_list_hh,'Th_values_hh','A2');
xlswrite(filename,headers_horz,'Th_times_hh','A1');
xlswrite(filename,thres_coords_orgina_list_hh,'Th_times_hh','A2');

%%
'end of processing'
%%

