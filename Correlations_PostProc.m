thresvalue_list_noZero = thresvalue_list(thresvalue_list~=0);
temp_list_noZero = temp_list(temp_list~=0);

[nyy, myy]= size(thres_coords_orgina_list);
thres_coords_orgina_list_OffBase = zeros(nyy, myy); % temes of AP Th parameters
HeatRate_list = zeros(myy, 2); DeltaTime = zeros(nyy, myy);
for yy=1:myy
    thres_coords_orgina_list_OffBase(:,yy) = (thres_coords_orgina_list(:,yy)-thres_coords_orgina_list(1,yy));
    DeltaTime(2:end,yy) = thres_coords_orgina_list(2:end,yy)-thres_coords_orgina_list(1:end-1,yy);
    Xx = thres_coords_orgina_list_OffBase(:,yy);
    X = Xx(Xx>=0).*(dt/10^6); %time
    Yy = temp_list(:,yy); 
    Y= Yy(Yy>0);% temperature
    B = horzcat(ones(length(X),1), X)\Y; % a+bX = Y ==> B = [a, b]
    if B(2)<0
        B(2)=0;
    end;
    HeatRate_list(yy,:) = B';
end;
DeltaTime(1,:) = DeltaTime(2,:);
thres_coords_orgina_list_OffBase = thres_coords_orgina_list_OffBase.*(dt/10^6);
thres_coords_orgina_list_noZero = thres_coords_orgina_list_OffBase(thres_coords_orgina_list_OffBase>=0); %Time sec
DeltaTime_No_Zero = DeltaTime(DeltaTime>0).*(dt/10^6); %sec
Freq = DeltaTime_No_Zero.^(-1); %Hz


[CorrListTH, CorrProbTH, meanCorrTH] = AverCorr(thresvalue_list, temp_list);
[RlalaTH, PlalaTH] = corrcoef(thresvalue_list_noZero, temp_list_noZero);

A = table2array(AP_Data_Table_Big(:,2:7));% Am, Am50, VRise max, VDecay max, MaxRise, MaxDecay
A = A(:,1:2:12);
[R_Am, P_Am] = corrcoef(A(:,1), temp_list_noZero);
[R_Am50, P_Am50] = corrcoef(A(:,2), temp_list_noZero);
[R_VRise, P_VRise] = corrcoef(A(:,3), temp_list_noZero);
[R_VDec, P_VDec] = corrcoef(A(:,4), temp_list_noZero);
[R_RiseMax, P_RiseMax] = corrcoef(A(:,5), temp_list_noZero);
[R_DecMax, P_DecMax] = corrcoef(A(:,6), temp_list_noZero);

filename = 'C:\MathLab_test_folder\MMW UWA uploads\Hi_Ca\CorrelationReport_MMW_Hi_Ca.xlsx';

xlswrite(filename,'R','Correlation','A1');
xlswrite(filename,R_Am,'Correlation','B1');
xlswrite(filename,'P','Correlation','D1');
xlswrite(filename,P_Am,'Correlation','E1');

xlswrite(filename,'R','Correlation','A3');
xlswrite(filename,R_Am50,'Correlation','B3');
xlswrite(filename,'P','Correlation','D3');
xlswrite(filename,P_Am50,'Correlation','E3');

xlswrite(filename,'R','Correlation','A5');
xlswrite(filename,R_VRise,'Correlation','B5');
xlswrite(filename,'P','Correlation','D5');
xlswrite(filename,P_VRise,'Correlation','E5');

xlswrite(filename,'R','Correlation','A7');
xlswrite(filename,R_VDec,'Correlation','B7');
xlswrite(filename,'P','Correlation','D7');
xlswrite(filename,P_VDec,'Correlation','E7');

xlswrite(filename,'R','Correlation','A9');
xlswrite(filename,R_RiseMax,'Correlation','B9');
xlswrite(filename,'P','Correlation','D9');
xlswrite(filename,P_RiseMax,'Correlation','E9');

xlswrite(filename,'R','Correlation','A11');
xlswrite(filename,R_DecMax,'Correlation','B11');
xlswrite(filename,'P','Correlation','D11');
xlswrite(filename,P_DecMax,'Correlation','E11');

HeadHorzN = {'Am',' ', 'Am50', ' ', 'RiseMax', ' ', 'DecayMax', ' ', 'MaxInCurr', ' ', 'MaxOutCurr'}';
xlswrite(filename,HeadHorzN,'Correlation','H1');

HeadHorzN = {'R_Th', ' ', 'P_Th', ' ', 'R_Th_list -------- P_Th_list'}';
xlswrite(filename,HeadHorzN,'CorrTH','A1');
xlswrite(filename,RlalaTH,'CorrTH','B1');
%xlswrite(filename,'P_Th','CorrTH','A3');
xlswrite(filename,PlalaTH,'CorrTH','B3');
%xlswrite(filename,'R_Th_list','CorrTH','A5');
xlswrite(filename,CorrListTH','CorrTH','A6');
%xlswrite(filename,'P_Th_list','CorrTH','A7');
xlswrite(filename,CorrProbTH','CorrTH','B6');

%%
thresvalue_list_noZero = thresvalue_list(thresvalue_list~=0);
temp_list_noZero = temp_list(temp_list~=0);

[nyy, myy]= size(thres_coords_orgina_list);
thres_coords_orgina_list_OffBase = zeros(nyy, myy); % temes of AP Th parameters
HeatRate_list = zeros(myy, 2); DeltaTime = zeros(nyy, myy);
for yy=1:myy
    thres_coords_orgina_list_OffBase(:,yy) = (thres_coords_orgina_list(:,yy)-thres_coords_orgina_list(1,yy));
    DeltaTime(2:end,yy) = thres_coords_orgina_list(2:end,yy)-thres_coords_orgina_list(1:end-1,yy);
    Xx = thres_coords_orgina_list_OffBase(:,yy);
    X = Xx(Xx>=0).*(dt/10^6); %time
    Yy = temp_list(:,yy); 
    Y= Yy(Yy>0);% temperature
    B = horzcat(ones(length(X),1), X)\Y; % a+bX = Y ==> B = [a, b]
    if B(2)<0
        B(2)=0;
    end;
    HeatRate_list(yy,:) = B';
end;
DeltaTime(1,:) = DeltaTime(2,:);
thres_coords_orgina_list_OffBase = thres_coords_orgina_list_OffBase.*(dt/10^6);
thres_coords_orgina_list_noZero = thres_coords_orgina_list_OffBase(thres_coords_orgina_list_OffBase>=0); %Time sec
DeltaTime_No_Zero = DeltaTime(DeltaTime>0).*(dt/10^6); %sec
Freq = DeltaTime_No_Zero.^(-1); %Hz

[CorrListTH, CorrProbTH, meanCorrTH] = AverCorr(thresvalue_list, temp_list);
[RlalaTH, PlalaTH] = corrcoef(thresvalue_list_noZero, temp_list_noZero);

A = table2array(AP_Data_Table_Big(:,2:7));% Am, Am50, VRise max, VDecay max, MaxRise, MaxDecay
A = A(:,1:2:12);
[R_Am, P_Am] = corrcoef(A(:,1), temp_list_noZero);
[R_Am50, P_Am50] = corrcoef(A(:,2), temp_list_noZero);
[R_VRise, P_VRise] = corrcoef(A(:,3), temp_list_noZero);
[R_VDec, P_VDec] = corrcoef(A(:,4), temp_list_noZero);
[R_RiseMax, P_RiseMax] = corrcoef(A(:,5), temp_list_noZero);
[R_DecMax, P_DecMax] = corrcoef(A(:,6), temp_list_noZero);

filename = 'C:\MathLab_test_folder\MMW UWA uploads\Upload\MMW\Test processing\CorelationReport_MMW.xlsx';

xlswrite(filename,'R','Correlation','A1');
xlswrite(filename,R_Am,'Correlation','B1');
xlswrite(filename,'P','Correlation','D1');
xlswrite(filename,P_Am,'Correlation','E1');

xlswrite(filename,'R','Correlation','A3');
xlswrite(filename,R_Am50,'Correlation','B3');
xlswrite(filename,'P','Correlation','D3');
xlswrite(filename,P_Am50,'Correlation','E3');

xlswrite(filename,'R','Correlation','A5');
xlswrite(filename,R_VRise,'Correlation','B5');
xlswrite(filename,'P','Correlation','D5');
xlswrite(filename,P_VRise,'Correlation','E5');

xlswrite(filename,'R','Correlation','A7');
xlswrite(filename,R_VDec,'Correlation','B7');
xlswrite(filename,'P','Correlation','D7');
xlswrite(filename,P_VDec,'Correlation','E7');

xlswrite(filename,'R','Correlation','A9');
xlswrite(filename,R_RiseMax,'Correlation','B9');
xlswrite(filename,'P','Correlation','D9');
xlswrite(filename,P_RiseMax,'Correlation','E9');

xlswrite(filename,'R','Correlation','A11');
xlswrite(filename,R_DecMax,'Correlation','B11');
xlswrite(filename,'P','Correlation','D11');
xlswrite(filename,P_DecMax,'Correlation','E11');

HeadHorzN = {'Am',' ', 'Am50', ' ', 'RiseMax', ' ', 'DecayMax', ' ', 'MaxInCurr', ' ', 'MaxOutCurr'}';
xlswrite(filename,HeadHorzN,'Correlation','H1');

HeadHorzN = {'R_Th', ' ', 'P_Th', ' ', 'R_Th_list -------- P_Th_list'}';
xlswrite(filename,HeadHorzN,'CorrTH','A1');
xlswrite(filename,RlalaTH,'CorrTH','B1');
%xlswrite(filename,'P_Th','CorrTH','A3');
xlswrite(filename,PlalaTH,'CorrTH','B3');
%xlswrite(filename,'R_Th_list','CorrTH','A5');
xlswrite(filename,CorrListTH','CorrTH','A6');
%xlswrite(filename,'P_Th_list','CorrTH','A7');
xlswrite(filename,CorrProbTH','CorrTH','B6');
%%
thresvalue_list_noZero = thresvalue_list(thresvalue_list~=0);
temp_list_noZero = temp_list(temp_list~=0);

[nyy, myy]= size(thres_coords_orgina_list);
thres_coords_orgina_list_OffBase = zeros(nyy, myy); % temes of AP Th parameters
HeatRate_list = zeros(myy, 2); DeltaTime = zeros(nyy, myy);
for yy=1:myy
    Yy = temp_list(:,yy); 
    Y= Yy(Yy>0);% temperature
    thres_coords_orgina_list_OffBase(:,yy) = (thres_coords_orgina_list(:,yy)-thres_coords_orgina_list(1,yy));
    DeltaTime(2:end,yy) = thres_coords_orgina_list(2:end,yy)-thres_coords_orgina_list(1:end-1,yy);
    
    Xx = thres_coords_orgina_list_OffBase(:,yy);
    %Xx = thres_coords_orgina_list_OffBase(1:length(Y),yy);
    
    X = Xx(Xx>=0).*(dt/10^6); %time
        B = horzcat(ones(length(X),1), X)\Y; % a+bX = Y ==> B = [a, b]
    if B(2)<0
        B(2)=0;
    end;
    HeatRate_list(yy,:) = B';
end;
DeltaTime(1,:) = DeltaTime(2,:);
thres_coords_orgina_list_OffBase = thres_coords_orgina_list_OffBase.*(dt/10^6);
thres_coords_orgina_list_noZero = thres_coords_orgina_list_OffBase(thres_coords_orgina_list_OffBase>=0); %Time sec
DeltaTime_No_Zero = DeltaTime(DeltaTime>0).*(dt/10^6); %sec
Freq = DeltaTime_No_Zero.^(-1); %Hz

[CorrListTH, CorrProbTH, meanCorrTH] = AverCorr(thresvalue_list, temp_list);
[RlalaTH, PlalaTH] = corrcoef(thresvalue_list_noZero, temp_list_noZero);

A = table2array(AP_Data_Table_Big(:,2:7));% Am, Am50, VRise max, VDecay max, MaxRise, MaxDecay
A = A(:,1:2:12);
[R_Am, P_Am] = corrcoef(A(:,1), temp_list_noZero);
[R_Am50, P_Am50] = corrcoef(A(:,2), temp_list_noZero);
[R_VRise, P_VRise] = corrcoef(A(:,3), temp_list_noZero);
[R_VDec, P_VDec] = corrcoef(A(:,4), temp_list_noZero);
[R_RiseMax, P_RiseMax] = corrcoef(A(:,5), temp_list_noZero);
[R_DecMax, P_DecMax] = corrcoef(A(:,6), temp_list_noZero);

%filename = 'C:\MathLab_test_folder\MMW UWA uploads\Upload\Thermal\Test processing\CorelationReport_Thermal_A.xlsx';
%filename = 'C:\MathLab_test_folder\MMW UWA uploads\Upload\Thermal\Test processing\CorelationReport_Thermal_B.xlsx';
%filename = 'C:\MathLab_test_folder\MMW UWA uploads\Upload\SB366791\Test_processing\\CorelationReport_SB366791.xlsx';
filename = 'C:\MathLab_test_folder\MMW UWA uploads\Upload\Nm\Test_processing\CorelationReport_Nm_test.xlsx';

xlswrite(filename,'R','Correlation','A1');
xlswrite(filename,R_Am,'Correlation','B1');
xlswrite(filename,'P','Correlation','D1');
xlswrite(filename,P_Am,'Correlation','E1');

xlswrite(filename,'R','Correlation','A3');
xlswrite(filename,R_Am50,'Correlation','B3');
xlswrite(filename,'P','Correlation','D3');
xlswrite(filename,P_Am50,'Correlation','E3');

xlswrite(filename,'R','Correlation','A5');
xlswrite(filename,R_VRise,'Correlation','B5');
xlswrite(filename,'P','Correlation','D5');
xlswrite(filename,P_VRise,'Correlation','E5');

xlswrite(filename,'R','Correlation','A7');
xlswrite(filename,R_VDec,'Correlation','B7');
xlswrite(filename,'P','Correlation','D7');
xlswrite(filename,P_VDec,'Correlation','E7');

xlswrite(filename,'R','Correlation','A9');
xlswrite(filename,R_RiseMax,'Correlation','B9');
xlswrite(filename,'P','Correlation','D9');
xlswrite(filename,P_RiseMax,'Correlation','E9');

xlswrite(filename,'R','Correlation','A11');
xlswrite(filename,R_DecMax,'Correlation','B11');
xlswrite(filename,'P','Correlation','D11');
xlswrite(filename,P_DecMax,'Correlation','E11');

HeadHorzN = {'Am',' ', 'Am50', ' ', 'RiseMax', ' ', 'DecayMax', ' ', 'MaxInCurr', ' ', 'MaxOutCurr'}';
xlswrite(filename,HeadHorzN,'Correlation','H1');

HeadHorzN = {'R_Th', ' ', 'P_Th', ' ', 'R_Th_list -------- P_Th_list'}';
xlswrite(filename,HeadHorzN,'CorrTH','A1');
xlswrite(filename,RlalaTH,'CorrTH','B1');
%xlswrite(filename,'P_Th','CorrTH','A3');
xlswrite(filename,PlalaTH,'CorrTH','B3');
%xlswrite(filename,'R_Th_list','CorrTH','A5');
xlswrite(filename,CorrListTH','CorrTH','A6');
%xlswrite(filename,'P_Th_list','CorrTH','A7');
xlswrite(filename,CorrProbTH','CorrTH','B6');

five = thresvalue_list(1:5,:);
five = five(:);
%%
thresvalue_list_noZero = thresvalue_list(thresvalue_list~=0);
temp_list_noZero = temp_list(temp_list~=0);

[nyy, myy]= size(thres_coords_orgina_list);
thres_coords_orgina_list_OffBase = zeros(nyy, myy); % temes of AP Th parameters
HeatRate_list = zeros(myy, 2); DeltaTime = zeros(nyy, myy);
for yy=1:myy
    thres_coords_orgina_list_OffBase(:,yy) = (thres_coords_orgina_list(:,yy)-thres_coords_orgina_list(1,yy));
    DeltaTime(2:end,yy) = thres_coords_orgina_list(2:end,yy)-thres_coords_orgina_list(1:end-1,yy);
    Xx = thres_coords_orgina_list_OffBase(:,yy);
    X = Xx(Xx>=0).*(dt/10^6); %time
    Yy = temp_list(:,yy); 
    Y= Yy(Yy>0);% temperature
    B = horzcat(ones(length(X),1), X)\Y; % a+bX = Y ==> B = [a, b]
    if B(2)<0
        B(2)=0;
    end;
    HeatRate_list(yy,:) = B';
end;
DeltaTime(1,:) = DeltaTime(2,:);
thres_coords_orgina_list_OffBase = thres_coords_orgina_list_OffBase.*(dt/10^6);
thres_coords_orgina_list_noZero = thres_coords_orgina_list_OffBase(thres_coords_orgina_list_OffBase>=0); %Time sec
DeltaTime_No_Zero = DeltaTime(DeltaTime>0).*(dt/10^6); %sec
Freq = DeltaTime_No_Zero.^(-1); %Hz

[CorrListTH, CorrProbTH, meanCorrTH] = AverCorr(thresvalue_list, temp_list);
[RlalaTH, PlalaTH] = corrcoef(thresvalue_list_noZero, temp_list_noZero);

A = table2array(AP_Data_Table_Big(:,2:7));% Am, Am50, VRise max, VDecay max, MaxRise, MaxDecay
A = A(:,1:2:12);
[R_Am, P_Am] = corrcoef(A(:,1), temp_list_noZero);
[R_Am50, P_Am50] = corrcoef(A(:,2), temp_list_noZero);
[R_VRise, P_VRise] = corrcoef(A(:,3), temp_list_noZero);
[R_VDec, P_VDec] = corrcoef(A(:,4), temp_list_noZero);
[R_RiseMax, P_RiseMax] = corrcoef(A(:,5), temp_list_noZero);
[R_DecMax, P_DecMax] = corrcoef(A(:,6), temp_list_noZero);



xlswrite(filename,'R','Correlation','A1');
xlswrite(filename,R_Am,'Correlation','B1');
xlswrite(filename,'P','Correlation','D1');
xlswrite(filename,P_Am,'Correlation','E1');

xlswrite(filename,'R','Correlation','A3');
xlswrite(filename,R_Am50,'Correlation','B3');
xlswrite(filename,'P','Correlation','D3');
xlswrite(filename,P_Am50,'Correlation','E3');

xlswrite(filename,'R','Correlation','A5');
xlswrite(filename,R_VRise,'Correlation','B5');
xlswrite(filename,'P','Correlation','D5');
xlswrite(filename,P_VRise,'Correlation','E5');

xlswrite(filename,'R','Correlation','A7');
xlswrite(filename,R_VDec,'Correlation','B7');
xlswrite(filename,'P','Correlation','D7');
xlswrite(filename,P_VDec,'Correlation','E7');

xlswrite(filename,'R','Correlation','A9');
xlswrite(filename,R_RiseMax,'Correlation','B9');
xlswrite(filename,'P','Correlation','D9');
xlswrite(filename,P_RiseMax,'Correlation','E9');

xlswrite(filename,'R','Correlation','A11');
xlswrite(filename,R_DecMax,'Correlation','B11');
xlswrite(filename,'P','Correlation','D11');
xlswrite(filename,P_DecMax,'Correlation','E11');

HeadHorzN = {'Am',' ', 'Am50', ' ', 'RiseMax', ' ', 'DecayMax', ' ', 'MaxInCurr', ' ', 'MaxOutCurr'}';
xlswrite(filename,HeadHorzN,'Correlation','H1');

HeadHorzN = {'R_Th', ' ', 'P_Th', ' ', 'R_Th_list -------- P_Th_list'}';
xlswrite(filename,HeadHorzN,'CorrTH','A1');
xlswrite(filename,RlalaTH,'CorrTH','B1');
%xlswrite(filename,'P_Th','CorrTH','A3');
xlswrite(filename,PlalaTH,'CorrTH','B3');
%xlswrite(filename,'R_Th_list','CorrTH','A5');
xlswrite(filename,CorrListTH','CorrTH','A6');
%xlswrite(filename,'P_Th_list','CorrTH','A7');
xlswrite(filename,CorrProbTH','CorrTH','B6');
%%
thresvalue_list_noZero = thresvalue_list(thresvalue_list~=0);
temp_list_noZero = temp_list(temp_list~=0);

[nyy, myy]= size(thres_coords_orgina_list);
thres_coords_orgina_list_OffBase = zeros(nyy, myy); % temes of AP Th parameters
HeatRate_list = zeros(myy, 2); DeltaTime = zeros(nyy, myy);
for yy=1:myy
    thres_coords_orgina_list_OffBase(:,yy) = (thres_coords_orgina_list(:,yy)-thres_coords_orgina_list(1,yy));
    DeltaTime(2:end,yy) = thres_coords_orgina_list(2:end,yy)-thres_coords_orgina_list(1:end-1,yy);
    Xx = thres_coords_orgina_list_OffBase(:,yy);
    X = Xx(Xx>=0).*(dt/10^6); %time
    Yy = temp_list(:,yy); 
    Y= Yy(Yy>0);% temperature
    B = horzcat(ones(length(X),1), X)\Y; % a+bX = Y ==> B = [a, b]
    if B(2)<0
        B(2)=0;
    end;
    HeatRate_list(yy,:) = B';
end;
DeltaTime(1,:) = DeltaTime(2,:);
thres_coords_orgina_list_OffBase = thres_coords_orgina_list_OffBase.*(dt/10^6);
thres_coords_orgina_list_noZero = thres_coords_orgina_list_OffBase(thres_coords_orgina_list_OffBase>=0); %Time sec
DeltaTime_No_Zero = DeltaTime(DeltaTime>0).*(dt/10^6); %sec
Freq = DeltaTime_No_Zero.^(-1); %Hz

[CorrListTH, CorrProbTH, meanCorrTH] = AverCorr(thresvalue_list, temp_list);
[RlalaTH, PlalaTH] = corrcoef(thresvalue_list_noZero, temp_list_noZero);

A = table2array(AP_Data_Table_Big(:,2:7));% Am, Am50, VRise max, VDecay max, MaxRise, MaxDecay
A = A(:,1:2:12);
[R_Am, P_Am] = corrcoef(A(:,1), temp_list_noZero);
[R_Am50, P_Am50] = corrcoef(A(:,2), temp_list_noZero);
[R_VRise, P_VRise] = corrcoef(A(:,3), temp_list_noZero);
[R_VDec, P_VDec] = corrcoef(A(:,4), temp_list_noZero);
[R_RiseMax, P_RiseMax] = corrcoef(A(:,5), temp_list_noZero);
[R_DecMax, P_DecMax] = corrcoef(A(:,6), temp_list_noZero);



xlswrite(filename,'R','Correlation','A1');
xlswrite(filename,R_Am,'Correlation','B1');
xlswrite(filename,'P','Correlation','D1');
xlswrite(filename,P_Am,'Correlation','E1');

xlswrite(filename,'R','Correlation','A3');
xlswrite(filename,R_Am50,'Correlation','B3');
xlswrite(filename,'P','Correlation','D3');
xlswrite(filename,P_Am50,'Correlation','E3');

xlswrite(filename,'R','Correlation','A5');
xlswrite(filename,R_VRise,'Correlation','B5');
xlswrite(filename,'P','Correlation','D5');
xlswrite(filename,P_VRise,'Correlation','E5');

xlswrite(filename,'R','Correlation','A7');
xlswrite(filename,R_VDec,'Correlation','B7');
xlswrite(filename,'P','Correlation','D7');
xlswrite(filename,P_VDec,'Correlation','E7');

xlswrite(filename,'R','Correlation','A9');
xlswrite(filename,R_RiseMax,'Correlation','B9');
xlswrite(filename,'P','Correlation','D9');
xlswrite(filename,P_RiseMax,'Correlation','E9');

xlswrite(filename,'R','Correlation','A11');
xlswrite(filename,R_DecMax,'Correlation','B11');
xlswrite(filename,'P','Correlation','D11');
xlswrite(filename,P_DecMax,'Correlation','E11');

HeadHorzN = {'Am',' ', 'Am50', ' ', 'RiseMax', ' ', 'DecayMax', ' ', 'MaxInCurr', ' ', 'MaxOutCurr'}';
xlswrite(filename,HeadHorzN,'Correlation','H1');

HeadHorzN = {'R_Th', ' ', 'P_Th', ' ', 'R_Th_list -------- P_Th_list'}';
xlswrite(filename,HeadHorzN,'CorrTH','A1');
xlswrite(filename,RlalaTH,'CorrTH','B1');
%xlswrite(filename,'P_Th','CorrTH','A3');
xlswrite(filename,PlalaTH,'CorrTH','B3');
%xlswrite(filename,'R_Th_list','CorrTH','A5');
xlswrite(filename,CorrListTH','CorrTH','A6');
%xlswrite(filename,'P_Th_list','CorrTH','A7');
xlswrite(filename,CorrProbTH','CorrTH','B6');
%%


%%


