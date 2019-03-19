clear all
clc
close all
fs_org = 900;
A2D_2_Volt = 5/(2^16);

c = 3e8;
f = 24.125e9;
Landa = c/f;
Coef = 4*pi/Landa;


Zone_Index = 9;

disp(['The subject is in Zone #',num2str(Zone_Index)])

% Reading the radar signal
Data_Address = 'C:\Users\Hamidreza\Desktop\response ISAR\Isar papers\Isar papers\Confidenceinterval_paper\Data\';
Subject = {'IG_'};  %,,'IG_', 'IN_', 'XZ_'
Activity = {'NB_'}; %,'MV_'
Posture_Location = {'LY_F_'}; % 'LY_F_','SI_F_','ST_C_'
for sub =1:1
    for act = 1:1
        for PL =1:1
            Filename = [Data_Address,Subject{sub},Activity{act},Posture_Location{PL},'A02_FR'];
            Data =csvread([Filename,'.csv'],1,0);
            
            % Reading Belt signal
            Data_Belt = csvread([Filename,'_BELT.csv'],1,0);
            Data_Belt(Data_Belt>10000) = 0;
            Data_Br = Data_Belt(:,2);
            % Data_Br = tsmovavg(Data_Br', 's',1000);
            
            Data_Hr = Data_Belt(:,1);
            fs_br = 680;%size(Data_Belt,1)/180;
            
            
            %%%%%% separating and downsampling the zone of interest
            Signal = downsample(Data(:,Zone_Index),10);
            fs = fs_org/10;
            t = (0:length(Signal)-1)/fs;
            
            Window_Length = 15;  %sec
            Overlap = 10; %sec
            K_end = floor(length(Signal)/(Overlap*fs))-1;  % number of segments
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for K =1:K_end
                
                T1(K) = (K-1)*Overlap;  %start time of segment in second
                
                
                T2(K) = T1(K)+Window_Length;  %stop time of segment in second
                
                disp(char(10))
                disp(['second ',num2str(T1(K)), ' to second ', num2str(T2(K))]);
                
                
                S1=T1(K)*fs+1;  %start time of segment in sample
                S1_br=(T1(K)-0)*fs_br+1;  %start time of segment in sample
                if K == K_end
                    S2 = length(Signal);
                    S2_br = length(Data_Br);
                else
                    S2=T2(K)*fs;   %stop time of segment in sample
                    S2_br = (T2(K)-0)*fs_br;
                end
                Signal_K = Signal(S1:S2);  %Segment of ineterest
                t_K = (0:length(Signal_K)-1)/fs;
                
                Signal_Br_RIP = Data_Br(S1_br:S2_br);
                Signal_Br_RIP= ( Signal_Br_RIP-mean( Signal_Br_RIP))/std( Signal_Br_RIP);
                
                Signal_K = detrend((Signal_K-mean(Signal_K)));  %deterend the signal
                
                
                
                [b,a] = butter(4,0.2/fs,'high'); % 0.1Hz highpass
                Signal_K = filter(b,a,detrend(Signal_K));
                [b,a] = butter(5,10/fs,'low'); % 5Hz lowpass
                Signal_K = filter(b,a,(Signal_K));
                
                Signal_K_t = Signal_K;  % keep a copy of bandpass filetered signal
                
                
                [x_br_fft(K),amp,phs]= Chirp_est( Signal_K,fs ,0.1,0.4,0);  % Estiamtion of breathing from radar using Chirp transform
                
                x_br_fft(K) = 60*x_br_fft(K);  %bpm
                
                %find the main component
                disp('Extracting first component ...')
                [WT_br,freq_br,wopt_br]=wft(Signal_K,fs,'f0',2,'fmin',0.1,'fmax',0.4,'Display','off','Plot','off');
                [tfsupp_br]=ecurve(WT_br,freq_br,wopt_br,'Display','off','Plot','off');
                [iamp1,iphi1,ifreq1] = rectfr(tfsupp_br,WT_br,freq_br,wopt_br);
                %construct the component in time
                recon_signal = iamp1.*cos(iphi1);
                %subtract it from the original signal
                Res = Signal_K'-recon_signal;
                %construct the breathing model
                model = recon_signal;
              
                
                x_br_TFR = mean(ifreq1);
               
                
                Breathing = model;
                
                
                %     [x_br(K),amp,phs]= Chirp_est( Breathing,fs ,0.1,0.4,0);  % Estiamtion of breathing using chirp transform
                
                RESID = Res;  %residuals to construct bootstrap samples
                RESID = RESID - mean(RESID);
                
                s = length(RESID);
                B = 100;   % number of constructed boostrap signals
                
                surr = createsurr(RESID,B);
                
                %%%%%%%%%%%  Bootstrap algorithm
                
                
                
                alpha = 0.05;   %  95% confidence interval
                for b = 1:B
                    
                    %%%%%%%%% Constructing a bootstrap sample %%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    x_star =Breathing +surr(b,:);  % constructed bootstrap sample
                    
                    [WT_br_star,freq_br_star,wopt_br_star]=wft( x_star,fs,'f0',2,'fmin',0.1,'fmax',0.4,'Display','off','Plot','off');
                    [tfsupp_br_star]=ecurve(WT_br_star,freq_br_star,wopt_br_star,'Display','off','Plot','off');
                    [iamp1_star,iphi1_star,ifreq1_star] = rectfr(tfsupp_br_star,WT_br_star,freq_br_star,wopt_br_star);
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %studentized statistic measure of this bootstrap sample
                    T_star(b) = (mean(ifreq1_star)-x_br_TFR)/std(ifreq1_star);
                    
                end
                
                T = sort(T_star);
                
                
                U = ceil(B- (B*alpha/2)+1);  %Upper T index
                L = ceil(B*alpha/2);  %Lower T index
                
                
                Up_est(K) = (T(U)*std(ifreq1)+mean(ifreq1))*60;  %% Upper limit estimation
                L_est(K) = (T(L)*std(ifreq1)+mean(ifreq1))*60;   %% Lower limit estimation
                Actuall_Est(K) = mean(ifreq1)*60; %Block bootstrap estimation
                %     [x_br_ref(K),amp,phs]= Chirp_est( Signal_Br_RIP,fs_br ,0.1,0.4,0);  % Estiamtion of breathing from belt
                [WT_br_ref,freq_br_ref,wopt_br_ref]=wft(downsample(Signal_Br_RIP,10),fs_br/10,'f0',1,'fmin',0.1,'fmax',0.4,'Display','off','Plot','off');
                [tfsupp_br_ref]=ecurve(WT_br_ref,freq_br_ref,wopt_br_ref,'Display','off','Plot','off');
                [iamp1_ref,iphi1_ref,ifreq1_ref] = rectfr(tfsupp_br_ref,WT_br_ref,freq_br_ref,wopt_br_ref);
                x_br_ref(K) = mean(ifreq1_ref)*60; %bpm
                disp(char(10))
                disp(['FFT_Estimation = ',num2str(x_br_fft(K)) ,'bpm ----- TFR_Estimation = ',num2str(Actuall_Est(K)) ,'bpm -----Refrence =' ,num2str( x_br_ref(K)), ' bpm'])
                disp([' 95% CI = (',num2str(L_est(K)),',',num2str(Up_est(K)),') bpm'])
                
            end
            FFTT = mean(x_br_fft)
            Estimated = mean(Actuall_Est)
            Actual =mean( x_br_ref)
            % calculating relative errors
            E_fft = 100*abs(x_br_fft-x_br_ref)./x_br_ref;
            E_TFR = 100*abs(Actuall_Est-x_br_ref)./x_br_ref;
            E_Up = 100*abs(Up_est-x_br_ref)./x_br_ref;
            E_Low = 100*abs(L_est-x_br_ref)./x_br_ref;
            
            % Saving th eresults
            Result_Address = 'C:\Users\Hamidreza\Desktop\response ISAR\Isar papers\Isar papers\Confidenceinterval_paper\Results\';
            Result_Filename = [Result_Address,Subject{sub},Activity{act},Posture_Location{PL},'A02_FR_result.mat'];
            
            save(Result_Filename, 'Up_est','L_est','Actuall_Est','x_br_ref','x_br_fft','E_fft','E_TFR','E_Up','E_Low');
            
            
            mean(E_fft)
            mean(E_TFR)
            mean(E_Low)
            mean(E_Up)
            
            
            std(E_fft)
            std(E_TFR)
            std(E_Low)
            std(E_Up)
            
        end
    end
end
