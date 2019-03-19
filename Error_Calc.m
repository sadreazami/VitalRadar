clear all
clc
close all




% Reading the Results
Result_Address = 'C:\Users\inejadgh\Desktop\Radar\ConfidenceInterval\Confidenceinterval_paper\Results\';
Subject = {'IN_','XZ_','IG_'};
Activity = {'NB_','MV_'};
Posture_Location = {'LY_F_','SI_F_','ST_C_'};
KK = 1;
for act = 1:2
    for PL =1:3
        
        L = 0;
        E_fft_ALL{KK} = [];
        E_TFR_ALL{KK} = [];
        E_Up_ALL{KK} = [];
        E_Low_ALL{KK} = [];
        Width_ALL{KK}= [];
        Outlier_fft{KK}=[];
        Outlier_TFR{KK}=[];
        Outlier_Up{KK}=[];
        Outlier_Low{KK}=[];
        Outlier_Width{KK}=[];
        for sub =1:3
            
            
            Result_Filename = [Result_Address,Subject{sub},Activity{act},Posture_Location{PL},'A02_FR_result.mat'];
            
            load(Result_Filename);
            L = L + length(x_br_ref);
            
            % calculating relative errors
            E_fft_ALL{KK} = [ E_fft_ALL{KK} abs(x_br_fft-x_br_ref)];
            E_TFR_ALL{KK} = [ E_TFR_ALL{KK} abs(Actuall_Est-x_br_ref)];
            E_Up_ALL{KK} = [ E_Up_ALL{KK} abs(Up_est-x_br_ref)];
            E_Low_ALL{KK} = [ E_Low_ALL{KK} abs(L_est-x_br_ref)];
            Width_ALL{KK} = [ Width_ALL{KK} abs(Up_est-L_est)];
            
            
            
        end
        
        [E_fft_ALL{KK},Outlier_fft{KK}] = RemoveOutLier(E_fft_ALL{KK});
        [E_TFR_ALL{KK},Outlier_TFR{KK}] = RemoveOutLier(E_TFR_ALL{KK});
        [E_Up_ALL{KK},Outlier_Up{KK}] = RemoveOutLier(E_Up_ALL{KK});
        [E_Low_ALL{KK},Outlier_Low{KK}] = RemoveOutLier(E_Low_ALL{KK});
        [Width_ALL{KK},Outlier_Width{KK}] = RemoveOutLier(Width_ALL{KK});
        disp(char(10))
        disp (['Activity =',Activity{act}(1:2),'       --------      Posture = ',Posture_Location{PL}(1:2)])
        disp(['Number of Samples = ', num2str(L)])
        disp(['E-FFT =',num2str(mean(E_fft_ALL{KK})),' $\pm$ ',num2str(std(E_fft_ALL{KK})),'(',num2str(length(Outlier_fft{KK})),')']);
        disp(['E-TFR =',num2str(mean(E_TFR_ALL{KK})),' $\pm$ ',num2str(std(E_TFR_ALL{KK})),'(',num2str(length(Outlier_TFR{KK})),')']);
        disp(['E-Low =',num2str(mean(E_Low_ALL{KK})),' $\pm$ ',num2str(std(E_Low_ALL{KK})),'(',num2str(length(Outlier_Low{KK})),')']);
        disp(['E-Up =',num2str(mean(E_Up_ALL{KK})),' $\pm$ ',num2str(std(E_Up_ALL{KK})),'(',num2str(length(Outlier_Up{KK})),')']);
        disp(['Width =',num2str(mean(Width_ALL{KK})),' $\pm$ ',num2str(std(Width_ALL{KK})),'(',num2str(length(Outlier_Width{KK})),')']);
        KK = KK+1;
    end
end
