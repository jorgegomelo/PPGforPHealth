clear 
close all


load('ppg_features.mat');
load('ecgs.mat')
load('ppgs.mat')


RESTmeasures = struct;
STRESSmeasures = struct;

% Get ECG peaks
fs = 1000;

ecgs_peaks = findEcgPeaks(ecgs,fs);

window = 3*60*fs; %Janela de 3 minutos
RESTpacient = 1;
STRESSpacient = 1;

for q=1:length(ppg_features)
    ecg_mean = {};
    ecg_sdnn = {};
    ecg_rmssd = {};
    ecg_sdsd = {};
    ecg_nn50 = {};
    ecg_pnn50 = {};
    ppg_peak_mean = {};
    ppg_peak_sdnn = {};
    ppg_peak_rmssd = {};
    ppg_peak_sdsd = {};
    ppg_peak_nn50 = {};
    ppg_peak_pnn50 = {};
    ppg_80_mean = {};
    ppg_80_sdnn = {};
    ppg_80_rmssd = {};
    ppg_80_sdsd = {};
    ppg_80_nn50 = {};
    ppg_80_pnn50 = {};
    ppg_50_mean = {};
    ppg_50_sdnn = {};
    ppg_50_rmssd = {};
    ppg_50_sdsd = {};
    ppg_50_nn50 = {};
    ppg_50_pnn50 = {};
    ppg_20_mean = {};
    ppg_20_sdnn = {};
    ppg_20_rmssd = {};
    ppg_20_sdsd = {};
    ppg_20_nn50 = {};
    ppg_20_pnn50 = {};
    ppg_deriv_mean = {};
    ppg_deriv_sdnn = {};
    ppg_deriv_rmssd = {};
    ppg_deriv_sdsd = {};
    ppg_deriv_nn50 = {};
    ppg_deriv_pnn50 = {};
    ppg_onset_mean = {};
    ppg_onset_sdnn = {};
    ppg_onset_rmssd = {};
    ppg_onset_sdsd = {};
    ppg_onset_nn50 = {};
    ppg_onset_pnn50 = {};
    
    ecg_ralh = {};
    ppg_ralh = {};
    
    counter = 1;
    for i=1:window/10:length(ecgs{q})-(0.9*window)
        ecg_peaks = ecgs_peaks{q}(ecgs_peaks{q}>=i & ecgs_peaks{q}<(i+window));
        diff_ecg_peaks = diff(ecg_peaks);
        ppg_peaks = cell2mat(ppg_features(q).ppg_peaks);
        diff_ppg_peaks = diff(ppg_peaks(ppg_peaks>=i & ppg_peaks<(i+window)));
        ppg_80s = cell2mat(ppg_features(q).ppg_80s);
        diff_ppg_80s = diff(ppg_80s(ppg_80s>=i & ppg_80s<(i+window)));
        ppg_50s = cell2mat(ppg_features(q).ppg_50s);
        diff_ppg_50s = diff(ppg_50s(ppg_50s>=i & ppg_50s<(i+window)));
        ppg_20s = cell2mat(ppg_features(q).ppg_20s);
        diff_ppg_20s = diff(ppg_20s(ppg_20s>=i & ppg_20s<(i+window)));
        ppg_derivs = cell2mat(ppg_features(q).ppg_derivs);
        diff_ppg_derivs = diff(ppg_derivs(ppg_derivs>=i & ppg_derivs<(i+window)));
        ppg_onsets = cell2mat(ppg_features(q).ppg_onsets);
        diff_ppg_onsets = diff(ppg_onsets(ppg_onsets>=i & ppg_onsets<(i+window)));

        ecg_mean{counter} = mean(diff_ecg_peaks);
        ppg_peak_mean{counter} = mean(diff_ppg_peaks);
        ppg_80_mean{counter} = mean(diff_ppg_80s);
        ppg_50_mean{counter} = mean(diff_ppg_50s);
        ppg_20_mean{counter} = mean(diff_ppg_20s);
        ppg_deriv_mean{counter} = mean(diff_ppg_derivs);
        ppg_onset_mean{counter} = mean(diff_ppg_onsets);
        
        ecg_sdnn{counter} = std(diff_ecg_peaks);
        ppg_peak_sdnn{counter} = std(diff_ppg_peaks);
        ppg_80_sdnn{counter} = std(diff_ppg_80s);
        ppg_50_sdnn{counter} = std(diff_ppg_50s);
        ppg_20_sdnn{counter} = std(diff_ppg_20s);
        ppg_deriv_sdnn{counter} = std(diff_ppg_derivs);
        ppg_onset_sdnn{counter} = std(diff_ppg_onsets);

        ecg_rmssd{counter} = rms(diff(diff_ecg_peaks));
        ppg_peak_rmssd{counter} = rms(diff(diff_ppg_peaks));
        ppg_80_rmssd{counter} = rms(diff(diff_ppg_80s));
        ppg_50_rmssd{counter} = rms(diff(diff_ppg_50s));
        ppg_20_rmssd{counter} = rms(diff(diff_ppg_20s));
        ppg_deriv_rmssd{counter} = rms(diff(diff_ppg_derivs));
        ppg_onset_rmssd{counter} = rms(diff(diff_ppg_onsets));
        
        ecg_sdsd{counter} = std(diff(diff_ecg_peaks));
        ppg_peak_sdsd{counter} = std(diff(diff_ppg_peaks));
        ppg_80_sdsd{counter} = std(diff(diff_ppg_80s));
        ppg_50_sdsd{counter} = std(diff(diff_ppg_50s));
        ppg_20_sdsd{counter} = std(diff(diff_ppg_20s));
        ppg_deriv_sdsd{counter} = std(diff(diff_ppg_derivs));
        ppg_onset_sdsd{counter} = std(diff(diff_ppg_onsets));
        
        ecg_nn50{counter} = nnz(diff(diff_ecg_peaks)>50); %milissegundos
        ppg_peak_nn50{counter} = nnz(diff(diff_ppg_peaks)>50);
        ppg_80_nn50{counter} = nnz(diff(diff_ppg_80s)>50);
        ppg_50_nn50{counter} = nnz(diff(diff_ppg_50s)>50);
        ppg_20_nn50{counter} = nnz(diff(diff_ppg_20s)>50);
        ppg_deriv_nn50{counter} = nnz(diff(diff_ppg_derivs)>50);
        ppg_onset_nn50{counter} = nnz(diff(diff_ppg_onsets)>50);
        
        ecg_pnn50{counter} = nnz(diff(diff_ecg_peaks)>50)/length(diff_ecg_peaks);
        ppg_peak_pnn50{counter} = nnz(diff(diff_ppg_peaks)>50)/length(diff_ppg_peaks);
        ppg_80_pnn50{counter} = nnz(diff(diff_ppg_80s)>50)/length(diff_ppg_80s);
        ppg_50_pnn50{counter} = nnz(diff(diff_ppg_50s)>50)/length(diff_ppg_50s);
        ppg_20_pnn50{counter} = nnz(diff(diff_ppg_20s)>50)/length(diff_ppg_20s);
        ppg_deriv_pnn50{counter} = nnz(diff(diff_ppg_derivs)>50)/length(diff_ppg_derivs);
        ppg_onset_pnn50{counter} = nnz(diff(diff_ppg_onsets)>50)/length(diff_ppg_onsets);
%         % frequency features
% 
        if (i+window)>length(ecgs{q})
        	[ecg_high_freq] = pwelch(ecgs{q}(i:end),[],[],[0.15, 0.4]);
            [ecg_low_freq] = pwelch(ecgs{q}(i:end),[],[],[0.04, 0.15]);
            ecg_ralh{counter} = ecg_low_freq/ecg_high_freq;

            [ppg_high_freq] = pwelch(ppgs{q}(i:end),[],[],[0.15, 0.4]);
            [ppg_low_freq] = pwelch(ppgs{q}(i:end),[],[],[0.04, 0.15]);
            ppg_ralh{counter} = ppg_low_freq/ppg_high_freq;
        else
        	[ecg_high_freq] = pwelch(ecgs{q}(i:i+window),[],[],[0.15, 0.4]);
            [ecg_low_freq] = pwelch(ecgs{q}(i:i+window),[],[],[0.04, 0.15]);
            ecg_ralh{counter} = ecg_low_freq/ecg_high_freq;

            [ppg_high_freq] = pwelch(ppgs{q}(i:i+window),[],[],[0.15, 0.4]);
            [ppg_low_freq] = pwelch(ppgs{q}(i:i+window),[],[],[0.04, 0.15]);
            ppg_ralh{counter} = ppg_low_freq/ppg_high_freq;
        end

        
        counter = counter + 1;
    end
    if (rem(q,2)~=0)
            RESTmeasures(RESTpacient).ecg_mean = ecg_mean';
            RESTmeasures(RESTpacient).ppg_peak_mean = ppg_peak_mean';
            RESTmeasures(RESTpacient).ppg_80_mean = ppg_80_mean';
            RESTmeasures(RESTpacient).ppg_50_mean = ppg_50_mean';
            RESTmeasures(RESTpacient).ppg_20_mean = ppg_20_mean';
            RESTmeasures(RESTpacient).ppg_deriv_mean = ppg_deriv_mean';
            RESTmeasures(RESTpacient).ppg_onset_mean = ppg_onset_mean';
            
            RESTmeasures(RESTpacient).ecg_sdnn = ecg_sdnn';
            RESTmeasures(RESTpacient).ppg_peak_sdnn = ppg_peak_sdnn';
            RESTmeasures(RESTpacient).ppg_80_sdnn = ppg_80_sdnn';
            RESTmeasures(RESTpacient).ppg_50_sdnn = ppg_50_sdnn';
            RESTmeasures(RESTpacient).ppg_20_sdnn = ppg_20_sdnn';
            RESTmeasures(RESTpacient).ppg_deriv_sdnn = ppg_deriv_sdnn';
            RESTmeasures(RESTpacient).ppg_onset_sdnn = ppg_onset_sdnn';
            
            RESTmeasures(RESTpacient).ecg_rmssd = ecg_rmssd';
            RESTmeasures(RESTpacient).ppg_peak_rmssd = ppg_peak_rmssd';
            RESTmeasures(RESTpacient).ppg_80_rmssd = ppg_80_rmssd';
            RESTmeasures(RESTpacient).ppg_50_rmssd = ppg_50_rmssd';
            RESTmeasures(RESTpacient).ppg_20_rmssd = ppg_20_rmssd';
            RESTmeasures(RESTpacient).ppg_deriv_rmssd = ppg_deriv_rmssd';
            RESTmeasures(RESTpacient).ppg_onset_rmssd = ppg_onset_rmssd';
            
            RESTmeasures(RESTpacient).ecg_sdsd = ecg_sdsd';
            RESTmeasures(RESTpacient).ppg_peak_sdsd = ppg_peak_sdsd';
            RESTmeasures(RESTpacient).ppg_80_sdsd = ppg_80_sdsd';
            RESTmeasures(RESTpacient).ppg_50_sdsd = ppg_50_sdsd';
            RESTmeasures(RESTpacient).ppg_20_sdsd = ppg_20_sdsd';
            RESTmeasures(RESTpacient).ppg_deriv_sdsd = ppg_deriv_sdsd';
            RESTmeasures(RESTpacient).ppg_onset_sdsd = ppg_onset_sdsd';
            
            RESTmeasures(RESTpacient).ecg_nn50 = ecg_nn50';
            RESTmeasures(RESTpacient).ppg_peak_nn50 = ppg_peak_nn50';
            RESTmeasures(RESTpacient).ppg_80_nn50 = ppg_80_nn50';
            RESTmeasures(RESTpacient).ppg_50_nn50 = ppg_50_nn50';
            RESTmeasures(RESTpacient).ppg_20_nn50 = ppg_20_nn50';
            RESTmeasures(RESTpacient).ppg_deriv_nn50 = ppg_deriv_nn50';
            RESTmeasures(RESTpacient).ppg_onset_nn50 = ppg_onset_nn50';
            
            RESTmeasures(RESTpacient).ecg_pnn50 = ecg_pnn50';
            RESTmeasures(RESTpacient).ecg_pnn50 = ecg_pnn50';
            RESTmeasures(RESTpacient).ppg_peak_pnn50 = ppg_peak_pnn50';
            RESTmeasures(RESTpacient).ppg_80_pnn50 = ppg_80_pnn50';
            RESTmeasures(RESTpacient).ppg_50_pnn50 = ppg_50_pnn50';
            RESTmeasures(RESTpacient).ppg_20_pnn50 = ppg_20_pnn50';
            RESTmeasures(RESTpacient).ppg_deriv_pnn50 = ppg_deriv_pnn50';
            RESTmeasures(RESTpacient).ppg_onset_pnn50 = ppg_onset_pnn50';
            
            RESTmeasures(RESTpacient).ecg_ralh = ecg_ralh';
            RESTmeasures(RESTpacient).ppg_ralh = ppg_ralh';
            
            RESTpacient = RESTpacient + 1;
    else
            STRESSmeasures(STRESSpacient).ecg_mean = ecg_mean';
            STRESSmeasures(STRESSpacient).ppg_peak_mean = ppg_peak_mean';
            STRESSmeasures(STRESSpacient).ppg_80_mean = ppg_80_mean';
            STRESSmeasures(STRESSpacient).ppg_50_mean = ppg_50_mean';
            STRESSmeasures(STRESSpacient).ppg_20_mean = ppg_20_mean';
            STRESSmeasures(STRESSpacient).ppg_deriv_mean = ppg_deriv_mean';
            STRESSmeasures(STRESSpacient).ppg_onset_mean = ppg_onset_mean';
            
            STRESSmeasures(STRESSpacient).ecg_sdnn = ecg_sdnn';
            STRESSmeasures(STRESSpacient).ppg_peak_sdnn = ppg_peak_sdnn';
            STRESSmeasures(STRESSpacient).ppg_80_sdnn = ppg_80_sdnn';
            STRESSmeasures(STRESSpacient).ppg_50_sdnn = ppg_50_sdnn';
            STRESSmeasures(STRESSpacient).ppg_20_sdnn = ppg_20_sdnn';
            STRESSmeasures(STRESSpacient).ppg_deriv_sdnn = ppg_deriv_sdnn';
            STRESSmeasures(STRESSpacient).ppg_onset_sdnn = ppg_onset_sdnn';
            
            STRESSmeasures(STRESSpacient).ecg_rmssd = ecg_rmssd';
            STRESSmeasures(STRESSpacient).ppg_peak_rmssd = ppg_peak_rmssd';
            STRESSmeasures(STRESSpacient).ppg_80_rmssd = ppg_80_rmssd';
            STRESSmeasures(STRESSpacient).ppg_50_rmssd = ppg_50_rmssd';
            STRESSmeasures(STRESSpacient).ppg_20_rmssd = ppg_20_rmssd';
            STRESSmeasures(STRESSpacient).ppg_deriv_rmssd = ppg_deriv_rmssd';
            STRESSmeasures(STRESSpacient).ppg_onset_rmssd = ppg_onset_rmssd';
            
            STRESSmeasures(STRESSpacient).ecg_sdsd = ecg_sdsd';
            STRESSmeasures(STRESSpacient).ppg_peak_sdsd = ppg_peak_sdsd';
            STRESSmeasures(STRESSpacient).ppg_80_sdsd = ppg_80_sdsd';
            STRESSmeasures(STRESSpacient).ppg_50_sdsd = ppg_50_sdsd';
            STRESSmeasures(STRESSpacient).ppg_20_sdsd = ppg_20_sdsd';
            STRESSmeasures(STRESSpacient).ppg_deriv_sdsd = ppg_deriv_sdsd';
            STRESSmeasures(STRESSpacient).ppg_onset_sdsd = ppg_onset_sdsd';
            
            STRESSmeasures(STRESSpacient).ecg_nn50 = ecg_nn50';
            STRESSmeasures(STRESSpacient).ppg_peak_nn50 = ppg_peak_nn50';
            STRESSmeasures(STRESSpacient).ppg_80_nn50 = ppg_80_nn50';
            STRESSmeasures(STRESSpacient).ppg_50_nn50 = ppg_50_nn50';
            STRESSmeasures(STRESSpacient).ppg_20_nn50 = ppg_20_nn50';
            STRESSmeasures(STRESSpacient).ppg_deriv_nn50 = ppg_deriv_nn50';
            STRESSmeasures(STRESSpacient).ppg_onset_nn50 = ppg_onset_nn50';
            
            STRESSmeasures(STRESSpacient).ecg_pnn50 = ecg_pnn50';
            STRESSmeasures(STRESSpacient).ecg_pnn50 = ecg_pnn50';
            STRESSmeasures(STRESSpacient).ppg_peak_pnn50 = ppg_peak_pnn50';
            STRESSmeasures(STRESSpacient).ppg_80_pnn50 = ppg_80_pnn50';
            STRESSmeasures(STRESSpacient).ppg_50_pnn50 = ppg_50_pnn50';
            STRESSmeasures(STRESSpacient).ppg_20_pnn50 = ppg_20_pnn50';
            STRESSmeasures(STRESSpacient).ppg_deriv_pnn50 = ppg_deriv_pnn50';
            STRESSmeasures(STRESSpacient).ppg_onset_pnn50 = ppg_onset_pnn50';
            
            STRESSmeasures(STRESSpacient).ecg_ralh = ecg_ralh';
            STRESSmeasures(STRESSpacient).ppg_ralh = ppg_ralh';
            
            STRESSpacient = STRESSpacient + 1;
    end
end

save('STRESSmeasures.mat','STRESSmeasures');
save('RESTmeasures.mat','RESTmeasures');
clear