clear
close all

cd ICTST
load('data.mat','file');
cd ..

%% Control Section
% Consider only the data between start button and stop button pressing by
% Erica
ecgs = {};
ppgs = {};
for q = 1:length(file)
    control = file(q).data.channel_2;
    press_idx = find(control>40000);
    find_start = press_idx(find(press_idx<100000)); %Search for start button before 100 sec. of total experiment
    find_stop = press_idx(find(press_idx>100000)); %Search for stop button after 100 sec. of total experiment;
    if (isempty(find_start))
        ecgs{q} =file(q).data.channel_3(1:round(mean(find_stop)));
        ppgs{q} = file(q).data.channel_1(1:round(mean(find_stop)));
    elseif (isempty(find_stop))
        ecgs{q} =file(q).data.channel_3(round(mean(find_start)):end);
        ppgs{q} = file(q).data.channel_1(round(mean(find_start)):end);
    else
        ecgs{q} =file(q).data.channel_3(round(mean(find_start)):round(mean(find_stop)));
        ppgs{q} = file(q).data.channel_1(round(mean(find_start)):round(mean(find_stop)));
    end
end
save('ecgs.mat','ecgs')
save('ppgs.mat','ppgs')
clear control press_idx find_start find_stop q
%% Filtering
% High and low pass filters with different parameters for PPG and ECG
fs = 1000; %Sampling Frequency (Hz)
filteredPpg = {};
for q = 1:length(file)
 
    %High pass filter
    fc = 1; %cut-off frequency in Hz
    wc = fc/(0.5*fs);
    order = 6;
    [b,a] = butter(order,wc,'high'); %hih-pass filter design
    filteredPpg{q} = filtfilt(b,a,ppgs{q}); %zero-phase filter for ppg
    filteredEcg{q} = filtfilt(b,a,ecgs{q}); %zero-phase filter for ecg
    
    %Low pass filter
    fcPPG = 18; %cut-off for PPG
    fcECG = 30; %cut-ogg for ECG
    wcPPG = fcPPG/(0.5*fs);
    wcECG = fcECG/(0.5*fs);
    order = 8;
    [bPPG,aPPG] = butter(order,wcPPG,'low');
    [bECG,aECG] = butter(order,wcECG,'low');
    filteredPpg{q} = filtfilt(bPPG,aPPG,filteredPpg{q});
    filteredEcg{q} = filtfilt(bECG,aECG,filteredEcg{q});

end
clear q fc wc order b a fcPPG fcECG wcPPG wcECG aPPG bPPG aECG bECG
%% Segmentation 
% From the peaks extracted with findEcgPeaks.m, we define ECG beat segment from
% 1/3 of the mean distance between peaks before each peak until 2/3 of this
% distance after each peak. 
seg_ecgs = {};
seg_ppg = {};
ecgs_peaks = findEcgPeaks(ecgs,fs);

%NN ECG interval (s)
nn_ecg = {};
for q=1:length(file)
    ecg_peaks = ecgs_peaks{q};
    num_peaks = length(ecg_peaks);
    num_points = length(filteredEcg{q});
    mean_distance = num_points/num_peaks;
    go_back = round(mean_distance/3);
    go_forward = round(mean_distance*2/3);
    
    % Searching for local minima in PPG with Minimum Prominence of 400.
    TF = islocalmin(filteredPpg{q},'MinProminence',400);
    ppg_local_mins_idx = 1:length(filteredPpg{q});
    ppg_local_mins_idx = ppg_local_mins_idx(TF);
    
    
    for j=2:num_peaks-1
        if (ecg_peaks(j)-go_back>0)
            splitEcg = filteredEcg{q}((ecg_peaks(j)-go_back):(ecg_peaks(j)+go_forward));
            splitEcg_index = (ecg_peaks(j)-go_back):(ecg_peaks(j)+go_forward);
            next_peaks = ppg_local_mins_idx(find(ppg_local_mins_idx>ecg_peaks(j)));
            if length(next_peaks)>1
                % Avoid the last segment of ECG
                % PPG segment is the curve that starts after ECG peak
                splitPpg = filteredPpg{q}(next_peaks(1)-50:next_peaks(2)+50);
                seg_ecgs{q,j,1} = ((ecg_peaks(j)-go_back):(ecg_peaks(j)+go_forward));
                seg_ecgs{q,j,2} = splitEcg;
                seg_ppg{q,j,1} = next_peaks(1)-50:next_peaks(2)+50;
                seg_ppg{q,j,2} = splitPpg;
                nn_ecg{q,j} = (ecg_peaks(j+1)-ecg_peaks(j))/fs;
            end
            
        end
    end
end
clear q j ecg_peaks num_peaks num_points mean_distance go_back go_forward TF ppg_local_mins_idx splitEcg splitEcg_index next_peaks splitPpg 
%% PPG derivatives
PPG_derivatives = {};
for q=1:length(file)
    len = length(filteredPpg{q});
    PPG_derivatives{q,1} = (nd5p(filteredPpg{q},1,len));
    wnd = round((31.3/1000)*fs);
    b = ones(1,wnd)/wnd;
    for i=2:3 %3rd derivative
        PPG_derivatives{q,i} = (nd5p(PPG_derivatives{q,i-1,:},1,len));
        PPG_derivatives{q,i} = filtfilt(b,1,PPG_derivatives{q,i,:});
    end
end
clear len wnd b
%% PPG features
ppg_features = struct;
for q=1:length(file)
    ppg_peaks = {};
    ppg_80s = {};
    ppg_50s = {};
    ppg_20s = {};
    ppg_derivs = {};
    ppg_onsets = {};
    for j=1:size(seg_ppg,2)
        if (~isempty(seg_ppg{q,j,2}))
            TF = islocalmin(seg_ppg{q,j,2});
            t = 1:length(seg_ppg{q,j,2});
            t = t(TF);
            min_idx = t(1);
            minima = seg_ppg{q,j,2}(min_idx);
            [f_max,locs] = findpeaks(seg_ppg{q,j,2}(min_idx:end),'MinPeakProminence',400);
            i = locs(1)+min_idx-1;
            ppg_peak = seg_ppg{q,j,1}(i);
            abs_signal = seg_ppg{q,j,2}+abs(minima);
            [~,idx] = min(abs(abs_signal(1:i)-0.8*(f_max(1)+abs(minima))));
            ppg_80 = seg_ppg{q,j,1}(idx);
            [~,idx] = min(abs(abs_signal(1:i)-0.5*(f_max(1)+abs(minima))));
            ppg_50 = seg_ppg{q,j,1}(idx);
            [~,idx] = min(abs(abs_signal(min_idx:i)-0.2*(f_max(1)+abs(minima))));
            ppg_20 = seg_ppg{q,j,1}(idx+min_idx-1);
            t = seg_ppg{q,j,1}(1:i);
            [~,idx] = max(PPG_derivatives{q,1}(t));
            ppg_deriv = seg_ppg{q,j,1}(idx);
            [~,idx] = max(PPG_derivatives{q,3}(t(1):ppg_20));
            ppg_onset = seg_ppg{q,j,1}(idx);
            ppg_peaks{j} = ppg_peak;
            ppg_80s{j} = ppg_80;
            ppg_50s{j} = ppg_50;
            ppg_20s{j} = ppg_20;
            ppg_derivs{j} = ppg_deriv;
            ppg_onsets{j} = ppg_onset;            
        end
    end
    ppg_features(q).ppg_peaks = ppg_peaks;
    ppg_features(q).ppg_80s = ppg_80s;
    ppg_features(q).ppg_50s = ppg_50s;
    ppg_features(q).ppg_20s = ppg_20s;
    ppg_features(q).ppg_derivs = ppg_derivs;
    ppg_features(q).ppg_onsets = ppg_onsets;
end

save('ppg_features.mat','ppg_features')

% Plot segments of PPG with features
for q=1:length(file)
    for j=1:size(seg_ppg,2)
        if (~isempty(seg_ppg{q,j,2}))
            plot(seg_ppg{q,j,1},seg_ppg{q,j,2})
            hold on
            plot(ppg_features(q).ppg_peaks{j},seg_ppg{q,j,2}(ppg_features(q).ppg_peaks{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold on
            plot(ppg_features(q).ppg_80s{j},seg_ppg{q,j,2}(ppg_features(q).ppg_80s{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold on
            plot(ppg_features(q).ppg_50s{j},seg_ppg{q,j,2}(ppg_features(q).ppg_50s{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold on
            plot(ppg_features(q).ppg_derivs{j},seg_ppg{q,j,2}(ppg_features(q).ppg_derivs{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold on
            plot(ppg_features(q).ppg_20s{j},seg_ppg{q,j,2}(ppg_features(q).ppg_20s{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold on
            plot(ppg_features(q).ppg_onsets{j},seg_ppg{q,j,2}(ppg_features(q).ppg_onsets{j}-(seg_ppg{q,j,1}(1)-1)),'o')
            hold off 
            pause;
        end
    end
end
