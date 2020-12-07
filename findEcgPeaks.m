function [ecgs_peaks] = findEcgPeaks(ecgs,fs)

    W  = 25*fs;  %.. window  = 5 seconds
    ecgs_peaks = {};
    %% Normalize ecgs
    for i=1:length(ecgs)
        ecgs{i} = normalize(ecgs{i});
    end

    %% Signal Processing
    energy_ecgs = ecgs;
    for i=1:length(ecgs)
        % Low-pass filter
        order = 4;
        wc = 15;
        fc = wc / (0.5 * fs);
        [b, a]=butter(order, fc,'low');
        e1 = filter(b, a, ecgs{i});

        % High-pass filter
        wc = 5;
        fc = wc / (0.5 * fs);
        [b,a] = butter(order, fc,'high');
        e2 = filter(b, a, e1);

        % Differentiation + Potentiation (squared-root)
        e3 = diff(e2);

        e4 = e3.^2;

        % Moving Average
        timeWindow = 0.2;
        N = timeWindow*fs;
        b = (1/N)*ones(N,1);
        a = 1;
        energy_ecgs{i} = filter(b, a, e4);

    end

    %% Find Peaks

    peak_idx = {};
    for i=1:size(energy_ecgs,2)
        minPeakHeight = 0.8*mean(energy_ecgs{i}); %Changed because 0.7*mean was too low threashold
        pause1 = 0.3*fs; %Max 150 batidas/min
        [~,locs] = findpeaks(energy_ecgs{i},'MinPeakHeight',minPeakHeight,'MinPeakDistance',pause1);
        peak_idx{i} = locs;
    end

    %% Solve Delay

    back = 0.2*fs;
    pred_peak_idx = {};

    for i=1:size(peak_idx,2)
        Is = [];
        ecg = abs(ecgs{i});
        for j=2:size(peak_idx{i})
            [~,I] = max(ecg(peak_idx{i}(j)-back:peak_idx{i}(j)));
            I = I + peak_idx{i}(j)-back;
            Is = [Is; I];
        end
        ecgs_peaks{i} = Is;

    end

%     %% Plot
% 
%     for i=1:size(file)
%         plot(file(i).data.channel_3)
%         hold on
%         plot(file(i).data.ecg_peaks,file(i).data.channel_3(file(i).data.ecg_peaks),'o')
%         hold off
%         pause;
%     end
end