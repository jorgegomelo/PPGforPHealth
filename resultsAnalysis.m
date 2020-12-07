clear 
close all

load('RESTmeasures.mat')
load('STRESSmeasures.mat')

 
%% Merge data
AllRestMat = [];
AllStressMat = [];
AllLabeledData = [];
for q=1:length(RESTmeasures)
    RESTtable = struct2table(RESTmeasures(q));
    RESTarray = table2array(RESTtable);
    RESTmat = cell2mat(RESTarray);
    STRESStable = struct2table(STRESSmeasures(q));
    STRESSarray = table2array(STRESStable);
    STRESSmat = cell2mat(STRESSarray);
    
    AllRestMat = [AllRestMat; RESTmat];
    AllStressMat = [AllStressMat; STRESSmat];
end
label_0 = zeros(size(AllRestMat,1),1);
label_1 = ones(size(AllStressMat,1),1);
AllLabeledData = [AllRestMat, label_0];
AllLabeledData = [AllLabeledData; [AllStressMat, label_1]];

%Remove nan lines

for i=1:size(AllLabeledData,1)
    if (i>size(AllLabeledData,1))
        continue;
    end
    TF = isnan(AllLabeledData(i,:));
    x = TF(TF==1);
    if (~isempty(x))
        AllLabeledData = [AllLabeledData(1:i-1,:); AllLabeledData(i+1:end,:)];
    end
end


%% Boxplot merged data

for i=1:7:size(AllRestMat,2)-6
    box_title = strcat('Results for', measures{round((i+6)/7)}, ' tool.');

    figure('visible','off');
    boxplot([AllRestMat(:,i:i+1),AllStressMat(:,i:i+1)],'Notch','on','Labels', {labels{1}, labels{2}, labels{8}, labels{9}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{2},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

    figure('visible','off');
    boxplot([AllRestMat(:,i),AllRestMat(:,i+2),AllStressMat(:,i),AllStressMat(:,i+2)],'Notch','on','Labels', {labels{1}, labels{3}, labels{8}, labels{10}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{3},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

    figure('visible','off');
    boxplot([AllRestMat(:,i),AllRestMat(:,i+3),AllStressMat(:,i),AllStressMat(:,i+3)],'Notch','on','Labels', {labels{1}, labels{4}, labels{8}, labels{11}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{4},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

    figure('visible','off');
    boxplot([AllRestMat(:,i),AllRestMat(:,i+4),AllStressMat(:,i),AllStressMat(:,i+4)],'Notch','on','Labels', {labels{1}, labels{5}, labels{8}, labels{12}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{5},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

    figure('visible','off');
    boxplot([AllRestMat(:,i),AllRestMat(:,i+5),AllStressMat(:,i),AllStressMat(:,i+5)],'Notch','on','Labels', {labels{1}, labels{6}, labels{8}, labels{13}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{6},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

    figure('visible','off');
    boxplot([AllRestMat(:,i),AllRestMat(:,i+6),AllStressMat(:,i),AllStressMat(:,i+6)],'Notch','on','Labels', {labels{1}, labels{7}, labels{8}, labels{14}})
    title(box_title)
    xlabel('Features')
    ylabel('Metric (ms)')
    cd BoxplotsMergedData
    saveas(gcf,strcat('REST_VS_STRESS_',labels{7},'_M_',measures{round((i+6)/7)},'.png'))
    cd ..

end
figure('visible','off');
boxplot([AllRestMat(:,43),AllRestMat(:,44),AllStressMat(:,43),AllStressMat(:,44)],'Notch','on','Labels', {'REST ECG RALH','REST PPG RALH', 'STRESS ECG RALH', 'STRESS PPG RALH'})
title(strcat('Results for RALH tool.'))
xlabel('Features')
cd BoxplotsMergedData
saveas(gcf,strcat('REST_VS_STRESS_M_RALH','.png'))
cd ..
close all

%% Correlation indexes and rmse

corr_ppg_ecg = {};
rmse_ppg_ecg = {};
p_value = {};
count = 1;
for i=1:7:size(AllRestMat,2)-6
    [R,P] = corrcoef(AllLabeledData(:,i:i+6));
    mseRef = mean(AllLabeledData(:,i));
    RMSE = [sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+1)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i))),...
        sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+2)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i))),...
        sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+3)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i))),...
        sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+4)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i))),...
        sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+5)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i))),...
        sqrt(mean(((AllLabeledData(:,i)-AllLabeledData(:,i+6)).^2)))/(max(AllLabeledData(:,i))-min(AllLabeledData(:,i)))];
    corr_ppg_ecg{count} = R;
    rmse_ppg_ecg{count} = RMSE;
    p_value{count} = P;
    count = count + 1;
end

[R,P] = corrcoef(AllLabeledData(:,43:44));
mseRef = mean(var(AllLabeledData(:,43),1));
RMSE = sqrt(mean(((AllLabeledData(:,43)-AllLabeledData(:,44)).^2)))/(max(AllLabeledData(:,43))-min(AllLabeledData(:,43)));

%% $Kolmogorov-Smirnov test

p_value_ks = {};
count = 1;
for i=1:7:size(AllLabeledData,2)-6
    [~,p_ecg] = kstest(AllLabeledData(:,i));
    [~,p_ppg_peak] = kstest(AllLabeledData(:,i+1));
    [~,p_ppg_80] = kstest(AllLabeledData(:,i+2));
    [~,p_ppg_50] = kstest(AllLabeledData(:,i+3));
    [~,p_ppg_20] = kstest(AllLabeledData(:,i+4));
    [~,p_ppg_deriv] = kstest(AllLabeledData(:,i+5));
    [~,p_ppg_onset] = kstest(AllLabeledData(:,i+6));
    p_value_ks{count,1} = p_ecg;
    p_value_ks{count,2} = p_ppg_peak;
    p_value_ks{count,3} = p_ppg_80;
    p_value_ks{count,4} = p_ppg_50;
    p_value_ks{count,5} = p_ppg_20;
    p_value_ks{count,6} = p_ppg_deriv;
    p_value_ks{count,7} = p_ppg_onset;
    count = count + 1;
end
[~,p_ecg_ralh] = kstest(AllLabeledData(:,43));
[~,p_ppg_ralh] = kstest(AllLabeledData(:,44));

%% Kruskal Wallis test

p_value = {};
count = 1;
for i=1:7:size(AllLabeledData,2)-6
    p_ecg = kruskalwallis(AllLabeledData(:,i),AllLabeledData(:,end),'off');
    p_ppg_peak = kruskalwallis(AllLabeledData(:,i+1),AllLabeledData(:,end),'off');
    p_ppg_80 = kruskalwallis(AllLabeledData(:,i+2),AllLabeledData(:,end),'off');
    p_ppg_50 = kruskalwallis(AllLabeledData(:,i+3),AllLabeledData(:,end),'off');
    p_ppg_20 = kruskalwallis(AllLabeledData(:,i+4),AllLabeledData(:,end),'off');
    p_ppg_deriv = kruskalwallis(AllLabeledData(:,i+5),AllLabeledData(:,end),'off');
    p_ppg_onset = kruskalwallis(AllLabeledData(:,i+6),AllLabeledData(:,end),'off');
    p_value{count,1} = p_ecg;
    p_value{count,2} = p_ppg_peak;
    p_value{count,3} = p_ppg_80;
    p_value{count,4} = p_ppg_50;
    p_value{count,5} = p_ppg_20;
    p_value{count,6} = p_ppg_deriv;
    p_value{count,7} = p_ppg_onset;
    count = count + 1;
end
    
p_ecg_ralh = kruskalwallis(AllLabeledData(:,43),AllLabeledData(:,end),'off');
p_ppg_ralh = kruskalwallis(AllLabeledData(:,44),AllLabeledData(:,end),'off');