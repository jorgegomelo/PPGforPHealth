clear
close all

cd ICTST
mat = dir('*.txt'); 
file = struct('data',cell(length(mat),1));
info = {};


for q = 1:length(mat) 
    
    index = [];
    DI = [];
    channel_1 = [];
    channel_2 = [];
    channel_3 = [];
    
    x = fopen(mat(q).name,'r');

    file2019=textscan(x,'%f\t%f\t%f\t%f\t%f','HeaderLines',3,'CollectOutput',1);
    if (contains(mat(q).name,'2018'))
        channel_1 = file2019{1,1}(:,4);
        channel_2 = file2019{1,1}(:,5);
        channel_3 = file2019{1,1}(:,3);
    else
        channel_1 = file2019{1,1}(:,3);
        channel_2 = file2019{1,1}(:,4);
        channel_3 = file2019{1,1}(:,5);
    end
    index = file2019{1,1}(:,1);
    DI = file2019{1,1}(:,2);
    
        
    fclose(x);
    data = struct('index',index,'DI',DI,'channel_1',channel_1,'channel_2',channel_2,'channel_3',channel_3);
    file(q).data = data;
end

save('data.mat', 'file');