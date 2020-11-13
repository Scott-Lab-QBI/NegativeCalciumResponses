load('FluoTraces_1.mat')

%DF with global min
min_DF=min(FluorescentTraces,[],2);
min_DF(min_DF==0)=prctile(FluorescentTraces(min_DF==0,:),25,2);
DF=(FluorescentTraces-min_DF)./min_DF;
DF(max(DF,[],2)>5,:)=5*(DF(max(DF,[],2)>5,:)./max(DF(max(DF,[],2)>5,:),[],2)); %Some baselines ~0 so gives crazy DF

par = tps_mlspikes('par'); %Set parameters based on GCaMP6f and estimates from the paper
par.a = 0.1;
par.tau = 2;
par.dt = 0.2;
par.pnonlin=[];
par.finetune.sigma=0.025;
par.display = 'none';
% temp=spk_demoGUI;
% par=temp.pest;
[P LL]=tps_mlspikes(dF_traces',par);

File_list=dir('FluoTraces*.mat');
File_list(11)=File_list(2);
for i=3:10
    File_list(i-1)=File_list(i);
end
File_list(10)=File_list(11);
File_list(11)=[];

P={};LL={};
for i=1:length(File_list)
    %load(strrep(File_list_NAOMI(i).name,'.mat','_DF.mat'));
    load(File_list(i).name,'FluorescentTraces');
    min_DF=min(FluorescentTraces,[],2);
    min_DF(min_DF==0)=prctile(FluorescentTraces(min_DF==0,:),25,2);
    min_DF(min_DF==0)=10;
    DF=(FluorescentTraces-min_DF)./min_DF;
    DF(max(DF,[],2)>5,:)=5*(DF(max(DF,[],2)>5,:)./max(DF(max(DF,[],2)>5,:),[],2)); %Some baselines ~0 so gives crazy DF
    DF(~isfinite(DF))=0;
    [P{i} LL{i}]=tps_mlspikes(DF',par);    
end
clearvars i FluorescentTraces min_DF DF
save('ML_spikes_results_datasetMix.mat','-v7.3');


temp=tps_mlspikes;