load('Inhibited_volume_100by50_10percent_poisson_1.mat','spikes','vol_out','spike_opts')
% spike_predef = spikes.somas;
% rand_idx = randperm(size(spike_predef,1),ceil(size(spike_predef,1)/(100/percent)));
% spike_predef(rand_idx,:)=ones(length(rand_idx),size(spike_predef,2))*7.6e-06;
% %spike_predef(rand_idx,1:500)=0;
% temp=poissrnd(1,length(rand_idx),size(spikes.somas,2));
% spike_predef(rand_idx,:)=spike_predef(rand_idx,:).*temp;
% for i=1:length(rand_idx)
%     idx_time_temp=randperm(size(spike_predef,2)-100,ceil(size(spike_predef,2)/400));
%     for time_idx=1:10
%         spike_predef(rand_idx(i),idx_time_temp(time_idx):idx_time_temp(time_idx)+randi([20 500]))=0;
%     end
% end
% spike_predef=spike_predef(:,1:size(spikes.somas,2));
% 
% [neur_act,spikes] = generateTimeTraces(spike_opts,spike_predef,vol_out.locs); 
% 

temp=poissrnd(1,size(spikes.somas,1),size(spikes.somas,2));
spike_predef = ones(size(spikes.somas))*7.6e-06;
firing_rate=logspace(-1,3,12);
idx_time_temp=linspace(500,size(spikes.somas,2)-500,10);
for i=1:length(firing_rate)
    temp=poissrnd(firing_rate(i)/100,size(spike_predef(1+200*(i-1):200+200*(i-1),:),1),size(spikes.somas,2));
    spike_predef(1+200*(i-1):200+200*(i-1),:)=spike_predef(1+200*(i-1):200+200*(i-1),:).*temp;
    for time_idx=1:10
        spike_predef(151+200*(i-1):200+200*(i-1),idx_time_temp(time_idx):idx_time_temp(time_idx)+300)=0;
    end    
end
spike_opts.burst_mean=1;
spike_opts.dendflag=0;

[neur_act,spikes] = generateTimeTraces(spike_opts,spike_predef,vol_out.locs); 

figure;imagesc(zscore(neur_act.soma(1:2400,:),1,2));

Neuronal_act=neur_act.soma(1:2400,:);
Spikes=spikes.somas(1:2400,:);

IdealSpikes_bin=zeros(size(Spikes,1),size(Spikes,2)/20);
for k=1:size(Spikes,2)/20
    IdealSpikes_bin(:,k)=sum(Spikes(:,(k*20)-19:k*20),2);
end

%% CaImAn

CaImAnSpikes={};
CaImAnSpikesResults=struct();
FluorescentTraces=Neuronal_act;
SpikeInfer=zeros(size(FluorescentTraces));
CleanFluo=zeros(size(FluorescentTraces));
for ij=1:size(FluorescentTraces,1)
    [temp2, temp, ~]=deconvolveCa(FluorescentTraces(ij,:)','ar1','foopsi');
    SpikeInfer(ij,:)=temp';
    CleanFluo(ij,:)=temp2';
end
CaImAnSpikes{1}=SpikeInfer;
CaImAnSpikes{2}=pdist2(SpikeInfer(:,2:end),IdealSpikes_bin(:,1:end-1),'correlation');
% CaImAnSpikes{3}=pdist2(conv2(1,gausswin(3),SpikeInfer),conv2(1,gausswin(3),IdealSpikes_bin),'correlation');
temp=diag(CaImAnSpikes{2});
% temp2=diag(CaImAnSpikes{3});

for i=1:length(firing_rate)
    CaImAnSpikesResults(i).Pos_Correl=1-temp(1+200*(i-1):150+200*(i-1));
    %     CaImAnSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CaImAnSpikesResults(1).mean_Poscorrel(i)=nanmean(CaImAnSpikesResults(i).Pos_Correl);
    %     CaImAnSpikesResults(1).mean_correl_Gauss(i)=nanmean(CaImAnSpikesResults(i).Correl_gauss);
    CaImAnSpikesResults(i).Neg_Correl=1-temp(151+200*(i-1):200+200*(i-1));
    %     CaImAnSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CaImAnSpikesResults(1).mean_Negcorrel(i)=nanmean(CaImAnSpikesResults(i).Neg_Correl);
end

figure;plot(CaImAnSpikesResults(1).mean_Poscorrel);hold on;plot(CaImAnSpikesResults(1).mean_Negcorrel);
% hold on;plot(CaImAnSpikesResults(1).mean_correl_Gauss);


roi_nb=6;
start=50;ending=150;
% figure;
% for i=1:10
%     subplot(2,5,i);
%     plot((Spikes(roi_nb+200*(i-1),start*20:ending*20))');xlim([0 200*20]);
% end

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1000, 1000],'visible','off');
xplot=3;yplot=4;
ha = tight_subplot(xplot,yplot,[.07 .07],[.07 .07],[.07 .07]);
for i=1:length(firing_rate)
    axes(ha(i));
    plot(zscore(FluorescentTraces(roi_nb+200*(i-1),start:ending))','k','LineWidth',2);hold on;
    plot(zscore(IdealSpikes_bin(roi_nb+200*(i-1),start:ending))','Color',[0.7 0 0],'LineWidth',2);hold on;
    title(strcat(num2str(firing_rate(i)),' Hz'));
    %plot(zscore(SpikeInfer(roi_nb+200*(i-1),start+1:ending+1))','Color',[0 0.5 0.5],'LineWidth',1);hold on;xlim([0 100]);
end
print(Fighandle,strcat('SimulatedFiringRates'),'-depsc','-r0');


%% CellSort

CellSortSpikes={};
CellSortSpikesResults=struct();

CellSortSpikes{1}=full(CellsortFindspikes(zscore(FluorescentTraces,1,2)',2,0.2,2,1));
CellSortSpikes{2}=pdist2(CellSortSpikes{1}(:,2:end),IdealSpikes_bin(:,1:end-1),'correlation');

Fighandle=figure;
set(Fighandle, 'Position', [10,10, 1000, 1000]);
xplot=3;yplot=4;
ha = tight_subplot(xplot,yplot,[.07 .07],[.07 .07],[.07 .07]);
for i=1:length(firing_rate)
    axes(ha(i));
    plot((FluorescentTraces(roi_nb+200*(i-1),start:ending))','k','LineWidth',2);hold on;
    plot(zscore(IdealSpikes_bin(roi_nb+200*(i-1),start:ending))','Color',[0.7 0 0],'LineWidth',2);hold on;
    title(strcat(num2str(firing_rate(i)),' Hz'));
%     plot(zscore(CellSortSpikes{1}(roi_nb+200*(i-1),start+1:ending+1))','Color',[0 0.5 0.5],'LineWidth',1);hold on;xlim([0 100]);
end
print(Fighandle,strcat('SimulatedFiringRates'),'-depsc','-r0');

temp=diag(CellSortSpikes{2});
% temp2=diag(CellSortSpikes{3});

for i=1:length(firing_rate)
    CellSortSpikesResults(i).Pos_Correl=1-temp(1+200*(i-1):150+200*(i-1));
    %     CellSortSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CellSortSpikesResults(1).mean_Poscorrel(i)=nanmean(CellSortSpikesResults(i).Pos_Correl);
    %     CellSortSpikesResults(1).mean_correl_Gauss(i)=nanmean(CellSortSpikesResults(i).Correl_gauss);
    CellSortSpikesResults(i).Neg_Correl=1-temp(151+200*(i-1):200+200*(i-1));
    %     CellSortSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CellSortSpikesResults(1).mean_Negcorrel(i)=nanmean(CellSortSpikesResults(i).Neg_Correl);
end

figure;plot(CellSortSpikesResults(1).mean_Poscorrel);hold on;plot(CellSortSpikesResults(1).mean_Negcorrel);


%% suite2p

save(strcat('FluoTraces.mat'),'FluorescentTraces');
load('Suite2pDeconv.mat')
Suite2pSpikes={};
Suite2pSpikesResults=struct();
Suite2pSpikes{1}=Spikes2;
Suite2pSpikes{2}=pdist2(Spikes2(:,2:end),IdealSpikes_bin(:,1:end-1),'correlation');

temp=diag(Suite2pSpikes{2});
% temp2=diag(Suite2pSpikes{3});

for i=1:length(firing_rate)
    Suite2pSpikesResults(i).Pos_Correl=1-temp(1+200*(i-1):150+200*(i-1));
    %     Suite2pSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    Suite2pSpikesResults(1).mean_Poscorrel(i)=nanmean(Suite2pSpikesResults(i).Pos_Correl);
    %     Suite2pSpikesResults(1).mean_correl_Gauss(i)=nanmean(Suite2pSpikesResults(i).Correl_gauss);
    Suite2pSpikesResults(i).Neg_Correl=1-temp(151+200*(i-1):200+200*(i-1));
    %     Suite2pSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    Suite2pSpikesResults(1).mean_Negcorrel(i)=nanmean(Suite2pSpikesResults(i).Neg_Correl);
end

figure;plot(Suite2pSpikesResults(1).mean_Poscorrel);hold on;plot(Suite2pSpikesResults(1).mean_Negcorrel);


%% MLspike

min_DF=min(FluorescentTraces,[],2);
min_DF(min_DF==0)=prctile(FluorescentTraces(min_DF==0,:),25,2);
min_DF(min_DF==0)=10;
DF=(FluorescentTraces-min_DF)./min_DF;
DF(max(DF,[],2)>5,:)=5*(DF(max(DF,[],2)>5,:)./max(DF(max(DF,[],2)>5,:),[],2)); %Some baselines ~0 so gives crazy DF
DF(~isfinite(DF))=0;
dF_traces=DF;save('DF.mat','dF_traces');
[P LL]=tps_mlspikes(DF',par);
P=P';

MLspikeSpikes={};
MLspikeSpikesResults=struct();
MLspikeSpikes{1}=P;
MLspikeSpikes{2}=pdist2(MLspikeSpikes{1}(:,2:end),IdealSpikes_bin(:,1:end-1),'correlation');

temp=diag(MLspikeSpikes{2});
% temp2=diag(MLspikeSpikes{3});

for i=1:length(firing_rate)
    MLspikeSpikesResults(i).Pos_Correl=1-temp(1+200*(i-1):150+200*(i-1));
    %     MLspikeSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    MLspikeSpikesResults(1).mean_Poscorrel(i)=nanmean(MLspikeSpikesResults(i).Pos_Correl);
    %     MLspikeSpikesResults(1).mean_correl_Gauss(i)=nanmean(MLspikeSpikesResults(i).Correl_gauss);
    MLspikeSpikesResults(i).Neg_Correl=1-temp(151+200*(i-1):200+200*(i-1));
    %     MLspikeSpikesResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    MLspikeSpikesResults(1).mean_Negcorrel(i)=nanmean(MLspikeSpikesResults(i).Neg_Correl);
end

figure;plot(MLspikeSpikesResults(1).mean_Poscorrel);hold on;plot(MLspikeSpikesResults(1).mean_Negcorrel);

%% CASCADE

load('predictions_DF.mat')

CascadeSpikes={};%num_temp=1;
CascadeSpikesResults=struct();
spike_rates(isnan(spike_rates))=0;
CascadeSpikes{1}=spike_rates(:,33:end-32);
CascadeSpikes{2}=pdist2(CascadeSpikes{1},IdealSpikes_bin(:,32:end-33),'correlation');
temp=diag(CascadeSpikes{2});

for i=1:length(firing_rate)
    CASCADEResults(i).Pos_Correl=1-temp(1+200*(i-1):150+200*(i-1));
    %     CASCADEResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CASCADEResults(1).mean_Poscorrel(i)=nanmean(CASCADEResults(i).Pos_Correl);
    %     CASCADEResults(1).mean_correl_Gauss(i)=nanmean(CASCADEResults(i).Correl_gauss);
    CASCADEResults(i).Neg_Correl=1-temp(151+200*(i-1):200+200*(i-1));
    %     CASCADEResults(i).Correl_gauss=1-temp2(1+200*(i-1):200+200*(i-1));
    CASCADEResults(1).mean_Negcorrel(i)=nanmean(CASCADEResults(i).Neg_Correl);
end

figure;plot(CASCADEResults(1).mean_Poscorrel);hold on;plot(CASCADEResults(1).mean_Negcorrel);


CascadeSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
CascadeSpikesResults(i).Pos_correl_max=nanmax(1-CascadeSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
CascadeSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
CascadeSpikesResults(i).mean_correl=nanmean(diag(CascadeSpikes{i,2}));