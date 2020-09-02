%% load up the files and make the raster plots

File_list_Naomi=dir("*.mat");
temp=zeros(1,length(File_list_Naomi));
for i=1:length(File_list_Naomi)
    test=regexp(File_list_Naomi(i).name,"\d\.mat$");
    if test
        temp(i)=1;
    end    
end
File_list_Naomi(~temp)=[];

IdealNegIdx={};
for i=1:length(File_list_Naomi)
    load(File_list_Naomi(i).name, 'rand_idx');  
    load(File_list_Naomi(i).name, 'idealTraces');
    [~,ideal_neg,~]=intersect(find(idealTraces(1:end-1,1)>0),rand_idx);
    IdealNegIdx{i}=ideal_neg;
end

IdealSpikes_cat=[];
IdealSpikes={};
for i=1:length(File_list_Naomi)
    load(File_list_Naomi(i).name, 'idealTraces');
    load(File_list_Naomi(i).name, 'spikes');
    temp=spikes.somas;
    if i==1
        IdealSpikes_cat=temp(idealTraces(1:end-1,1)>0,:);        
    else
        IdealSpikes_cat=vertcat(temp(idealTraces(1:end-1,1)>0,:),IdealSpikes_cat);        
    end
    IdealSpikes{i}=temp(idealTraces(1:end-1,1)>0,:);    
end
clearvars i temp spikes idealTraces

bg_idxC=zeros(1,10);
for i=1:10
    if i==1
        bg_idxC(i)=size(IdealSpikes{i},1)+1;
    else
        bg_idxC(i)=size(IdealSpikes{i},1)+bg_idxC(i-1)+1;
    end
end

IdealSpikes_cat_bin=zeros(size(IdealSpikes_cat,1),size(IdealSpikes_cat,2)/20);
for i=1:size(IdealSpikes_cat,2)/20
    IdealSpikes_cat_bin(:,i)=sum(IdealSpikes_cat(:,(i*20)-19:i*20),2);
end

IdealSpikes_bin={};
for i=1:length(IdealSpikes)
    temp=IdealSpikes{i};
    temp_bin=zeros(size(IdealSpikes{i},1),size(IdealSpikes{i},2)/20);
    for k=1:size(IdealSpikes_cat,2)/20
        temp_bin(:,k)=sum(temp(:,(k*20)-19:k*20),2);
    end
    IdealSpikes_bin{i}=temp_bin;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(IdealSpikes_cat_bin,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\IdealSpikes'),'-dsvg','-r0');

imagesc(IdealSpikes_cat_bin>0);colormap gray
for i=bg_idxC
    rectangle('FaceColor','r','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\IdealSpikes_binary'),'-dsvg','-r0');
clearvars i Fighandle j k test temp spks Spikes bg_idxC

%% test deconvolution algorithms on the ideal trace to directly compare to the spike

CellSortSpikes={};
CellSortSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    CellSortSpikes{i}=pdist2(full(CellsortFindspikes(IdealResponses{1}(1:end-1,:),2,0.2,20,1))',IdealSpikes_bin{i},'correlation');
    CellSortSpikesResults(i).Neg_correl_max=max(1-CellSortSpikes{i}(IdealNegIdx{i},:),[],2);
    temp=diag(CellSortSpikes{i});
    CellSortSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CellSortSpikesResults(i).Pos_correl_max=max(1-CellSortSpikes{i}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CellSortSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CellSortSpikesResults(i).mean_correl=mean(diag(CellSortSpikes{i}));
end

CaImAnSpikes={};
CaImAnSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    FluorescentTraces=IdealResponses{i}(1:end-1,:);
    SpikeInfer=zeros(size(FluorescentTraces));
    for ij=1:size(FluorescentTraces,1)
        SpikeInfer(ij,:)=deconvolveCa(FluorescentTraces(ij,:)','ar1','foopsi')';
    end
    CaImAnSpikes{i,1}=SpikeInfer;
    CaImAnSpikes{i,2}=pdist2(SpikeInfer,IdealSpikes_bin{i},'correlation');
    CaImAnSpikesResults(i).Neg_correl_max=max(1-CaImAnSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(CaImAnSpikes{i,2});
    CaImAnSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CaImAnSpikesResults(i).Pos_correl_max=max(1-CaImAnSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CaImAnSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CaImAnSpikesResults(i).mean_correl=mean(diag(CaImAnSpikes{i,2}));
end

%% Suite2p is python based so need to run in python

for i=1:length(IdealSpikes_bin)
    FluorescentTraces=IdealResponses{i}(1:end-1,:);
    save(strcat('FluoTraces_',num2str(i),'.mat'),'FluorescentTraces');
end

Suite2pSpikes={};
Suite2pSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    %load(strcat('Suite2pDeconv_',num2str(i),'.mat'));
    load(strcat('Suite2pDeconv',num2str(i),'b.mat'));
    Suite2pSpikes{i,1}=Spikes2;
    Suite2pSpikes{i,2}=pdist2(Spikes2,IdealSpikes_bin{i},'correlation');
    Suite2pSpikesResults(i).Neg_correl_max=max(1-Suite2pSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(Suite2pSpikes{i,2});
    Suite2pSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    Suite2pSpikesResults(i).Pos_correl_max=max(1-Suite2pSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    Suite2pSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    Suite2pSpikesResults(i).mean_correl=mean(diag(Suite2pSpikes{i,2}));
end

%% Same with the DL approach (CASCADE)
%Had to go from z-scored to DF data, we used a global mean to generate the
%DF/F to avoid the moving window issue from Fig.1

load('C:\Data\Inhibited neurons\NAOMI\predictions_df_traces_all.mat')
spike_rates_all=spike_rates;
CascadeSpikes={};num_temp=1;
CascadeSpikesResults=struct();
for i=1:length(IdealSpikes_bin)
    spike_rates=spike_rates_all(num_temp:num_temp-1+size(IdealSpikes_bin{i},1),:);
    spike_rates(isnan(spike_rates))=0;
    CascadeSpikes{i,1}=spike_rates;
    CascadeSpikes{i,2}=pdist2(CascadeSpikes{i,1},IdealSpikes_bin{i},'correlation');
    CascadeSpikesResults(i).Neg_correl_max=nanmax(1-CascadeSpikes{i,2}(IdealNegIdx{i},:),[],2);
    temp=diag(CascadeSpikes{i,2});
    CascadeSpikesResults(i).Neg_correl=1-temp(IdealNegIdx{i});
    CascadeSpikesResults(i).Pos_correl_max=nanmax(1-CascadeSpikes{i,2}(setdiff(1:end,IdealNegIdx{i}),:),[],2);
    CascadeSpikesResults(i).Pos_correl=1-temp(setdiff(1:end,IdealNegIdx{i}));
    CascadeSpikesResults(i).mean_correl=nanmean(diag(CascadeSpikes{i,2}));
    num_temp=num_temp+size(IdealSpikes_bin{i},1);
end


%%


MeanCorrel=[];
for i=1:length(IdealSpikes_bin)
    MeanCorrel(i,3)=mean(CellSortSpikesResults(i).Pos_correl);
    MeanCorrel(i,1)=mean(CaImAnSpikesResults(i).Pos_correl);
    MeanCorrel(i,2)=mean(Suite2pSpikesResults(i).Pos_correl);
    MeanCorrel(i,4)=mean(CascadeSpikesResults(i).Pos_correl);
    MeanCorrel(i,8)=mean(CellSortSpikesResults(i).Neg_correl);
    MeanCorrel(i,6)=mean(CaImAnSpikesResults(i).Neg_correl);
    MeanCorrel(i,7)=mean(Suite2pSpikesResults(i).Neg_correl);
    MeanCorrel(i,9)=mean(CascadeSpikesResults(i).Neg_correl);
end


%% representative deconv traces

x=linspace(0,size(idealTraces,2)/5,size(idealTraces,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 600]);
plot(x,zscore(CaImAnSpikes{1}(1,:)),'color',[166 33 255]/255,'LineWidth',2);
hold on;plot(x,zscore(CaImAnSpikes{1}(2,:)),'color',[89 255 0]/255,'LineWidth',2);xlim([0 40]);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','CaImAnSpikesPosNeg'),'-dsvg','-r0');
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','CaImAnSpikesPosNeg'),'-depsc','-r0');

x=linspace(0,size(idealTraces,2)/5,size(idealTraces,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 600]);
plot(x,zscore(Suite2pSpikes{1}(1,:)),'color',[166 33 255]/255,'LineWidth',2);
hold on;plot(x,zscore(Suite2pSpikes{1}(2,:)),'color',[89 255 0]/255,'LineWidth',2);xlim([0 40]);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','Suite2pSpikesPosNeg'),'-dsvg','-r0');
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','Suite2pSpikesPosNeg'),'-depsc','-r0');

x=linspace(0,size(idealTraces,2)/5,size(idealTraces,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 600]);
plot(x,zscore(PCAICASpikes{1}(1,:)),'color',[166 33 255]/255,'LineWidth',2);
hold on;plot(x,zscore(PCAICASpikes{1}(2,:)),'color',[89 255 0]/255,'LineWidth',2);xlim([0 40]);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','PCAICASpikesPosNeg'),'-dsvg','-r0');
print(Fighandle,strcat('I:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','PCAICASpikesPosNeg'),'-depsc','-r0');

x=linspace(0,size(idealTraces,2)/5,size(idealTraces,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 600]);
plot(x,zscore(CascadeSpikes{1}(1,:)),'color',[166 33 255]/255,'LineWidth',2);
hold on;plot(x,zscore(CascadeSpikes{1}(2,:)),'color',[89 255 0]/255,'LineWidth',2);xlim([0 40]);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('R:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','CascadePosNeg'),'-dsvg','-r0');
print(Fighandle,strcat('R:\INHIB2020-Q2046\Manuscript\Draft\Figure 4\','CascadePosNeg'),'-depsc','-r0');


%% binarized spikes

CellSortSpikes_cat=[];
for i=1:length(IdealSpikes_bin)
    if i==1
        CellSortSpikes_cat=full(CellsortFindspikes(IdealResponses{1}(1:end-1,:),2,0.2,20,1))';
    else
        CellSortSpikes_cat=vertcat(CellSortSpikes_cat,full(CellsortFindspikes(IdealResponses{1}(1:end-1,:),2,0.2,20,1))');
    end    
end

CaImAnSpikes_cat=[];
for i=1:length(IdealSpikes_bin)
    if i==1
        CaImAnSpikes_cat=CaImAnSpikes{i,1};
    else
        CaImAnSpikes_cat=vertcat(CaImAnSpikes_cat,CaImAnSpikes{i,1});
    end    
end


Suite2pSpikes_cat=[];
for i=1:length(IdealSpikes_bin)
    if i==1
        Suite2pSpikes_cat=Suite2pSpikes{i,1};
    else
        Suite2pSpikes_cat=vertcat(Suite2pSpikes_cat,Suite2pSpikes{i,1});
    end    
end


CascadeSpikes_cat=[];
for i=1:length(IdealSpikes_bin)
    if i==1
        CascadeSpikes_cat=CascadeSpikes{i,1};
    else
        CascadeSpikes_cat=vertcat(CascadeSpikes_cat,CascadeSpikes{i,1});
    end    
end


bg_idxC=zeros(1,10);
for i=1:10
    if i==1
        bg_idxC(i)=size(IdealSpikes{i},1)+1;
    else
        bg_idxC(i)=size(IdealSpikes{i},1)+bg_idxC(i-1)+1;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(CaImAnSpikes_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAnIdealSpikes'),'-dsvg','-r0');

imagesc(CaImAnSpikes_cat>mean(CaImAnSpikes_cat(:)));colormap gray
for i=bg_idxC
    rectangle('FaceColor','r','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAnIdealSpikes_binary'),'-dsvg','-r0');
%clearvars i Fighandle j k test temp spks Spikes bg_idxC

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(Suite2pSpikes_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2pIdealSpikes'),'-dsvg','-r0');

imagesc(Suite2pSpikes_cat>mean(Suite2pSpikes_cat(:)));colormap gray
for i=bg_idxC
    rectangle('FaceColor','r','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2pIdealSpikes_binary'),'-dsvg','-r0');
%clearvars i Fighandle j k test temp spks Spikes bg_idxC

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(CellSortSpikes_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CellSortIdealSpikes'),'-dsvg','-r0');

imagesc(CellSortSpikes_cat>mean(CellSortSpikes_cat(:)));colormap gray
for i=bg_idxC
    rectangle('FaceColor','r','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CellSortIdealSpikes_binary'),'-dsvg','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(CascadeSpikes_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('R:\INHIB2020-Q2046\Manuscript\Draft\SuppFig\CascadeIdealSpikes'),'-dsvg','-r0');

imagesc(CascadeSpikes_cat>mean(CascadeSpikes_cat(:)));colormap gray
for i=bg_idxC
    rectangle('FaceColor','r','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('R:\INHIB2020-Q2046\Manuscript\Draft\SuppFig\CascadeIdealSpikes_binary'),'-dsvg','-r0');
clearvars i Fighandle j k test temp spks Spikes bg_idxC
