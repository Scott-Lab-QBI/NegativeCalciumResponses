%% Loading all the data Fig.2
File_list_NAOMI=dir("*.mat");
temp=zeros(1,length(File_list_NAOMI));
for i=1:length(File_list_NAOMI)
    test=regexp(File_list_NAOMI(i).name,"\d\.mat$");
    if test
        temp(i)=1;
    end
end
File_list_NAOMI(~temp)=[];

idx_freq=zeros(1,length(File_list_NAOMI));
for i=1:length(File_list_NAOMI)
    test=regexp(File_list_NAOMI(i).name,"freq(.+)\.mat$",'tokens');
    idx_freq(i)=str2num(test{1}{1});
end

%Load the ideal responses
IdealResponses={};counter=1;
IdealResponses_cat_freq={};
IdealResponses_cat_all=[];
for freq_nb=unique(idx_freq)
    IdealResponses_cat=[];
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        load(File_list_NAOMI(freq_id).name, 'idealTraces');
        load(File_list_NAOMI(freq_id).name, 'neur_act');
        if i==1
            IdealResponses_cat=idealTraces(idealTraces(:,1)>0,:);
        else
            IdealResponses_cat=vertcat(idealTraces(idealTraces(:,1)>0,:),IdealResponses_cat);
        end
        IdealResponses{i,counter}=idealTraces(idealTraces(:,1)>0,:);
        FullResponses{i}=neur_act.soma;
        i=i+1;
    end
    IdealResponses_cat_all=vertcat(IdealResponses_cat_all,IdealResponses_cat);
    IdealResponses_cat_freq{counter}=IdealResponses_cat;
    counter=counter+1;
end

IdealNegIdx={};counter=1;
for freq_nb=unique(idx_freq)
    IdealResponses_cat=[];
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        load(File_list_NAOMI(freq_id).name, 'rand_idx');
        load(File_list_NAOMI(freq_id).name, 'idealTraces');
        [~,ideal_neg,~]=intersect(find(idealTraces(:,1)>0),rand_idx);
        IdealNegIdx{i,counter}=ideal_neg;
        i=i+1;
    end
    counter=counter+1;
end

bg_idx=zeros(4,5);counter=1;
for freq_nb=unique(idx_freq)
    IdealResponses_cat=[];
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        if i==1
            bg_idx(counter,i)=size(IdealResponses{i,counter},1);
        else
            bg_idx(counter,i)=size(IdealResponses{i,counter},1)+bg_idx(counter,i-1);
        end
        i=i+1;
    end
    counter=counter+1;
end

PurpGreen = zeros(100,3);
PurpGreen(1:33,[1 3])=repmat(flip([0:1/32:1]),2,1)';
PurpGreen(33:end,2)=[0:1/67:1];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 1200]);
xplot=4;yplot=1;
ha=tight_subplot(yplot,xplot);
for freq_nb=1:length(unique(idx_freq))
    axes(ha(freq_nb));
    imagesc(zscore(IdealResponses_cat_freq{freq_nb},1,2),[-3 6]);colormap(PurpGreen)
    for i=bg_idx(freq_nb,:)
        rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
    end
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\IdealResponsesFreq'),'-dsvg','-r0');

%Load the ideal responses
CaImAnResponses={};counter=1;
CaImAnResponses_cat_freq={};
CaImAnResponses_cat_all=[];
for freq_nb=unique(idx_freq)
    CaImAnResponses_cat=[];
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        load(strrep(File_list_NAOMI(freq_id).name,'.mat','_output_analysis_matlab.mat'), 'DenoisedTraces');
        if i==1
            CaImAnResponses_cat=DenoisedTraces(max(DenoisedTraces,[],2)>0,:);
        else
            CaImAnResponses_cat=vertcat(DenoisedTraces(max(DenoisedTraces,[],2)>0,:),CaImAnResponses_cat);
        end
        CaImAnResponses_cat=vertcat(zeros(1,size(DenoisedTraces,2)),CaImAnResponses_cat);
        CaImAnResponses{i,counter}=DenoisedTraces(max(DenoisedTraces,[],2)>0,:);
        i=i+1;
    end
    CaImAnResponses_cat_all=vertcat(CaImAnResponses_cat_all,CaImAnResponses_cat);
    CaImAnResponses_cat_freq{counter}=CaImAnResponses_cat;
    counter=counter+1;
end

bg_idx_CaImAn=zeros(4,5);counter=1;
for freq_nb=unique(idx_freq)
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        if i==1
            bg_idx_CaImAn(counter,i)=size(CaImAnResponses{i,counter},1);
        else
            bg_idx_CaImAn(counter,i)=size(CaImAnResponses{i,counter},1)+bg_idx_CaImAn(counter,i-1);
        end
        i=i+1;
    end
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 1200]);
xplot=4;yplot=1;
ha=tight_subplot(yplot,xplot);
for freq_nb=1:length(unique(idx_freq))
    axes(ha(freq_nb));
    imagesc(zscore(CaImAnResponses_cat_freq{freq_nb},1,2),[-3 6]);colormap(PurpGreen)
    for i=bg_idx_CaImAn(freq_nb,:)
        rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
    end
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAnResponsesFreq'),'-dsvg','-r0');


%Load the Suite2p responses
Suite2pResponses={};counter=1;
Suite2pResponses_cat_freq={};
Suite2pResponses_cat_all=[];
for freq_nb=unique(idx_freq)
    Suite2pResponses_cat=[];
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        load(strcat('C:\Data\Inhibited neurons\NAOMI\FreqTest\',num2str(freq_id),'\suite2p\plane0\Fall.mat'),'F');
        if i==1
            Suite2pResponses_cat=F;
        else
            Suite2pResponses_cat=vertcat(F,Suite2pResponses_cat);
        end
        Suite2pResponses{i,counter}=F;
        Suite2pResponses_cat=vertcat(zeros(1,size(DenoisedTraces,2)),Suite2pResponses_cat);
        i=i+1;
    end
    Suite2pResponses_cat_all=vertcat(Suite2pResponses_cat_all,Suite2pResponses_cat);
    Suite2pResponses_cat_freq{counter}=Suite2pResponses_cat;
    counter=counter+1;
end

bg_idx_Suite2p=zeros(4,5);counter=1;
for freq_nb=unique(idx_freq)
    idx_temp=find(idx_freq==freq_nb);
    i=1;
    for freq_id=idx_temp
        if i==1
            bg_idx_Suite2p(counter,i)=size(Suite2pResponses{i,counter},1);
        else
            bg_idx_Suite2p(counter,i)=size(Suite2pResponses{i,counter},1)+bg_idx_Suite2p(counter,i-1);
        end
        i=i+1;
    end
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1800, 1200]);
xplot=4;yplot=1;
ha=tight_subplot(yplot,xplot);
for freq_nb=1:length(unique(idx_freq))
    axes(ha(freq_nb));
    imagesc(zscore(Suite2pResponses_cat_freq{freq_nb},1,2),[-3 6]);colormap(PurpGreen)
    for i=bg_idx_Suite2p(freq_nb,:)
        rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
    end
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2pResponsesFreq'),'-dsvg','-r0');

%Load data CellSort
File_list=dir('*PCAICA.mat');
PCAICAResponses_cat=[];
PCAICAResponses={};
for i=1:length(File_list)
    load(File_list(i).name, 'PCA_ICA_results');
    if i==1
        PCAICAResponses_cat=PCA_ICA_results.Cell_sig;
    else
        PCAICAResponses_cat=vertcat(PCA_ICA_results.Cell_sig,PCAICAResponses_cat);
    end
    PCAICAResponses_cat=vertcat(zeros(1,size(PCA_ICA_results.Cell_sig,2)),PCAICAResponses_cat);
    PCAICAResponses{i}=PCA_ICA_results.Cell_sig;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(PCAICAResponses_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\PCAICAResponsesNoInhib'),'-dsvg','-r0');

%% Correlate ROIs

% %Correlate images CaImAn and suite2p
%
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 450, 450]);
% ax = axes('position',[0 0 1 1]);
% set(ax,'xtick',[])
% set(ax,'ytick',[])
% axis(ax,'equal')
%
% Image_correlations_Suite2p_full=cell(1,length(File_list_NAOMI));
% Image_correlations_Suite2p_full_mean=nan(1,length(File_list_NAOMI));
% Image_correlations_CaImAn_full=cell(1,length(File_list_NAOMI));
% Image_correlations_CaImAn_full_mean=nan(1,length(File_list_NAOMI));
% for i=1:length(File_list_NAOMI)
%     load(File_list_NAOMI(i).name, 'idealTraces');
%     load(File_list_NAOMI(i).name, 'comps');
%     load(strcat('C:\Data\Inhibited neurons\NAOMI\NewDatasets\',num2str(i),'\suite2p\plane0\Fall.mat'),'stat');
%     load(strcat('C:\Data\Inhibited neurons\NAOMI\NewDatasets\',num2str(i),'\suite2p\plane0\Fall.mat'),'ops');
%     load(strrep(File_list_NAOMI(i).name,'.mat','_output_analysis_matlab.mat'), 'ROIs_dense');
%     load(strrep(File_list_NAOMI(i).name,'.mat','_output_analysis_matlab.mat'), 'idx_components');
%     comps=comps(:,:,(idealTraces(1:size(idealTraces,1)-1,1)>0));
%     ROIs_Suite2p=zeros(ops.Lx,ops.Ly,length(stat));
%     for ij=1:length(stat)
%         temp=[stat{ij}.ypix;stat{ij}.xpix];
%         for ik=1:length(temp)
%             ROIs_Suite2p(temp(1,ik)+1,temp(2,ik)+1,ij)=stat{ij}.lam(ik);
%         end
%     end
%     temp=nan(size(comps,3),size(ROIs_Suite2p,3));
%     for ij=1:size(comps,3)
%         for jk=1:size(ROIs_Suite2p,3)
%             temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs_Suite2p(:,:,jk)).*ops.meanImg);
%         end
%     end
%     Image_correlations_Suite2p_full{i}=temp;
%     Image_correlations_Suite2p_full_mean(i)=nanmean(max(temp,[],1));
%     Correlation_image=imread(strrep(File_list_NAOMI(i).name,'.mat','_mean.tiff'));
%     ROIs=reshape(ROIs_dense,[size(Correlation_image) size(ROIs_dense,2)]);
%     temp=nan(size(comps,3),size(ROIs,3));
%     for ij=1:size(comps,3)
%         for jk=1:size(ROIs,3)
%             temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs(:,:,jk)));
%         end
%     end
%     Image_correlations_CaImAn_full{i}=temp;
%     Image_correlations_CaImAn_full_mean(i)=nanmean(max(temp,[],1));
%     Filename=strsplit(File_list_NAOMI(i).name,'.mat');
%     Filename=Filename{1};
%     imagesc(max(comps,[],3),[0 100]);colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Ideal_',Filename,'NoInhib'),'-dsvg','-r0');
%     imagesc(max(ROIs,[],3),[0 0.2]);colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAn_',Filename,'NoInhib'),'-dsvg','-r0');
%     imagesc(nanmax(ROIs_Suite2p.*ops.meanImg,[],3));colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename,'NoInhib'),'-dsvg','-r0');
%     imagesc(nanmax(ROIs_Suite2p,[],3));colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename,'bNoInhib'),'-dsvg','-r0');
%     imagesc(nanmax(ROIs_Suite2p.*ops.meanImgE,[],3));colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename,'cNoInhib'),'-dsvg','-r0');
% end
% close all;
%
% % Correlate images CellSort
% Image_correlations_Cellsort_full=cell(1,length(File_list));
% Image_correlations_Cellsort_full_mean=nan(1,length(File_list));
% File_list=dir('*PCAICA.mat');
% for i=1:length(File_list)
%     load(File_list(i).name, 'PCA_ICA_results');
%     load(strcat(File_list_NAOMI(i).name), 'idealTraces');
%     load(strcat(File_list_NAOMI(i).name), 'comps');
%     comps=comps(:,:,(idealTraces(1:size(idealTraces,1)-1,1)>0));
%     ROIs=PCA_ICA_results.ROIs;
%     Filename=strsplit(File_list_NAOMI(i).name,'.mat');
%     Filename=Filename{1};
%     imagesc(squeeze(max(ROIs,[],1)));colormap gray
%     print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CellSort_',Filename,'NoInhib'),'-dsvg','-r0');
%     temp=nan(size(comps,3),size(ROIs,1));
%     for ij=1:size(comps,3)
%         for jk=1:size(ROIs,1)
%             temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs(jk,:,:)));
%         end
%     end
%     Image_correlations_Cellsort_full{i}=temp;
%     Image_correlations_Cellsort_full_mean(i)=mean(max(temp,[],1));
% end
% close all;


%% Correlation between ideal responses and algorithms

CaImAn_correl=cell(4,5);
Suite2p_correl=cell(4,5);
PCAICA_correl=cell(4,5);
for freq_nb=1:4
    for freq_id=1:5
        CaImAn_correl{freq_nb,freq_id}=pdist2(CaImAnResponses{freq_id,freq_nb},IdealResponses{freq_id,freq_nb}(1:size(IdealResponses{freq_id,freq_nb},1)-1,:),'correlation');
        Suite2p_correl{freq_nb,freq_id}=pdist2(Suite2pResponses{freq_id,freq_nb},IdealResponses{freq_id,freq_nb}(1:size(IdealResponses{freq_id,freq_nb},1)-1,:),'correlation');
        PCAICA_correl{freq_nb,freq_id}=pdist2(PCAICAResponses{freq_id,freq_nb},IdealResponses{freq_id,freq_nb}(1:size(IdealResponses{freq_id,freq_nb},1)-1,:),'correlation');
    end
end

%% Proportion of ideal responses identified

CaImanResults=struct();Threshold=0.5;
Prop_pos=zeros(4,5);
Prop_pos_unique=zeros(4,5);
Prop_neg=zeros(4,5);
Prop_neg_unique=zeros(4,5);
for freq_nb=1:4
    for freq_id=1:5       
        [max_corel idx_corel]=max(1-CaImAn_correl{freq_nb,freq_id}(:,IdealNegIdx{freq_id,freq_nb}),[],1);
        temp=max_corel>Threshold;
        CaImanResults(freq_nb,freq_id).Neg_correl=[max_corel; idx_corel];
        Prop_neg(freq_nb,freq_id)=sum(temp)/length(max_corel);
        Prop_neg_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
        [max_corel idx_corel]=max(1-CaImAn_correl{freq_nb,freq_id}(:,setdiff(1:end,IdealNegIdx{freq_id,freq_nb})),[],1);
        temp=max_corel>Threshold;
        CaImanResults(freq_nb,freq_id).Pos_correl=[max_corel; idx_corel];
        Prop_pos(freq_nb,freq_id)=sum(temp)/length(max_corel);
        Prop_pos_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
    end
end
CaImanResults(1).Prop_pos=Prop_pos;
CaImanResults(1).Prop_pos_unique=Prop_pos_unique;
CaImanResults(1).Prop_neg=Prop_neg;
CaImanResults(1).Prop_neg_unique=Prop_neg_unique;

PCAICAResults=struct();Threshold=0.5;
Prop_pos=zeros(4,5);
Prop_pos_unique=zeros(4,5);
Prop_neg=zeros(4,5);
Prop_neg_unique=zeros(4,5);
for freq_nb=1:4
    for freq_id=1:5
        [max_corel idx_corel]=max(1-PCAICA_correl{freq_nb,freq_id}(:,IdealNegIdx{freq_id,freq_nb}),[],1);
        temp=max_corel>Threshold;
        PCAICAResults(freq_nb,freq_id).Neg_correl=[max_corel; idx_corel];
        if size(temp,1)==0
            Prop_neg(freq_nb,freq_id)=0;
            Prop_neg_unique(freq_nb,freq_id)=0;
        else
            Prop_neg(freq_nb,freq_id)=sum(temp)/length(max_corel);
            Prop_neg_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
        end
        [max_corel idx_corel]=max(1-PCAICA_correl{freq_nb,freq_id}(:,setdiff(1:end,IdealNegIdx{freq_id,freq_nb})),[],1);
        temp=max_corel>Threshold;
        PCAICAResults(freq_nb,freq_id).Pos_correl=[max_corel; idx_corel];
        if size(temp,1)==0
            Prop_pos(freq_nb,freq_id)=0;
            Prop_pos_unique(freq_nb,freq_id)=0;
        else
            Prop_pos(freq_nb,freq_id)=sum(temp)/length(max_corel);
            Prop_pos_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
        end
    end
end
PCAICAResults(1).Prop_pos=Prop_pos;
PCAICAResults(1).Prop_pos_unique=Prop_pos_unique;
PCAICAResults(1).Prop_neg=Prop_neg;
PCAICAResults(1).Prop_neg_unique=Prop_neg_unique;

Suite2pResults=struct();Threshold=0.5;
Prop_pos=zeros(4,5);
Prop_pos_unique=zeros(4,5);
Prop_neg=zeros(4,5);
Prop_neg_unique=zeros(4,5);
for freq_nb=1:4
    for freq_id=1:5
        [max_corel idx_corel]=max(1-Suite2p_correl{freq_nb,freq_id}(:,IdealNegIdx{freq_id,freq_nb}),[],1);
        temp=max_corel>Threshold;
        Suite2pResults(freq_nb,freq_id).Neg_correl=[max_corel; idx_corel];
        if size(temp,1)==0
            Prop_neg(freq_nb,freq_id)=0;
            Prop_neg_unique(freq_nb,freq_id)=0;
        else
            Prop_neg(freq_nb,freq_id)=sum(temp)/length(max_corel);
            Prop_neg_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
        end
        [max_corel idx_corel]=max(1-Suite2p_correl{freq_nb,freq_id}(:,setdiff(1:end,IdealNegIdx{freq_id,freq_nb})),[],1);
        temp=max_corel>Threshold;
        Suite2pResults(freq_nb,freq_id).Pos_correl=[max_corel; idx_corel];
        if size(temp,1)==0
            Prop_pos(freq_nb,freq_id)=0;
            Prop_pos_unique(freq_nb,freq_id)=0;
        else
            Prop_pos(freq_nb,freq_id)=sum(temp)/length(max_corel);
            Prop_pos_unique(freq_nb,freq_id)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
        end
    end
end
Suite2pResults(1).Prop_pos=Prop_pos;
Suite2pResults(1).Prop_pos_unique=Prop_pos_unique;
Suite2pResults(1).Prop_neg=Prop_neg;
Suite2pResults(1).Prop_neg_unique=Prop_neg_unique;

clearvars i j k temp max_corel idx_corel ax bg_idxC comps DenoisedTraces F Filename ideal_neg
clearvars ij ik idx_components idx_temp iscell neur_act ops rand_idx ROIs_dense ROIs ROIs_Suite2p stat test

%Concatenate data for plotting in Prism
PrismTemp=[CaImanResults(1).Prop_pos_unique;Suite2pResults(1).Prop_pos_unique;
    PCAICAResults(1).Prop_pos_unique];

MeanIdealCorrel=[];
for i=1:length(IdealResponses)
    MeanIdealCorrel(i,1)=mean(CaImanResults(i).Pos_correl(1,:));
    MeanIdealCorrel(i,2)=mean(Suite2pResults(i).Pos_correl(1,:));
    MeanIdealCorrel(i,3)=mean(PCAICAResults(i).Pos_correl(1,:));
end