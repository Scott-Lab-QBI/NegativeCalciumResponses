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

%Load the ideal responses
IdealResponses_cat=[];
IdealResponses={};
for i=1:length(File_list_NAOMI)
    load(File_list_NAOMI(i).name, 'idealTraces');
    load(File_list_NAOMI(i).name, 'neur_act');
    if i==1
        IdealResponses_cat=idealTraces(idealTraces(:,1)>0,:);        
    else
        IdealResponses_cat=vertcat(idealTraces(idealTraces(:,1)>0,:),IdealResponses_cat);        
    end
    IdealResponses{i}=idealTraces(idealTraces(:,1)>0,:);
    FullResponses{i}=neur_act.soma;
end

IdealNegIdx={};
for i=1:length(File_list_NAOMI)
    load(File_list_NAOMI(i).name, 'rand_idx');  
    load(File_list_NAOMI(i).name, 'idealTraces');
    [~,ideal_neg,~]=intersect(find(idealTraces(:,1)>0),rand_idx);
    IdealNegIdx{i}=ideal_neg;
end

bg_idx=zeros(1,10);
for i=1:10
    if i==1
        bg_idx(i)=size(IdealResponses{i},1);
    else
        bg_idx(i)=size(IdealResponses{i},1)+bg_idx(i-1);
    end
end


PurpGreen = zeros(100,3);
PurpGreen(1:33,[1 3])=repmat(flip([0:1/32:1]),2,1)';
PurpGreen(33:end,2)=[0:1/67:1];

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(IdealResponses_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idx
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\IdealResponses'),'-dsvg','-r0');

Fighandle=figure;
colormap(PurpGreen);caxis([-3 6]);colorbar;
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\IdealResponses_colorbar'),'-dsvg','-r0');

%Load CaImAn responses

CaImAnResponses_cat=[];
CaImAnResponses={};
for i=1:length(File_list_NAOMI)
    load(strrep(File_list_NAOMI(i).name,'.mat','_output_analysis_matlab.mat'), 'DenoisedTraces');
    if i==1
        CaImAnResponses_cat=DenoisedTraces(max(DenoisedTraces,[],2)>0,:);
    else
        CaImAnResponses_cat=vertcat(DenoisedTraces(max(DenoisedTraces,[],2)>0,:),CaImAnResponses_cat);
    end
    CaImAnResponses_cat=vertcat(zeros(1,size(DenoisedTraces,2)),CaImAnResponses_cat);
    CaImAnResponses{i}=DenoisedTraces(max(DenoisedTraces,[],2)>0,:);
end

bg_idxC=zeros(1,10);
for i=1:10
    if i==1
        bg_idxC(i)=size(CaImAnResponses{i},1)+1;
    else
        bg_idxC(i)=size(CaImAnResponses{i},1)+bg_idxC(i-1)+1;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(CaImAnResponses_cat,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAnResponses'),'-dsvg','-r0');

%Load suite2p responses

Suite2p_traces=[];
for i=1:length(File_list_NAOMI)    
    load(strcat('C:\Data\Inhibited neurons\NAOMI\Suite2P\dir',num2str(i),'\suite2p\plane0\Fall.mat'),'F');        
    if i==1
        Suite2p_traces=F;
    else
        Suite2p_traces=vertcat(F,Suite2p_traces);
    end
    Suite2p_traces=vertcat(zeros(1,size(F,2)),Suite2p_traces);
    Suite2p_trace{i}=F;
end
close all;

bg_idxC=zeros(1,10);
for i=1:10
    if i==1
        bg_idxC(i)=size(Suite2p_trace{i},1)+1;
    else
        bg_idxC(i)=size(Suite2p_trace{i},1)+bg_idxC(i-1)+1;
    end
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 600, 1200]);
imagesc(zscore(Suite2p_traces,1,2),[-3 6]);colormap(PurpGreen)
for i=bg_idxC
    rectangle('FaceColor','w','Position',[0 i-0.5 1000 1],'EdgeColor','none');
end
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2pResponses'),'-dsvg','-r0');

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
print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\PCAICAResponses'),'-dsvg','-r0');

%% Correlate ROIs

%Correlate images CaImAn and suite2p

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 450, 450]);
ax = axes('position',[0 0 1 1]);
set(ax,'xtick',[])
set(ax,'ytick',[])
axis(ax,'equal')

Image_correlations_Suite2p_full=cell(1,length(File_list_NAOMI));
Image_correlations_Suite2p_full_mean=nan(1,length(File_list_NAOMI));
Image_correlations_CaImAn_full=cell(1,length(File_list_NAOMI));
Image_correlations_CaImAn_full_mean=nan(1,length(File_list_NAOMI));
for i=1:length(File_list_NAOMI)
    load(File_list_NAOMI(i).name, 'idealTraces');
    load(File_list_NAOMI(i).name, 'comps');
    load(strcat('C:\Data\Inhibited neurons\NAOMI\Suite2P\dir',num2str(i),'\suite2p\plane0\Fall.mat'),'stat');    
    load(strcat('C:\Data\Inhibited neurons\NAOMI\Suite2P\dir',num2str(i),'\suite2p\plane0\Fall.mat'),'ops');
    load(strrep(File_list_NAOMI(i).name,'.mat','_output_analysis_matlab.mat'), 'ROIs_dense');
    load(strrep(File_list_NAOMI(i).name,'.mat','_output_analysis_matlab.mat'), 'idx_components');
    comps=comps(:,:,(idealTraces(1:size(idealTraces,1)-1,1)>0));
    ROIs_Suite2p=zeros(ops.Lx,ops.Ly,length(stat));
    for ij=1:length(stat)
        temp=[stat{ij}.ypix;stat{ij}.xpix];
        for ik=1:length(temp)
            ROIs_Suite2p(temp(1,ik)+1,temp(2,ik)+1,ij)=stat{ij}.lam(ik);
        end
    end
    temp=nan(size(comps,3),size(ROIs_Suite2p,3));
    for ij=1:size(comps,3)
        for jk=1:size(ROIs_Suite2p,3)
            temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs_Suite2p(:,:,jk)).*ops.meanImg);  
        end
    end
    Image_correlations_Suite2p_full{i}=temp;
    Image_correlations_Suite2p_full_mean(i)=nanmean(max(temp,[],1));
    temp=nan(size(comps,3),size(ROIs_Suite2p,3));
    ROIs=reshape(ROIs_dense,[size(Correlation_image) size(ROIs_dense,2)]);
    for ij=1:size(comps,3)
        for jk=1:size(ROIs_Suite2p,3)
            temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs(:,:,jk)));  
        end
    end
    Image_correlations_CaImAn_full{i}=temp;
    Image_correlations_CaImAn_full_mean(i)=nanmean(max(temp,[],1));
    Filename=strsplit(File_list_NAOMI(i).name,'.mat');
    Filename=Filename{1};
    imagesc(max(comps,[],3),[0 100]);colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Ideal_',Filename),'-dsvg','-r0');
    imagesc(max(ROIs,[],3),[0 0.2]);colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CaImAn_',Filename),'-dsvg','-r0');        
    imagesc(nanmax(ROIs_Suite2p.*ops.meanImg,[],3));colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename),'-dsvg','-r0');    
    imagesc(nanmax(ROIs_Suite2p,[],3));colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename,'b'),'-dsvg','-r0');    
    imagesc(nanmax(ROIs_Suite2p.*ops.meanImgE,[],3));colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\Suite2p_img_',Filename,'c'),'-dsvg','-r0');    
end
close all;

% Correlate images CellSort
Image_correlations_Cellsort_full=cell(1,length(File_list));
Image_correlations_Cellsort_full_mean=nan(1,length(File_list));
File_list=dir('*PCAICA.mat');
for i=1:length(File_list)
    load(File_list(i).name, 'PCA_ICA_results');    
    load(strcat(File_list_Naomi(i).folder,'\',File_list_Naomi(i).name), 'idealTraces');
    load(strcat(File_list_Naomi(i).folder,'\',File_list_Naomi(i).name), 'comps');
    load(strrep(strrep(File_list(i).name,'Inhibited','C:\Data\Inhibited neurons\NAOMI\Inhibited'),'_PCAICA.mat','.mat'), 'comps');
    load(strrep(strrep(File_list(i).name,'Inhibited','C:\Data\Inhibited neurons\NAOMI\Inhibited'),'_PCAICA.mat','.mat'), 'idealTraces');
    comps=comps(:,:,(idealTraces(1:size(idealTraces,1)-1,1)>0));
    ROIs=PCA_ICA_results.ROIs;
    Filename=strsplit(File_list_NAOMI(i).name,'.mat');
    Filename=Filename{1};
    imagesc(squeeze(max(ROIs,[],1)));colormap gray
    print(Fighandle,strcat('C:\Data\Inhibited neurons\Figures\CellSort_',Filename),'-dsvg','-r0');    
    temp=nan(size(comps,3),size(ROIs,1));
    for ij=1:size(comps,3)
        for jk=1:size(ROIs,1)
            temp(ij,jk)=corr2(squeeze(comps(:,:,ij)),squeeze(ROIs(jk,:,:)));  
        end
    end
    Image_correlations_Cellsort_full{i}=temp;
    Image_correlations_Cellsort_full_mean(i)=mean(max(temp,[],1));
end
close all;


%% Correlation between ideal responses and algorithms

CaImAn_correl=cell(2,10);
Suite2p_correl=cell(2,10);
PCAICA_correl=cell(2,10);
for i=1:length(IdealResponses)    
    CaImAn_correl{1,i}=pdist2(CaImAnResponses{i},IdealResponses{i}(1:size(IdealResponses{i},1)-1,:),'correlation');    
    CaImAn_correl{2,i}=pdist2(CaImAnResponses{i},IdealResponses{i}(IdealNegIdx{i},:),'correlation');
    Suite2p_correl{1,i}=pdist2(Suite2p_trace{i},IdealResponses{i}(1:size(IdealResponses{i},1)-1,:),'correlation');    
    Suite2p_correl{2,i}=pdist2(Suite2p_trace{i},IdealResponses{i}(IdealNegIdx{i},:),'correlation');
    PCAICA_correl{1,i}=pdist2(PCAICAResponses{i},IdealResponses{i}(1:size(IdealResponses{i},1)-1,:),'correlation');    
    PCAICA_correl{2,i}=pdist2(PCAICAResponses{i},IdealResponses{i}(IdealNegIdx{i},:),'correlation');
end

%% Proportion of ideal responses identified

CaImanResults=struct();Threshold=0.5;
Prop_neg=zeros(1,length(CaImAn_correl));
Prop_pos=zeros(1,length(CaImAn_correl));
Prop_neg_unique=zeros(1,length(CaImAn_correl));
Prop_pos_unique=zeros(1,length(CaImAn_correl));
for i=1:length(CaImAn_correl)
    [max_corel idx_corel]=max(1-CaImAn_correl{1,i}(:,IdealNegIdx{i}),[],1);
    CaImanResults(i).Neg_correl=[max_corel; idx_corel];
    temp=max_corel>Threshold;    
    Prop_neg(i)=sum(temp)/length(max_corel);
    Prop_neg_unique(i)=length(unique(idx_corel(temp)))/length(max_corel);
    [max_corel idx_corel]=max(1-CaImAn_correl{1,i}(:,setdiff(1:end,IdealNegIdx{i})),[],1);
    temp=max_corel>Threshold;    
    CaImanResults(i).Pos_correl=[max_corel; idx_corel];
    Prop_pos(i)=sum(temp)/length(max_corel);
    Prop_pos_unique(i)=length(unique(categorical(idx_corel(temp))))/length(max_corel);
end
CaImanResults(1).Prop_pos=Prop_pos;
CaImanResults(1).Prop_neg_unique=Prop_neg_unique;
CaImanResults(1).Prop_neg=Prop_neg;
CaImanResults(1).Prop_pos_unique=Prop_pos_unique;

Suite2pResults=struct();
Prop_neg=zeros(1,length(Suite2p_correl));
Prop_pos=zeros(1,length(Suite2p_correl));
Prop_neg_unique=zeros(1,length(Suite2p_correl));
Prop_pos_unique=zeros(1,length(Suite2p_correl));
for i=1:length(Suite2p_correl)
    [max_corel idx_corel]=max(1-Suite2p_correl{1,i}(:,IdealNegIdx{i}),[],1);
    Suite2pResults(i).Neg_correl=[max_corel; idx_corel];
    temp=max_corel>Threshold;    
    Prop_neg(i)=sum(temp)/length(max_corel);
    Prop_neg_unique(i)=length(unique(idx_corel(temp)))/length(max_corel);
    [max_corel idx_corel]=max(1-Suite2p_correl{1,i}(:,setdiff(1:end,IdealNegIdx{i})),[],1);
    temp=max_corel>Threshold;    
    Suite2pResults(i).Pos_correl=[max_corel; idx_corel];
    Prop_pos(i)=sum(temp)/length(max_corel);
    Prop_pos_unique(i)=length(unique(idx_corel(temp)))/length(max_corel);
end
Suite2pResults(1).Prop_pos=Prop_pos;
Suite2pResults(1).Prop_neg_unique=Prop_neg_unique;
Suite2pResults(1).Prop_neg=Prop_neg;
Suite2pResults(1).Prop_pos_unique=Prop_pos_unique;

PCAICAResults=struct();
Prop_neg=zeros(1,length(PCAICA_correl));
Prop_pos=zeros(1,length(PCAICA_correl));
Prop_neg_unique=zeros(1,length(PCAICA_correl));
Prop_pos_unique=zeros(1,length(PCAICA_correl));
for i=1:length(PCAICA_correl)
    [max_corel idx_corel]=max(1-PCAICA_correl{1,i}(:,IdealNegIdx{i}),[],1);
    PCAICAResults(i).Neg_correl=[max_corel; idx_corel];
    temp=max_corel>Threshold;    
    Prop_neg(i)=sum(temp)/length(max_corel);
    Prop_neg_unique(i)=length(unique(idx_corel(temp)))/length(max_corel);
    [max_corel idx_corel]=max(1-PCAICA_correl{1,i}(:,setdiff(1:end,IdealNegIdx{i})),[],1);
    temp=max_corel>Threshold;    
    PCAICAResults(i).Pos_correl=[max_corel; idx_corel];
    Prop_pos(i)=sum(temp)/length(max_corel);
    Prop_pos_unique(i)=length(unique(idx_corel(temp)))/length(max_corel);
end
PCAICAResults(1).Prop_pos=Prop_pos;
PCAICAResults(1).Prop_neg_unique=Prop_neg_unique;
PCAICAResults(1).Prop_neg=Prop_neg;
PCAICAResults(1).Prop_pos_unique=Prop_pos_unique;
clearvars i j k temp max_corel idx_corel ax bg_idxC comps DenoisedTraces F Filename ideal_neg
clearvars ij ik idx_components idx_temp iscell neur_act ops rand_idx ROIs_dense ROIs ROIs_Suite2p stat test

%Concatenate data for plotting in Prism
PrismTemp=[[CaImanResults(1).Prop_pos_unique CaImanResults(1).Prop_neg_unique];[Suite2pResults(1).Prop_pos_unique Suite2pResults(1).Prop_neg_unique];
[PCAICAResults(1).Prop_pos_unique PCAICAResults(1).Prop_neg_unique]];

MeanIdealCorrel=[];
for i=1:length(IdealResponses)
    MeanIdealCorrel(i,1)=mean(CaImanResults(i).Pos_correl(1,:));
    MeanIdealCorrel(i,2)=mean(Suite2pResults(i).Pos_correl(1,:));
    MeanIdealCorrel(i,3)=mean(PCAICAResults(i).Pos_correl(1,:));
    MeanIdealCorrel(i,4)=mean(CaImanResults(i).Neg_correl(1,:));
    MeanIdealCorrel(i,5)=mean(Suite2pResults(i).Neg_correl(1,:));
    MeanIdealCorrel(i,6)=mean(PCAICAResults(i).Neg_correl(1,:));
end