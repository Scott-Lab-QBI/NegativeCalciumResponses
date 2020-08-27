%Data is from Favre-Bulle, I. A., G. Vanwalleghem, M. A. Taylor, H. Rubinsztein-Dunlop and E. K. Scott (2018). 
%"Cellular-Resolution Imaging of Vestibular Processing across the Larval Zebrafish Brain." 
%Curr Biol 28(23): 3711-3722 e3713

%ZS2 are the z-scored responses, and the K-Means result
%H0 and W0 are the nnmf output for reproducibility
%Stimuli represent the multiple laser powers and their timings, convolved with a GCaMP spike
%idx_rsq are all the ZS2 traces that passed our 0.1 r2 threshold following linear regression with the stimuli times

colors=[0.14 1 0.14; 0.7 0.4 1;0 0.6 0.6]*256;
x = linspace(0.25,246/4,246);
ZS_AVG2=zeros(size(ZS2,1),246);
parfor idx_ZS=1:size(ZS2,1)
    start=30;
    AVG=[];
    for i=1:3
        AVG(i,:)=ZS2(idx_ZS,start:start+40);
        start=start+40;
    end    
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:6
        for i=1:3
            temp(i,:)=ZS2(idx_ZS,start:start+40);
            start=start+40;
        end        
        temp=mean(temp,1);
        temp=temp-min(temp);       
        AVG=[AVG temp];
    end
    ZS_AVG2(idx_ZS,:)=AVG;   
end

ZS2_rsq=ZS2(idx_rsq,:);

KmeansMeanClusters=zeros(length(GoodBetas_merge),size(ZS2_rsq,2));
for Beta=1:length(GoodBetas_merge)    
    idx_temp2=find(idxKmeans_final_goodmemberInBrain_merge(idx_rsq)==GoodBetas_merge(Beta));    
    KmeansMeanClusters(Beta,:)=mean(ZS2_rsq(idx_temp2,:),1);    
end

figure;plot(KmeansMeanClusters');
opt = statset('MaxIter',1000,'Display','final','UseParallel',1);
[W0,H0] = nnmf(ZS2_rsq,15,'replicates',5,'options',opt,'algorithm','mult');

NMFMeanClusters=zeros(size(H0,1),size(ZS2_rsq,2));
for Beta=1:size(H0,1)    
    temp=W0(:,Beta);
    idx_temp2=find(temp>(mean(temp)+3*std(temp)));    
    NMFMeanClusters(Beta,:)=mean(ZS2_rsq(idx_temp2,:),1);    
end

CorrelNMF_Kmeans=1-pdist2(KmeansMeanClusters,NMFMeanClusters,'correlation');
CorrelNMF_Kmeans_max={};
[CorrelNMF_Kmeans_max{1} CorrelNMF_Kmeans_max{2}]=max(CorrelNMF_Kmeans,[],2);

H0_AVG=zeros(size(NMFMeanClusters,1),246);
parfor idx_ZS=1:size(NMFMeanClusters,1)
    start=30;
    AVG=[];
    for i=1:3
        AVG(i,:)=NMFMeanClusters(idx_ZS,start:start+40);
        start=start+40;
    end    
    AVG=mean(AVG,1);
    AVG=AVG-min(AVG);
    j=1;
    for j=2:6
        for i=1:3
            temp(i,:)=NMFMeanClusters(idx_ZS,start:start+40);
            start=start+40;
        end        
        temp=mean(temp,1);
        temp=temp-min(temp);       
        AVG=[AVG temp];
    end
    H0_AVG(idx_ZS,:)=AVG;   
end

for Beta=1:length(GoodBetas_merge)
    Fighandle=figure;
    set(Fighandle, 'Position', [10, 10, 400, 350]);
    idx_temp2=find(idxKmeans_final_goodmemberInBrain_merge==GoodBetas_merge(Beta));    
    for k=0:5
        meanToPlot=mean(ZS_AVG2(idx_temp2,4+(k*41):40+(k*41)),1);
        meanToPlot2=H0_AVG(CorrelNMF_Kmeans_max{2}(Beta),4+(k*41):40+(k*41));        
        hold on;
        plot(x(4+(k*40):40+(k*40)),meanToPlot-mean(meanToPlot(1:4)),'color',colors(Beta,:)/256,'LineWidth',3);axis([0 60 -2 5]);
        hold on;plot(x(4+(k*40):40+(k*40)),meanToPlot2-mean(meanToPlot2(1:4)),'k','LineWidth',2);
        if (k<5)
            rectangle('EdgeColor','none','FaceColor',[0.25, 0.25, 0.25, 0.4-k*0.06],'Position',[x(9+(k*40)) -3 1 8]);
        end
    end    
    set(gca,'FontSize',14);set(gca,'fontname','arial')
    print(Fighandle,strcat('_WB_Cluster-NMF_',num2str(Beta)),'-dsvg','-r0');
    print(Fighandle,strcat('_WB_Cluster-NMF_',num2str(Beta)),'-depsc','-r0');    
end

PurpGreen = zeros(100,3);
PurpGreen(1:33,[1 3])=repmat(flip([0:1/32:1]),2,1)';
PurpGreen(33:end,2)=[0:1/67:1];

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 400]);
imagesc(NMFMeanClusters,[-3 6]);colormap(PurpGreen);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('_WB_Cluster-NMF'),'-dsvg','-r0');
print(Fighandle,strcat('_WB_Cluster-NMF'),'-depsc','-r0');

Fighandle=figure;
set(Fighandle, 'Position', [10, 10, 1000, 100]);
imagesc(KmeansMeanClusters([1 3 2],:),[-3 6]);colormap(PurpGreen);set(gca,'FontSize',14);set(gca,'fontname','arial')
print(Fighandle,strcat('_WB_Cluster-KM'),'-dsvg','-r0');
print(Fighandle,strcat('_WB_Cluster-KM'),'-depsc','-r0');