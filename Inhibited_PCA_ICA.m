filelist=dir('*.tif');
for File=1:length(filelist)
    PCA_ICA_results=[];
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],210,1,'C:\Temp\CROP');    
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [],  0.8,size(mixedsig,1),[],1e-6,5000);
    PCA_ICA_results_cat(File).ica_sig=ica_sig;
    PCA_ICA_results_cat(File).ica_filters=ica_filters;
    [PCA_ICA_results_cat(File).ica_segments, PCA_ICA_results_cat(File).segmentlabel, PCA_ICA_results_cat(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 50, 1);
    [PCA_ICA_results_cat(File).ROIs,PCA_ICA_results_cat(File).segmentabel,PCA_ICA_results_cat(File).centroids]  = CellsortSegmentation(ica_filters, 1, 2, 20, 0);    
    PCA_ICA_results_cat(File).Cell_sig = CellsortApplyFilter(filelist(File).name, PCA_ICA_results_cat(File).ROIs,[],movm,0);
    PCA_ICA_results.ROIs=PCA_ICA_results_cat(File).ROIs;
    PCA_ICA_results.Cell_sig=PCA_ICA_results_cat(File).Cell_sig;
    save(strcat(filelist(File).name(1:length(filelist(File).name)-4),'_PCAICA.mat'),'PCA_ICA_results','-v7.3');
end
clearvars mixedsig ica_sig ica_filters ica_A numiter mixedsig mixedfilters CovEvals covtrace movm movtm

