clear
close all


dirdata='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/raw/zfromBoris_GeneRescaled/';
[mypic,map]=imread([dirdata,'esg.tif']);
mypic=double(mypic);

ming=min(min(mypic));maxg=max(max(mypic));
clims = [ming maxg];

figure
imagesc(mypic,clims)
set(gca,'xtick',[]);set(gca,'ytick',[])

% filtering
sig=10;
mypic_filtered=imgaussfilt(mypic,sig);

figure
figure
imagesc(mypic_filtered,clims)
set(gca,'xtick',[]);set(gca,'ytick',[])

load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')
figure
imagesc(data_gene(9).gene)
set(gca,'xtick',[]);set(gca,'ytick',[])

