clear all
close all
%% Preprocess data

% path to get data
mydir    ='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_test/';
dirmask  ='test_2017_sept_6th_2_mask_without BH1/';
dirdata  ='Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';
%% 1 gene or 1 phenotype = 1 matrix

% code used to import data
% data transformed from grayscale to double (im2double) ==> implies a rescaling to 0-1
ara   =importimage([mydir,dirdata],'ara_12apf_5.tif');
BH1   =importimage([mydir,dirdata],'BH1_12apf_6.tif');
bi    =importimage([mydir,dirdata],'bi_12apf_2.tif');
caup  =importimage([mydir,dirdata],'caup_12apf_2.tif');
DIAP1 =importimage([mydir,dirdata],'DIAP1_12apf_3.tif');
esg   =importimage([mydir,dirdata],'esg_12apf_4.tif');
exd   =importimage([mydir,dirdata],'exd_12apf_6.tif');
eyg   =importimage([mydir,dirdata],'eyg_12apf_4.tif');
h     =importimage([mydir,dirdata],'h_12apf_5.tif');
hh    =importimage([mydir,dirdata],'hh_12apf_2.tif');
hth   =importimage([mydir,dirdata],'hth_12apf_4.tif');
lgs   =importimage([mydir,dirdata],'lgs_12apf_3.tif');
mirr  =importimage([mydir,dirdata],'mirr_12apf_8.tif');
Pax   =importimage([mydir,dirdata],'Pax_12apf_1.tif');
pct   =importimage([mydir,dirdata],'ptc_12apf_2.tif');
salm  =importimage([mydir,dirdata],'salm_12apf_3.tif');
sd    =importimage([mydir,dirdata],'sd_12apf_7.tif');
Ush   =importimage([mydir,dirdata],'Ush_12apf_2.tif');
Sr    =importimage([mydir,dirdata],'Sr_12apf_2.tif');
usp   =importimage([mydir,dirdata],'usp_12apf_3.tif');
wg    =importimage([mydir,dirdata],'wg_12apf_7.tif');


%MASK
mask_notum=imread([mydir,dirmask,'Mask_Notum_2.tif']);
mask_notum=im2double(mask_notum);

% PHYSIOLOGICAL MAP
% several issues with those maps - need Boris input:
% - difficulty image from Maria: cell border is not one color but a gradient
% - image from Boris is not at the right size / postionning
% I need :
% - gray scale image without macrocheats or text,
% - image 2987 by 3779 with notum well positonned
% - to be able to identify the background and the cell border easily
delam =imread([mydir,dirdata,'z_Delamination_12h.tif']);
prolif =imread([mydir,dirdata,'z_Proliferation_12h.tif']);
%prolif = rgb2gray(prolif);
save('Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data')

%}
%%  RESIZE IMAGE
%
clear 
load('Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data')

scaling=0.08;
ara   =imresize(ara,scaling,'bicubic');
BH1   =imresize(BH1,scaling,'bicubic');
bi    =imresize(bi,scaling,'bicubic');
caup  =imresize(caup,scaling,'bicubic');
DIAP1 =imresize(DIAP1,scaling,'bicubic');
esg   =imresize(esg,scaling,'bicubic');
exd   =imresize(exd,scaling,'bicubic');
eyg   =imresize(eyg,scaling,'bicubic');
h     =imresize(h,scaling,'bicubic');
hh    =imresize(hh,scaling,'bicubic');
hth   =imresize(hth,scaling,'bicubic');
lgs   =imresize(lgs,scaling,'bicubic');
mirr  =imresize(mirr,scaling,'bicubic');
Pax   =imresize(Pax,scaling,'bicubic');
pct   =imresize(pct,scaling,'bicubic');
salm  =imresize(salm,scaling,'bicubic');
sd    =imresize(sd,scaling,'bicubic');
Sr    =imresize(Sr,scaling,'bicubic');
Ush   =imresize(Ush,scaling,'bicubic');
usp   =imresize(usp,scaling,'bicubic');
wg    =imresize(wg,scaling,'bicubic');
mask_notum=imresize(mask_notum,scaling);
delam =imresize(delam,scaling);
prolif =imresize(prolif,scaling);
%}

save('Patt_CellProp_Rescaled_2017_sept_06th_test_resized_data')

%% gaussian filtering : 
clear 
close all
load('Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data')

sig = 10;

ara   =imgaussfilt(ara,sig);
imagesc(BH1)
BH1   =imgaussfilt(BH1,sig);
figure;imagesc(BH1)
figure;imagesc(bi)
bi    =imgaussfilt(bi,sig);
figure;imagesc(bi)
caup  =imgaussfilt(caup,sig);
DIAP1 =imgaussfilt(DIAP1,sig);
esg   =imgaussfilt(esg,sig);
eyg   =imgaussfilt(eyg,sig);
exd   =imgaussfilt(exd,sig);
h     =imgaussfilt(h,sig);
hh    =imgaussfilt(hh,sig);
hth   =imgaussfilt(hth,sig);
lgs   =imgaussfilt(lgs,sig);
mirr  =imgaussfilt(mirr,sig);
Pax   =imgaussfilt(Pax,sig);
pct   =imgaussfilt(pct,sig);
salm  =imgaussfilt(salm,sig);
sd    =imgaussfilt(sd,sig);
Sr    =imgaussfilt(Sr,sig);
Ush   =imgaussfilt(Ush,sig);
usp   =imgaussfilt(usp,sig);
wg    =imgaussfilt(wg,sig);

%delam =imgaussfilt(delam,sig);
%prolif=imgaussfilt(prolif,sig);
save('Patt_CellProp_Rescaled_2017_sept_06th_test_filtered_data')

%% gaussian filtering + dichotomize : 
clear 
close all
load('Patt_CellProp_Rescaled_2017_sept_06th_test_filtered_data')

edges =[0,0.3,0.6,1];
values=[0,(edges(2)+edges(3))/2,1];
close all

gene=h;
imagesc(gene)
test=discretize(gene,edges,values);
figure ;imagesc(test)

ara   =discretize(ara,edges,values);
BH1   =discretize(BH1,edges,values);
bi    =discretize(bi,edges,values);
caup  =discretize(caup,edges,values);
DIAP1 =discretize(DIAP1,edges,values);
esg   =discretize(esg,edges,values);
eyg   =discretize(eyg,edges,values);
exd   =discretize(exd,edges,values);
h     =discretize(h,edges,values);
hh    =discretize(hh,edges,values);
hth    =discretize(hth,edges,values);
lgs   =discretize(lgs,edges,values);
mirr  =discretize(mirr,edges,values);
Pax   =discretize(Pax,edges,values);
pct   =discretize(pct,edges,values);
salm  =discretize(salm,edges,values);
sd    =discretize(sd,edges,values);
Sr    =discretize(Sr,edges,values);
wg    =discretize(wg,edges,values);
Ush   =discretize(Ush,edges,values);
usp   =discretize(usp,edges,values);
clear edges values test


save('Patt_CellProp_Rescaled_2017_sept_06th_test_filtered_discretized')

function y = importimage(dirpath,file)
y =imread([dirpath,file]);
% im2double converts uint8 to double + normalize the values between 0 and 1
y=im2double(y);
y=y/max(y(:));
end
function gene=plotdicho(gene,plotyn,edges)
gene = discretize(gene,edges)
if plotyn
    figure('visible','on')
    subplot(1,2,1)
    imagesc(gene)
    gene(gene<0.25)=0;
    subplot(1,2,2)
    imagesc(gene)
    
else
    gene(gene<0.25)=0;
end
end
