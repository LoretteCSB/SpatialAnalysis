%% code to plot dynamically linear combination of genes and play with
%%sliders to change value coefficients

% Future ameliorations
% ajouter R2
% ajouter push button pour choisir delamination ou proliferation
% Demander la formule a l'utilisateur/nb de parametres libre
% faire une image RGB pour superposer les cartes
% carte Boris sur une intensit? de noir


%% Parameters
close all
clear all
clc % clear command line to avoid funy behavior (sometimes command line keeps printing lines continuoulsy
mins=-1;maxs=1;

%% 1 gene or 1 phenotype = 1 matrix
%load('Patt_CellProp_Rescaled_2017_sept_06th_FormulaTitle_raw_data')

%
% code used to import data
% path to get data
mydir    ='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_FormulaTitle/';
dirmask  ='FormulaTitle_2017_sept_6th_2_mask_without BH1/';
dirdata  ='Patt_CellProp_Rescaled_2017_sept_06th_FormulaTitle_raw_data/';

%% 1 gene or 1 phenotype = 1 matrix
%load('Patt_CellProp_Rescaled_2017_sept_06th_FormulaTitle_raw_data')

%
% code used to import data
% path to get data
mydir    ='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_test/';
dirmask  ='test_2017_sept_6th_2_mask_without BH1/';
dirdata  ='Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';
%dirphysio='ImageGivenByBoris/';
% data transformed from grayscale to double (im2double) ==> implies a rescaling to 0-1
BH1   =importimage([mydir,dirdata],'BH1_12apf_6.tif');
bi    =importimage([mydir,dirdata],'bi_12apf_2.tif');
DIAP1 =importimage([mydir,dirdata],'DIAP1_12apf_3.tif');
esg   =importimage([mydir,dirdata],'esg_12apf_4.tif');
eyg   =importimage([mydir,dirdata],'eyg_12apf_4.tif');
h     =importimage([mydir,dirdata],'h_12apf_5.tif');
hth   =importimage([mydir,dirdata],'hth_12apf_4.tif');
sd    =importimage([mydir,dirdata],'sd_12apf_7.tif');
Sr    =importimage([mydir,dirdata],'Sr_12apf_2.tif');
wg    =importimage([mydir,dirdata],'wg_12apf_7.tif');


%% Dichotomize variables
%{
imagesc(hth.*bi)
title('hth*bi')
figure
imagesc(eyg.*bi)
title('eyg*bi')

plotyn=0;
BH1=plotdic(BH1,plotyn);
bi=plotdic(bi,plotyn);
DIAP1=plotdic(DIAP1,plotyn);
esg=plotdic(esg,plotyn);
eyg=plotdic(eyg,plotyn);
h=plotdic(h,plotyn);
hth=plotdic(hth,plotyn);
sd=plotdic(sd,plotyn);
Sr=plotdic(Sr,plotyn);
wg=plotdic(wg,plotyn);
figure
imagesc(hth.*bi)
title('hth*bi dicho')
figure
imagesc(eyg.*bi)
title('eyg*bi dicho')
%}
mask_notum=imread([mydir,dirmask,'Mask_Notum_2.tif']);
mask_notum=im2double(mask_notum);
% several issues with those maps - need Boris input:
% - difficulty image from Maria: cell border is not one color but a gradient
% - image from Boris is not at the right size / postionning
% I need :

% - gray scale image without macrocheats or text,
% - image 2987 by 3779 with notum well positonned
% - to be able to identify the background and the cell border easily
%delam  = importimage([mydir,dirphysio],'Delamination_BIGwt2_nDelMax=4_0015_NEG.png');
%prolif = importimage([mydir,dirphysio],'Proliferation_BIGwt2_nDivMax=7_0015_NEG.png');
delam =imread([mydir,dirdata,'z_Delamination_12h.tif']);
prolif =imread([mydir,dirdata,'z_Proliferation_12h.tif']);

%}
%%  RESIZE IMAGE
%
scaling=0.10;
BH1   =imresize(BH1,scaling,'bicubic');
bi    =imresize(bi,scaling,'bicubic');
DIAP1 =imresize(DIAP1,scaling,'bicubic');
esg   =imresize(esg,scaling,'bicubic');
eyg   =imresize(eyg,scaling,'bicubic');
h     =imresize(h,scaling,'bicubic');
hth   =imresize(hth,scaling,'bicubic');
sd    =imresize(sd,scaling,'bicubic');
Sr    =imresize(Sr,scaling,'bicubic');
wg    =imresize(wg,scaling,'bicubic');
mask_notum=imresize(mask_notum,scaling);
delam =imresize(delam,scaling);
prolif =imresize(prolif,scaling);
%}

%%


LargMask  =(mask_notum==0);%.*(wg==0).*(bi==0).*(hth==0).*(eyg==0).*(Sr==0);
IndBckgrnd=find(LargMask);

%% Define function : Y=a1*Gene1 + a2*Gene2...
% functions used in slider for update (callback)

%function for delamination
myfunc = @(a1,a2,a3,a4,a5,a6,a7) a1*esg+a2*sd+a3*bi.*h+a4*h.*esg
%text for figure title
frmt='%.2f';
FormulaTitle=@(a1,a2,a3,a4,a5,a6,a7) strcat(num2str(a1,frmt),'esg +  ',num2str(a2,frmt),'bi + ',num2str(a3,frmt),'bi*h  - ',num2str(a4,frmt),'esg*h')

%% PLOT
a1=0;a2=0;a3=0;a4=0;a5=0;a6=0;a7=0;
alpha1=1;alpha2=1;
[r,c]=size(BH1);aspectratioimage =c/r;scale = 0.91;

% initialize figure
l=560*1.8 ;
w=460;
%hFig=figure('Name','mon image','Visible', 'on','units','pixel','Position',[400 300 560*1.5 420])
hFig=figure('Name','mon image','Visible', 'on','units','pixel','Position',[300 300 l w])

%'normalized','Position',[0.3 0.4 0.7 0.5])
C = [0 0 0; 1 0 0];
hAxe=gca;
imagesc(cat(3,zeros(r,c),delam,zeros(r,c)))
set(hAxe,'xtick',[]); set(hAxe,'ytick',[]);
posFig=hFig.Position
posAx=hAxe.Position;


set(hAxe,'units','pixel','Position',[40 25 420*scale*aspectratioimage   420*scale])
hAxe.Units='normalized';

txtTitle = uicontrol(hFig,'style','text','Units','normalized','position', [0.01,0.96,0.6,0.05],'String','Delamination','Fontsize',14)
mytitle=  @(a1,a2,a3,a4,a5,a6,a7) uicontrol(hFig,'style','text','Units','normalized','position',[0.01,0.91,0.6,0.05],'String',FormulaTitle(a1,a2,a3,a4,a5,a6,a7),'Fontsize',14);

mytitle(0,0,0,0,0,0,0);

% opacity alpha1
slid_blend1 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.60,0.83,0.1,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.1 0.10],'callback' , 'alpha1= get(slid_blend1,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
set(slid_blend1,'Value',get(slid_blend1,'Max'))
txtalpha1m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.56,0.9,0.04,0.03],'String',num2str(mins));
txtalpha1M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.7,0.90,0.04,0.03],'String',num2str(maxs));
txtalpha1Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.62,0.94,0.08,0.033],'String','Blending Genes ');
% opacity alpha2
slid_blend2 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.80,0.83,0.1,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.1 0.10],'callback' , 'alpha2= get(slid_blend2,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
set(slid_blend2,'Value',get(slid_blend2,'Max'))
txtalpha2m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.76,0.90,0.04,0.03],'String',num2str(mins));
txtalpha2M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.90,0.90,0.04,0.03],'String',num2str(maxs));
txtalpha2Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.82,0.94,0.08,0.033],'String','Blending Delam');

% slider 1: esg
slid1 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.73,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a1= get(slid1,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt1m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.80,0.04,0.03],'String',num2str(mins));
txt1M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.80,0.04,0.03],'String',num2str(maxs));
txt1Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.84,0.04,0.033],'String','esg');
set(slid1,'Value',0)

% slider 2: sd
slid2 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.63,0.22,0.1],'Min' ,mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a2= get(slid2,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt2m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.7,0.04,0.04],'String',num2str(mins));
txt2M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.7,0.04,0.04],'String',num2str(maxs));
txt2Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.74,0.04,0.04],'String','sd');
set(slid2,'Value',0)

% slider 3: bi*h

slid3 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.53,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a3= get(slid3,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt3m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.6,0.04,0.04],'String',num2str(mins));
txt3M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.6,0.04,0.04],'String',num2str(maxs));
txt3Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.64,0.04,0.04],'String','bi*h');
set(slid3,'Value',0)

% slider 4: esg*h
slid4 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.43,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a4= get(slid4,''value'');ImResult=Modif5(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt4m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.5,0.04,0.04],'String',num2str(mins));
txt4M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.5,0.04,0.04],'String',num2str(maxs));
txt4Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.54,0.04,0.04],'String','esg*h');
set(slid4,'Value',0)


%}
function y = importimage(dirpath,file)
y =imread([dirpath,file]);
% im2double converts uint8 to double + normalize the values between 0 and 1
y=im2double(y);
y=y/max(y(:));
end

%
function gene=plotdic(gene,plotyn)
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

%}