close all
clear all

clear
close all


%% Ajouter pushbutton pour avoir R2 / pas R2
%% Demander la formule a l'utilisateur?
% nb de parametres libre
% ajouter la normalisation de la carte 
% loader les donnees 

%% Parameters


mydir=  '/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_test/';
dirmask='test_2017_sept_6th_2_mask_without BH1/';
dirdata='Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';

%% Import Gene Expression Image (image with gene expression and phenotype) and organize data
% ==> matri X : 1 column = 1 gene

%rawfiles
k = struct2table(dir([mydir,dirdata,'*.tif']));
filenames = k.name(k.bytes>0);

%split gene expression and phenotype (files starting with z)
ix =regexp(filenames,'z*');
ix=arrayfun(@(x) length(x{:}),ix);
zfilenames= filenames(ix==1);%phenotypic image
filenames= filenames(ix==0);%genetic expression image
clear k ix


%% Create matrix with 1 column = 1 gene
global BH1
global esg
% global BH1,esg
[BH1,map]=imread([mydir,dirdata,'BH1_12apf_6.tif']);
[esg,map]=imread([mydir,dirdata,'esg_12apf_4.tif']);
BH1=double(BH1);esg=double(esg);


global a1 a2 x
a1=1;a2=1;
x=1:10;
%plot(x,myfunc(a1,a2))
fig1= figure
%ylim([0 100])
set(gca,'xtick',[])
set(gca,'ytick',[])

slid1=uicontrol(fig1,'style','slider','position', ...
    [100,50,150,20],'Min' , -1 , 'Max' , 1 , ...
    'callback' , 'a1 = get(slid1 , ''value'' ) ;Modif1(a1)');
%slid2=uicontrol(fig1,'style','slider','position', ...
%    [100,50,150,20] ,'Min' , -1 , 'Max' , 1 , ...
%    'callback' , 'a2 = get(slid1 , ''value'' ) ;Modif2(a2)');




