%% code to plot dynamically linear combination of genes and play with
%%sliders to change value coefficients

% Future am?liorations
% ajouter R2
% ajouter push button pour choisir delamination ou proliferation
% Demander la formule a l'utilisateur/nb de parametres libre

%% Parameters
close all
clear all

% specify if delamination or proliferation map
Delamination01=1;

%min max slider
mins=0;maxs=1;

%% 1 gene or 1 phenotype = 1 matrix
load('Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data')

%{
% code used to import data
% path to get data
mydir=  '/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_test/';
dirmask='test_2017_sept_6th_2_mask_without BH1/';
dirdata='Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';
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

mask_notum=imread([mydir,dirmask,'Mask_Notum_2.tif']);
masknotum=im2double(masknotum)
% several issues with those map - need Boris input
%delam  = importimage([mydir,dirphysio],'Delamination_BIGwt2_nDelMax=4_0015_NEG.png',0);
%prolif = importimage([mydir,dirphysio],'Proliferation_BIGwt2_nDivMax=7_0015_NEG.png',0);

%}
%% Define function : Y=a1*Gene1 + a2*Gene2...
% functions used in slider for update (callback)

if Delamination01==0
    %proliferation
    myfunc = @(a1,a2,a3,a4,a5,a6,a7) a1*wg-a2*bi-a3*hth-a4*eyg-a5*Sr+a6*eyg.*bi+a7*hth.*bi
    mytitle=  @(a1,a2,a3,a4,a5,a6,a7) title(strcat('a1= ',num2str(a1),'  a2= ',num2str(a2),'  a3=',num2str(a3),'  a4=',num2str(a4),'  a5=',num2str(a5),'  a6=',num2str(a6),'  a7=',num2str(a7)),'Fontsize',14);
    IndBckgrnd=find((mask_notum==0).*(wg==0).*(bi==0).*(hth==0).*(eyg==0).*(Sr==0));
else
    %delamination
    myfunc = @(a1,a2,a3,a4,a5,a6,a7) a1*esg-a2*sd+a3*bi.*h+a4*esg.*h
    mytitle=  @(a1,a2,a3,a4,a5,a6,a7) title(strcat('a1= ',num2str(a1),'  a2= ',num2str(a2),'  a3=',num2str(a3),'  a4=',num2str(a4)))
    IndBckgrnd=find((mask_notum==0).*(esg==0).*(sd==0).*(bi==0).*(h==0));
end
%% PLOT

% initialize figure
hFig=figure('Name','mon image','units','normalized','Position',[0.3 0.4 0.4 0.5])
C = [0 0 0; 1 0 0];
hAxe=gca;
set(hAxe,'xtick',[]); set(hAxe,'ytick',[]);
posFig=hFig.Position;
posAxe=hAxe.Position;
aspectratioimage = posAxe(4)/posAxe(3);scale = 0.65;
set(hAxe,'units','normalized','Position',[0.0500   , 0.100  ,  scale,    scale*aspectratioimage])

% figure title
if Delamination01==0
txtTitle = uicontrol(hFig,'style','text','Units','normalized','position', [0.1,0.9,0.5,0.1],'String','Proliferation : a1*wg - a2*bi - a3*hth - a4*eyg - a5*Sr + a6*eyg.*bi + a7*hth.*bi','Fontsize',14)
else
txtTitle = uicontrol(hFig,'style','text','Units','normalized','position', [0.1,0.9,0.5,0.1],'String','Delamination : a1*wg - a2*bi - a3*hth - a4*eyg - a5*Sr + a6*eyg.*bi + a7*hth.*bi','Fontsize',14)
end

a1=0;a2=0;a3=0;a4=0;a5=0;a6=0;a7=0;
% slider 1
slid1 = uicontrol(hFig,'style','slider','Units','normalized','position', [0.75,0.2,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a1= get(slid1,''value'');ImResult=Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt1m = uicontrol(hFig,'style','text','Units','normalized','position', [0.71,0.26,0.04,0.04],'String',num2str(mins))
txt1M = uicontrol(hFig,'style','text','Units','normalized','position', [0.97,0.26,0.04,0.04],'String',num2str(maxs))
txt1Title = uicontrol(hFig,'style','text','Units','normalized','position', [0.85,0.30,0.04,0.04],'String','a1')

% slider 2
slid2 = uicontrol(hFig,'style','slider','Units','normalized','position', [0.75,0.4,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a2= get(slid2,''value'');ImResult=Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt2m = uicontrol(hFig,'style','text','Units','normalized','position', [0.71,0.46,0.04,0.04],'String',num2str(mins))
txt2M = uicontrol(hFig,'style','text','Units','normalized','position', [0.97,0.46,0.04,0.04],'String',num2str(maxs))
txt2Title = uicontrol(hFig,'style','text','Units','normalized','position', [0.85,0.50,0.04,0.04],'String','a2')

% slider 3
slid3 = uicontrol(hFig,'style','slider','Units','normalized','position', [0.75,0.6,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a3= get(slid3,''value'');ImResult=Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt3m = uicontrol(hFig,'style','text','Units','normalized','position', [0.71,0.66,0.04,0.04],'String',num2str(mins))
txt3M = uicontrol(hFig,'style','text','Units','normalized','position', [0.97,0.66,0.04,0.04],'String',num2str(maxs))
txt3Title = uicontrol(hFig,'style','text','Units','normalized','position', [0.85,0.70,0.04,0.04],'String','a3')

% slider 4
slid4 = uicontrol(hFig,'style','slider','Units','normalized','position', [0.75,0.8,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a4= get(slid4,''value'');ImResult=Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');
txt4m = uicontrol(hFig,'style','text','Units','normalized','position', [0.71,0.86,0.04,0.04],'String',num2str(mins))
txt4M = uicontrol(hFig,'style','text','Units','normalized','position', [0.97,0.86,0.04,0.04],'String',num2str(maxs))
txt4Title = uicontrol(hFig,'style','text','Units','normalized','position', [0.85,0.90,0.04,0.04],'String','a4')

%slider 5


function y = importimage(dirpath,file)
y =imread([dirpath,file]);
% im2double converts uint8 to double + normalize the values between 0 and 1
y=im2double(y);
end

function fun=Modif2(fun,IndBckgrnd)
% set background to NaN
fun(IndBckgrnd)= NaN;
%normalize to [0,1]
fun = (fun -min(fun(:)))./(max(fun(:))-min(fun(:)));
imagesc(fun)
set(gca,'xtick',[])
set(gca,'ytick',[])

end

%{
% doesn't work because it update slid and not slide1, slide2... would problably need
to get adress of each handl
function [slid]=createSlider(hFig,pos,mins,maxs,IndBckgrnd,myfunc,mytitle,para)
global a1
global a2
global a3
global a4
global a5
global a6
global a7
%global a1, a2, a3, a4, a5, a6, a7
slid=uicontrol(hFig,'style','slider','Units','normalized','position', pos,'Min' , mins , 'Max' , maxs , ...
    'callback' , 'assignin(''base'',para, get(slid,''value''));Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');

%slid=uicontrol(hFig,'style','slider','Units','normalized','position', pos,'Min' , mins , 'Max' , maxs , ...
%    'callback' , 'assignin(''base'',para, get(slid,''value''));Modif2(myfunc(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd);mytitle(a1,a2,a3,a4,a5,a6,a7)');
%txtm=uicontrol(hFig,'style','text','Units','normalized','position', [pos(1)-0.4,pos(2)-0.05,0.05,0.05],'String',num2str(mins))
%txtM=uicontrol(hFig,'style','text','Units','normalized','position', [0.70,0.2,0.1,0.1],'String',num2str(maxs))
%txtTitle=uicontrol(hFig,'style','text','Units','normalized','position', [0.70,0.2,0.1,0.1],'String',num2str(txtTitle))

end
%}