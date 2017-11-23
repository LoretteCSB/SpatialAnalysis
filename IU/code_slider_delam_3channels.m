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
clear global
%clc % clear command line to avoid funy behavior (sometimes command line keeps printing lines continuoulsy
%min max slider
mins=0;maxs=1;

%% 1 gene or 1 phenotype = 1 matrix
% Preprocessing realized in file "PreProcessData.m"
% choose if you want raw data / resized data or Filtered Data

% load('Patt_CellProp_Rescaled_2017_sept_06th_FormulaTitle_raw_data')% raw image
% load('Patt_CellProp_Rescaled_2017_sept_06th_test_resized_data') %resize + bicubiv
%
load('Patt_CellProp_Rescaled_2017_sept_06th_test_filtered_data')% gaussian filter
%load('Patt_CellProp_Rescaled_2017_sept_06th_test_filtered_discretized')
%load('/Users/Lorette/Google Drive/YohannsBellaiche/z_var/Gene_processed_SR_22nov2017.mat')

%% Which mask to apply
LargMask  =(mask_notum==0);%.*(wg==0).*(bi==0).*(hth==0).*(eyg==0).*(Sr==0);
LargMask  =zeros(size(mask_notum));
IndBckgrnd=find(LargMask);


%% Define function : Y=a1*Gene1 + a2*Gene2...
% functions used in slider for update (callback)

%function for delamination
myfunc_Promo = @(a1,a2,a3,a4,a5,a6,a7) a1*esg+a2*esg.*hth+a3*hth;
myfunc_Inhib = @(a1,a2,a3,a4,a5,a6,a7) a5*sd+a6*h+a7*bi;

%text for figure title
frmt='%.2f';
text=@(a1,a2,a3,a4,a5,a6,a7) strcat(...
    num2str(a1,frmt),' esg +  ',num2str(a2,frmt),'esg*hth +',num2str(a3,frmt),'hth',...
    num2str(a5,frmt),' sd  ',num2str(a6,frmt),'h  ',num2str(a7,frmt),'bi  ');
%text_Promo=@(a1,a2,a3,a4,a5,a6,a7) strcat('Promotes: ',num2str(a1,frmt),' esg +  ',num2str(a2,frmt),'esg*hth');% + ',num2str(a3,frmt),'esg*h');  + ',num2str(a4,frmt),'hth*bi');
%text_Inhib=@(a1,a2,a3,a4,a5,a6,a7) strcat('Inhibits: ',num2str(a5,frmt),' sd  ',num2str(a6,frmt),'bi*h  ');

%% PLOT
% initialize parameter
a1=0;a2=0;a3=0;a4=0;a5=0;a6=0;a7=0;
alpha1=1;alpha2=1;alpha3=1;

% initialize figure
[r,c]=size(BH1);aspectratioimage =c/r;scale = 0.91;
l=560*1.8 ;
w=460;

hFig=figure('Name','Delamination','Visible', 'on','units','pixel','Position',[300 300 l w]);

hAxe=gca;
imagesc(cat(3,delam,zeros(r,c),zeros(r,c)));%LargMask==0)
set(hAxe,'xtick',[]); set(hAxe,'ytick',[]);
posFig=hFig.Position;
posAx=hAxe.Position;


set(hAxe,'units','pixel','Position',[40 25 420*scale*aspectratioimage   420*scale]);
hAxe.Units='normalized';

mytitle=  @(a1,a2,a3,a4,a5,a6,a7) uicontrol(hFig,'style','text','Units','normalized','position',[0.01,0.95,0.45,0.05],'String',text(a1,a2,a3,a4,a5,a6,a7),'Fontsize',13);
mytitle(0,0,0,0,0,0,0);

% opacity alpha1
slid_blend1 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.50,0.83,0.1,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.1 0.10],'callback' , 'alpha1= get(slid_blend1,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
set(slid_blend1,'Value',get(slid_blend1,'Max'))
txtalpha1m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.46,0.9,0.04,0.03],'String',num2str(mins));
txtalpha1M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.6,0.90,0.04,0.03],'String',num2str(maxs));
txtalpha1Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.52,0.94,0.08,0.033],'String','Delam ');
% opacity alpha2
slid_blend2 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.66,0.83,0.1,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.1 0.10],'callback' , 'alpha2= get(slid_blend2,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
set(slid_blend2,'Value',get(slid_blend2,'Max'))
txtalpha2m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.62,0.90,0.04,0.03],'String',num2str(mins));
txtalpha2M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.76,0.90,0.04,0.03],'String',num2str(maxs));
txtalpha2Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.68,0.94,0.08,0.033],'String','Promotor');

% opacity alpha3
slid_blend3 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.80,0.83,0.1,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.1 0.10],'callback' , 'alpha3= get(slid_blend3,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);;mytitle2(a1,a2,a3,a4,a5,a6,a7)');
set(slid_blend3,'Value',get(slid_blend3,'Max'))
txtalpha2m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.76,0.90,0.04,0.03],'String',num2str(mins));
txtalpha2M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.90,0.90,0.04,0.03],'String',num2str(maxs));
txtalpha2Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.82,0.94,0.08,0.033],'String','Inhibitor');

%% Promotion
% slider 1: esg
slid1 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.73,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a1= get(slid1,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt1m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.80,0.04,0.03],'String',num2str(mins));
txt1M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.80,0.04,0.03],'String',num2str(maxs));
txt1Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.84,0.04,0.033],'String','esg');

% slider 2: esg*h
slid2 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.63,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a2= get(slid2,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt2m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.7,0.04,0.04],'String',num2str(mins));
txt2M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.7,0.04,0.04],'String',num2str(maxs));
txt2Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.74,0.04,0.04],'String','esg*h');

% slider 3: 
slid3 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.53,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10],'callback' , 'a3= get(slid3,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt3m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.6,0.04,0.04],'String',num2str(mins));
txt3M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.6,0.04,0.04],'String',num2str(maxs));
txt3Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.64,0.04,0.04],'String','hth');
%{
%slider 4: 
slid4 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.43,0.22,0.1],'Min' , mins , 'Max' , maxs , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a4= get(slid4,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt4m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.5,0.04,0.04],'String',num2str(mins));
txt4M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.5,0.04,0.04],'String',num2str(maxs));
txt4Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.54,0.04,0.04],'String','-');
%}
%% Inhibits 


%slider 5: sd
slid5 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.33,0.22,0.1],'Min' , -1 , 'Max' , 0 , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a5= get(slid5,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt5m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.4,0.04,0.04],'String',num2str(-1));
txt5M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.4,0.04,0.04],'String',num2str(0));
txt5Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.44,0.04,0.04],'String','sd');
set(slid5,'Value',get(slid5,'Max'))


%slider 6: bi*h
slid6 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.23,0.22,0.1],'Min' , -1 , 'Max' , 0 , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a6= get(slid6,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt6m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.3,0.04,0.04],'String',num2str(-1));
txt6M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.3,0.04,0.04],'String',num2str(0));
txt6Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.34,0.04,0.04],'String','h');
set(slid6,'Value',get(slid6,'Max'))

%
% slider 7: eyg
slid7 = uicontrol(hFig,'style','slider','Units','normalized','position',  [0.70,0.13,0.22,0.1],'Min' , -1 , 'Max' , 0 , ...
    'SliderStep',[0.05 0.10], 'callback' , 'a7= get(slid7,''value'');ImResult=Modif6(myfunc_Promo(a1,a2,a3,a4,a5,a6,a7),myfunc_Inhib(a1,a2,a3,a4,a5,a6,a7),IndBckgrnd,delam,alpha1,alpha2,alpha3);mytitle(a1,a2,a3,a4,a5,a6,a7);');
txt7m = uicontrol(hFig,'style','text','Units','normalized','position',    [0.66,0.2,0.04,0.04],'String',num2str(-1));
txt7M = uicontrol(hFig,'style','text','Units','normalized','position',    [0.92,0.2,0.04,0.04],'String',num2str(0));
txt7Title = uicontrol(hFig,'style','text','Units','normalized','position',[0.80,0.24,0.04,0.04],'String','bi');
set(slid7,'Value',get(slid7,'Max'))

%}

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
