% Should we transform some gene (log, sqrt...)
% What is the basic correlation between genes/genes, and gene/physio
% can I define a threshold to dichotimize variables?
% How much are the pattern overlapping between 2 genes?
% Univariate prediction

clear
close all

dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesUniv/';

load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


%% Mask 
% mask based on intersection notum 
minAR = min(AreaRatio_vect'); % 
mask_intersect_notum = (minAR==1)+0; 

% mask based on area ratio for physio
mask1 = (data_phys.ar ==1)+0;
%% Comparison impact mask

%{


%}
fig = figure
subplot(2,1,1)
h=histogram(Phys_vect(mask1(:)==1,1),'Normalization','pdf');
poscenter = h.BinEdges(1:end-1) + diff(h.BinEdges) / 2;
edges = h.BinEdges;
y=h.BinCounts;
h2=histogram(Phys_vect(mask_intersect_notum(:)==1,1),edges,'Normalization','pdf');
y2=h2.BinCounts;
plot(poscenter,y,'k','LineWidth',1.5)
hold on
plot(poscenter,y2,'r','LineWidth',1.5)  
legend({'ar physio==1','intersect notum : min(ar gene) ==1'})
title('Histogram Division ')

subplot(2,1,2)
h=histogram(Phys_vect(mask1(:)==1,2),'Normalization','pdf');
poscenter = h.BinEdges(1:end-1) + diff(h.BinEdges) / 2;
edges = h.BinEdges;
y=h.BinCounts;
h2=histogram(Phys_vect(mask_intersect_notum(:)==1,2),edges,'Normalization','pdf');
y2=h2.BinCounts;
plot(poscenter,y,'k','LineWidth',1.5)
hold on
plot(poscenter,y2,'r','LineWidth',1.5)  
legend({'ar physio==1','intersect notum : min(ar gene) ==1'})
title('Histogram Delamination ')

saveas(fig,[dirout,'Compare_hist_MaskAR1_MaxIntersectNotum.pdf'])

clear h h2 fig y y2
%%
h=figure
subplot(2,1,1)
imagesc(data_phys.ar ==1)
set(gca,'xtick',[],'ytick',[])
title('Mask area ratio ==1')

subplot(2,1,2)
imagesc(data_phys.ar >=0.95)
set(gca,'xtick',[],'ytick',[])
title('Mask area ratio >= 0.95')

saveas(h,'ChoiceMaskPhysio.pdf')
clear h

%% Visualiser les images

h=figure 
subplot(2,2,1)
imagesc(data_phys.div)
title('Mean nb of divisions per hour per box')
set(gca,'xtick',[],'ytick',[])

subplot(2,2,2)
histogram(data_phys.div,0:0.5:10)
title('0.5 bin')

subplot(2,2,3)
imagesc(data_phys.delam)
title('Mean nb of delamination per hour per box')
set(gca,'xtick',[],'ytick',[])

subplot(2,2,4)
histogram(data_phys.delam,0:0.05:2)
title('0.05 bin')

saveas(h,[dirout,'Fig_ImageHistPhysio.pdf'])
clear h




%% Visualiser les boxplots
mask = (data_phys.ar ==1)+0;

h=figure 

subplot(2,1,1)
boxplot(data_phys.div(mask==1))
title('Mean nb of division per hour per box')
cutoff_outlier_div=prctile(data_phys.div((mask==1)),75)+1.5*iqr(data_phys.div((mask==1)));% q3+1.5(q3-q1)
text(1.05,cutoff_outlier_div,['outlier ',num2str(round(cutoff_outlier_div,2))]);

subplot(2,1,2)
boxplot(data_phys.delam((mask==1)))
title('Mean nb of delamination per hour per box')
cutoff_outlier_delam=prctile(data_phys.delam((mask==1)),75)+1.5*iqr(data_phys.delam((mask==1)));
text(1.05,cutoff_outlier_delam,['outlier ',num2str(round(cutoff_outlier_delam,2))]);

saveas(h,[dirout,'BoxPlotPhysio.pdf'])

clear h

%% where are the outliers

idx_cutoff_div   = (data_phys.div  >= cutoff_outlier_div  ) + 0 ;
idx_cutoff_delam = (data_phys.delam>= cutoff_outlier_delam) + 0 ;

h=figure

subplot(2,1,1)
imagesc(mask+idx_cutoff_div)
title(['Regions with the nb div >',num2str(round(cutoff_outlier_div,2))])

subplot(2,1,2)
imagesc(mask+idx_cutoff_delam)
title(['Regions with the nb delam of delam >',num2str(round(cutoff_outlier_delam,2))])

saveas(h,[dirout,'MapOutliersY.pdf'])

clear h idx*
%% 
%Recoder
histogram(discretize(Phys_vect(mask_phys,1),[0,1,2,10]))
[Div_4cl,b]= discretize(Phys_vect(mask_phys,1),3)

histogram(Div_4cl)


figure 
subplot(2,2,1)

imagesc(mask_phys+1*(data_phys.div<0.5))
title('Cell division')

%%
t=sort(Phys_vect)
mymethod = 'moving';%'lowess'
ys = smooth(t(:,1),t(:,2),0.1,mymethod);
hold on
plot(t(:,1),ys)
xlabel('Cell division')
ylabel('Cell delamination')