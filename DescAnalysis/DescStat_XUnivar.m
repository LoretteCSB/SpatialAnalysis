% Should we transform some gene (log, sqrt...)
% What is the basic correlation between genes/genes, and gene/physio
% can I define a threshold to dichotimize variables?
% How much are the pattern overlapping between 2 genes?
% Univariate prediction

clear
close all
load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


for f=1:length(data_gene)
    stat(f).name  = data_gene(f).name;
    stat(f).n_box_sup60   = sum(sum(data_gene(f).ar>0.60&data_gene(f).ar<1));
    stat(f).n_box_sup70   = sum(sum(data_gene(f).ar>0.70&data_gene(f).ar<1));
    stat(f).n_box_sup80   = sum(sum(data_gene(f).ar>0.80&data_gene(f).ar<1));
    stat(f).n_box_sup90   = sum(sum(data_gene(f).ar>0.90&data_gene(f).ar<1));

    stat(f).n_box_inf1 = sum(sum(data_gene(f).ar<1&data_gene(f).ar>0));
    stat(f).n_box_eq1     = sum(sum(data_gene(f).ar==1));
    stat(f).min   = min(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).q1    = prctile(data_gene(f).gene(data_gene(f).ar>0.90),25);
    stat(f).q5    = prctile(data_gene(f).gene(data_gene(f).ar>0.90),5);
    stat(f).q25   = prctile(data_gene(f).gene(data_gene(f).ar>0.90),25);
    stat(f).q75   = prctile(data_gene(f).gene(data_gene(f).ar>0.90),75);
    stat(f).q95   = prctile(data_gene(f).gene(data_gene(f).ar>0.90),95);
    stat(f).max   = max(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).median= median(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).moy   = mean(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).std   = std(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).skew  = skewness(data_gene(f).gene(data_gene(f).ar>0.90));
    stat(f).kurtosis = kurtosis(data_gene(f).gene(data_gene(f).ar>0.90));

end
stat=struct2table(stat);

%% Which mask should I take?
% All box with AreaRatio =1 across all genes?
% All box with AreaRatio > x across all genes?
 
minAR = min(AreaRatio_vect'); % 
nbox100= sum(minAR==1)% 653
nbox90= sum(minAR>0.9)% 673
nbox80= sum(minAR>0.8)% 690
nbox70= sum(minAR>0.7)% 701

mask = minAR==1; 
%% Visualiser les images
MyImage(data_gene,list_gene)

%% comparaison histogram avec mask et sans mask
% si je prends l'intersection de tous les notums, est-ce que cela change
% qqch ? ma distribution des valeurs?

MyHist(Gene_vect,AreaRatio_vect,mask,list_gene)


% commentaire: 
% - certaines distrib sont bimodales : pas de relation ?vidente avec pattern
% - certaines distrib ont des queues larges ==> quid normalization (?
% v?rifier avec bixplot)

labels={'BH1','DIAP1','Delta','Pax','Sr','Ush','bi','caup','esg','exd','eyg','fz3','h','hh','hth','lgs','mirr','ptc','salm','sd','tsh','usp','wg'}

h=figure
boxplot(Gene_vect,labels,'LabelOrientation','inline')
title('Genetic Expression - before normalization')
saveas(h,'/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/BoxPlotGeneExpr.pdf')
%% Normalization

% normalize by max 1
Gene_nrmlz_max= Gene_vect./max(Gene_vect);
% normalize by 0.99 quantile and set value above 1 to 1 
Gene_nrmlz_99= Gene_vect./prctile(Gene_vect,99);
Gene_nrmlz_99(Gene_nrmlz_99>1,:)=1;



%% visualization image
function MyImage(data_gene,list_gene)
ncol=2;nrow=2;p=0;
p=0;figure;nfig=1

for v=1:length(data_gene)
    p=p+1;
    if(p>ncol*nrow)
        p=1;figure;nfig=nfig+1;
    end
    subplot(nrow,ncol,p)
    imagesc(data_gene(v).gene);
    title(list_gene(v));
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    switch nfig
        case num2cell(0:9)
            index='a';
        case num2cell(10:19)
            index='b';
        case num2cell(20:29)
            index='c';
        case num2cell(30:39)
            index='d';
        case num2cell(40:49)
            index='e';
        otherwise
            nfig='f';
    end
    dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesUniv/';
    fileGraph='Image';

    print(strcat(dirout,fileGraph,index,num2str(round(nfig,0))),'-dpdf')
    
    
end

dirnow=pwd;
cd(dirout);
exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ','Fig_',fileGraph,'.pdf ',fileGraph,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraph,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
end


% Comparison histogram obtained when we focus on mask (intersection)

function MyHist(Gene_vect,AreaRatio_vect,mask,list_gene)
[~,ngene]=size(Gene_vect)
ncol=2;nrow=3;p=0;
p=0;figure;nfig=1

for v=1:ngene
    p=p+1;
    if(p>ncol*nrow)
        p=1;figure;nfig=nfig+1;
    end
    subplot(nrow,ncol,p)
    
    h=histogram(Gene_vect(AreaRatio_vect(:,v)==1,v),'Normalization','pdf');
    poscenter = h.BinEdges(1:end-1) + diff(h.BinEdges) / 2;
    edges = h.BinEdges;
    y=h.BinCounts;
    h2=histogram(Gene_vect(mask==1,v),edges,'Normalization','pdf');
    y2=h2.BinCounts;
    plot(poscenter,y,'k','LineWidth',1.5)
    hold on

    plot(poscenter,y2,'r','LineWidth',1.5)    
  %  legend('All1','After mask')

    title(list_gene(v));
   
    hold off
    switch nfig
        case num2cell(0:9)
            index='a';
        case num2cell(10:19)
            index='b';
        case num2cell(20:29)
            index='c';
        case num2cell(30:39)
            index='d';
        case num2cell(40:49)
            index='e';
        otherwise
            nfig='f';
    end
    dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesUniv/';
    fileGraph='Compare_Hist_Gene1_MaskIntersection1_red';

    print(strcat(dirout,fileGraph,index,num2str(round(nfig,0))),'-dpdf')
    
    
end

dirnow=pwd;
cd(dirout);
exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ','Fig_',fileGraph,'.pdf ',fileGraph,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraph,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
end


function MyImageNormlz(Gene_nrmlz_max,Gene_nrmlz_99,list_gene)

ncol=2;nrow=2;p=0;
p=0;figure;nfig=1

for v=1:length(data_gene)
    p=p+1;
    if(p>ncol*nrow)
        p=1;figure;nfig=nfig+1;
    end
    subplot(nrow,1,p)
    imagesc(data_gene(v).gene);
    title(list_gene(v));
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    switch nfig
        case num2cell(0:9)
            index='a';
        case num2cell(10:19)
            index='b';
        case num2cell(20:29)
            index='c';
        case num2cell(30:39)
            index='d';
        case num2cell(40:49)
            index='e';
        otherwise
            nfig='f';
    end
    dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesUniv/';
    fileGraph='Image';

    print(strcat(dirout,fileGraph,index,num2str(round(nfig,0))),'-dpdf')
    
    
end

dirnow=pwd;
cd(dirout);
exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ','Fig_',fileGraph,'.pdf ',fileGraph,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraph,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
end