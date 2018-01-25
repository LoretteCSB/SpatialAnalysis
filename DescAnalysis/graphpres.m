clear
close all
load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesBiv/';
fileGraphinterm='intBivariatePlot';


minAR = min(AreaRatio_vect'); %
mask = minAR==1;
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
ncol=2;nrow=3;
fileGraphinterm ='intbivariateplot';
for ii=1:1
    p=0;figure;nfig=1
    fileGraph=['BivariatePlot_cutoff_',list_phys{ii}];
    if ii==1
        mylim=[0,4];
    else 
        mylim=[0,1.5]
    end
    for jj= 5:5
        
        x=Gene_vect(mask==1,jj);
        y=Phys_vect(mask==1,ii);
        
        %order point for plotting
        [x,is]=sort(x);
        y = y(is);
        
       
        scatter(x,y)
        
        mymethod = 'moving';%'lowess'
        ys = smooth(x,y,0.1,mymethod);
        %{
        hold on
        plot(x,ys)
        %xlabel(list_gene(jj))
        %ylabel()
        plot([cutoff_gene(jj),cutoff_gene(jj)],[0 ,8])
        %ylim(mylim)
        hold off
         %}       
    end
    
end
clear ii jj x y file* ncol nrow p nfig index exec*
