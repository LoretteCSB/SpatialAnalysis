%
clear
close all
load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StudyThreshold/';
minAR = min(AreaRatio_vect'); %
mask = minAR==1;


%% Binarize data (rough cutoff selected manually)
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
for g=1:length(list_gene)
    Gene_bin(:,g)=Gene_vect(:,g)>cutoff_gene(g)+0;
end

%% Plot gene binarized
%{
p=0;figure;nfig=1
    
ncol=2;nrow=3;
fileGraph='BinarizedGene';
fileGraphinterm ='binarizedGene';
for ii=1:length(list_gene)
    % organize subplot
    p=p+1;
    if(p>ncol*nrow)
        p=1;figure;nfig=nfig+1;
    end
    subplot(nrow,ncol,p)
    imagesc(reshape(Gene_bin(:,ii),size(data_phys.div)))
    set(gca,'xtick',[]);set(gca,'ytick',[])
    title(list_gene(ii))
    
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
        case num2cell(50:59)
            index='f';
        otherwise
            nfig='g';
    end
    
    print(strcat(dirout,fileGraphinterm,index,num2str(round(nfig,0))),'-dpdf')
    
end
dirnow=pwd;
cd(dirout);

exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'.pdf ',fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
close all
%}