clear
close all
load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesBiv/';
fileGraphinterm='intBivariatePlot';


minAR = min(AreaRatio_vect'); %
mask = minAR==1;
%{
rho = corr(Gene_vect(mask==1,:));
rho=tril(rho,-1);

% Top 3 corr>0
CompareImage(data_gene,list_gene,'exd',1);
CompareImage(data_gene,list_gene,'hth',2);
CompareImage(data_gene,list_gene,'ptc',3);
CompareImage(data_gene,list_gene,'DIAP1',4);
CompareImage(data_gene,list_gene,'lgs',5);
CompareImage(data_gene,list_gene,'caup',6);


saveas(gcf,strcat(dirout,'TopPosCorrel.pdf'))


% Top 3 corr<0
CompareImage(data_gene,list_gene,'Ush',1);
CompareImage(data_gene,list_gene,'caup',2);
CompareImage(data_gene,list_gene,'eyg',3);
CompareImage(data_gene,list_gene,'salm',4);
CompareImage(data_gene,list_gene,'lgs',5);
CompareImage(data_gene,list_gene,'Ush',6);

saveas(gcf,strcat(dirout,'TopNegCorrel.pdf'))
%}
%% bivariate plot (triangular)
%
fileGraph='BivariatePlot';

ncol=2;nrow=3;p=0;
p=0;figure;nfig=1
fileGraph=[fileGraph,'_all'];
for ii=1:length(list_gene)
    for jj= ii+1:length(list_gene)
        %(ii+1)
        x=Gene_vect(mask==1,ii);
        y=Gene_vect(mask==1,jj);
        
        %order point for plotting
        [x,is]=sort(x);
        y = y(is);
        
        % organize subplot
        p=p+1;
        if(p>ncol*nrow)
            p=1;figure;nfig=nfig+1;
        end
        subplot(nrow,ncol,p)
        scatter(x,y,'.')
        mymethod = 'moving';%'lowess'
        ys = smooth(x,y,0.1,mymethod);
        hold on
        plot(x,ys)
        xlabel(list_gene(ii))
        ylabel(list_gene(jj))
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
            case num2cell(50:59)
                index='f';
            otherwise
                nfig='f';
        end
        
        print(strcat(dirout,fileGraphinterm,index,num2str(round(nfig,0))),'-dpdf')
        
    end
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
%% bivariate plot (full block)
%{
ncol=2;nrow=3;

for ii=1:length(list_gene)
    p=0;figure;nfig=1
    for jj= 1:length(list_gene)
       
        x=Gene_vect(mask==1,ii);
        y=Gene_vect(mask==1,jj);
        
        %order point for plotting
        [x,is]=sort(x);
        y = y(is);
        
        % organize subplot
        p=p+1;
        if(p>ncol*nrow)
            p=1;figure;nfig=nfig+1;
        end
        subplot(nrow,ncol,p)
        scatter(x,y,'.')
        mymethod = 'moving';%'lowess'
        ys = smooth(x,y,0.1,mymethod);
        hold on
        plot(x,ys)
        xlabel(list_gene(ii))
        ylabel(list_gene(jj))
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
            case num2cell(50:59)
                index='f';
            otherwise
                nfig='g';
        end
        
        print(strcat(dirout,fileGraphinterm,index,num2str(round(nfig,0))),'-dpdf')
        
    end
    dirnow=pwd;
    cd(dirout);
 
    exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'_',list_gene{ii},'.pdf ',fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    cd(dirnow);
    close all
end

%}
%% Look up at some perculiar gene relationship
size_grid=size(data_gene(1).gene);


gene='Pax';cutoff=0.08;
h=PlotBinarizedGene(Gene_vect,list_gene,mask,gene,cutoff,size_grid)

gene='exd';cutoff=0.09;
h1=PlotBinarizedGene(Gene_vect,list_gene,mask,gene,cutoff,size_grid)


gene='BH1';cutoff=0.05;
h2=PlotBinarizedGene(Gene_vect,list_gene,mask,gene,cutoff,size_grid)

gene='eyg';cutoff=0.03;
h3=PlotBinarizedGene(Gene_vect,list_gene,mask,gene,cutoff,size_grid)

gene='esg';cutoff=0.025;
h4=PlotBinarizedGene(Gene_vect,list_gene,mask,gene,cutoff,size_grid)




%%

CompareImage(data_gene,list_gene,'usp',1);
CompareImage(data_gene,list_gene,'Ush',2);
CompareImage(data_gene,list_gene,'lgs',3);
CompareImage(data_gene,list_gene,'Ush',4);
CompareImage(data_gene,list_gene,'lgs',5);
CompareImage(data_gene,list_gene,'caup',6);

saveas(gcf,strcat(dirout,'Closeup_scatterlinked.pdf'))




function CompareImage(data_gene,list_gene,gene,pos)
subplot(3,2,pos)
ix_gene =ismember(list_gene,gene);
imagesc(data_gene(ix_gene).gene)
title(gene)
end

function h=PlotBinarizedGene(Gene_vect,list_gene,mask,GeneOfInterest,cutoff,size_grid)
ix_gene=ismember(list_gene,GeneOfInterest);
binrzd_gene= mask'+mask'.*(Gene_vect(:,ix_gene==1)>cutoff)+0;
h=figure
imagesc(reshape(binrzd_gene,size_grid));
set(gca,'xtick',[])
set(gca,'ytick',[])
title(GeneOfInterest)


end