clear
close all
load('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')%,'data*')


dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/StatDesBiv/';
fileGraphinterm='intBivariatePlot';


minAR = min(AreaRatio_vect'); %
mask = minAR==1;

%% bivariate plot (full block)

%{
ncol=2;nrow=3;
fileGraphinterm ='bivariateplot';
for ii=1:2
    p=0;figure;nfig=1
    fileGraph=['BivariatePlot',list_phys{ii}];
    
    for jj= 1:length(list_gene)
        
        x=Gene_vect(mask==1,jj);
        y=Phys_vect(mask==1,ii);
        
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
        xlabel(list_gene(jj))
        ylabel(list_phys(ii))
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
    
    exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'_',list_phys{ii},'.pdf ',fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    cd(dirnow);
    close all
end
clear ii jj x y file* ncol nrow p nfig index exec*

%}

%% box plot 2 cat
%{
ncol=2;nrow=3;
fileGraphinterm ='tempboxplot';

DelamBin = Phys_vect(:,2)>0;
p=0;figure;nfig=1
fileGraph=['BoxPlot_Delam01_gene'];

for jj= 1:length(list_gene)
    
    x=Gene_vect(mask==1,jj);
    y=DelamBin(mask==1);
    
    % organize subplot
    p=p+1;
    if(p>ncol*nrow)
        p=1;figure;nfig=nfig+1;
    end
    subplot(nrow,ncol,p)
    boxplot(x,y)
    set(gca,'XTickLabel',{'no delam','>0'})
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

exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'_','.pdf ',fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
close all

clear ii jj x y file* ncol nrow p nfig index exec*
%}

%% bivariate plot with delamination categorize in 3 classes
%{
y=discretize(Phys_vect(mask==1,2),[0,0.00000001,0.5,3],[0,0.16,0.67]);
%median(Phys_vect(Phys_vect(:,2)>0&Phys_vect(:,2)<=0.5,2)) 
ncol=2;nrow=3;
fileGraphinterm ='intbivariateplot';

p=0;figure;nfig=1
fileGraph=['BivariatePlot_DelamDiscrete'];

for jj= 1:length(list_gene)
    
    x=Gene_vect(mask==1,jj);
   
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
    
    xlabel(list_gene(jj))
    ylabel('delam')
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

exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'_','.pdf ',fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
eval(exec_Terminal);
cd(dirnow);
close all
clear ii jj x y file* ncol nrow p nfig index exec*

%}

%% bivariate plot with threshold
%{
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
ncol=2;nrow=3;
fileGraphinterm ='intbivariateplot';
for ii=1:2
    p=0;figure;nfig=1
    fileGraph=['BivariatePlot_cutoff_',list_phys{ii}];
    if ii==1
        mylim=[0,4];
    else 
        mylim=[0,1.5]
    end
    for jj= 1:length(list_gene)
        
        x=Gene_vect(mask==1,jj);
        y=Phys_vect(mask==1,ii);
        
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
        xlabel(list_gene(jj))
        ylabel(list_phys(ii))
        plot([cutoff_gene(jj),cutoff_gene(jj)],[0 ,8])
        ylim(mylim)
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
    
    exec_Terminal=['! "/System/Library/Automator/Combine PDF Pages.action/Contents/Resources/join.py" -o ',fileGraph,'.pdf ',fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    exec_Terminal=['! rm ', fileGraphinterm,'*.pdf'];
    eval(exec_Terminal);
    cd(dirnow);
    close all
end
clear ii jj x y file* ncol nrow p nfig index exec*

%}
%% bivariate plot with threshold
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
ncol=2;nrow=3;
close all
fileGraphinterm ='intbivariateplot';
for ii=1:1
    p=0;figure;nfig=1
    fileGraph=['BivariatePlot_cutoff_',list_phys{ii}];
    if ii==1
        mylim=[0,4];
    else 
        mylim=[0,1.5]
    end
    for jj= 1:length(list_gene)
        
        x=Gene_vect(mask==1,jj);
        y=Phys_vect(mask==1,ii);
        
        %order point for plotting
        [x,is]=sort(x);
        y = y(is);
        
        figure
        scatter(x,y)
        
        mymethod = 'moving';%'lowess'
        ys = smooth(x,y,0.1,mymethod);
        hold on
        plot(x,ys,'LineWidth',2)
        xlabel(list_gene(jj))
        ylabel(list_phys(ii))
       % plot([cutoff_gene(jj),cutoff_gene(jj)],[0 ,8])
        ylim(mylim)
        hold off
        
    end
    
end
clear ii jj x y file* ncol nrow p nfig index exec*