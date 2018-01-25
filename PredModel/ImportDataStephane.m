% Import et pr?-analyse des 24 genes

clear all

% locate data
TopFolder = '/Users/Lorette/Google Drive/YohannsBellaiche/stephane/AOA_Output'
listSubFolder  = dir(TopFolder);

% select gene folder (not automatic)
FolderPhysio = listSubFolder(2,:);
FolderGene = listSubFolder([1,3:24],:);

clear listSubFolder

% Import Gene expression data
for f=1:length(FolderGene)
    load([TopFolder,'/',FolderGene(f).name,'/AOA_backup.mat'],'gene','AreaRatios');  
    pos_undscr = regexp(FolderGene(f).name,'_');
    data_gene(f).name = FolderGene(f).name(pos_undscr+1:end);
    data_gene(f).gene = gene;
    data_gene(f).ar   = AreaRatios;
end
clear f gene AreaRatios pos_undscr 

% Import Boris's physiological data
load([TopFolder,'/',FolderPhysio.name,'/AOA_backup.mat'],'dnD','dnA','AreaRatios');  
%{
data_phys(1).name='div';
data_phys(1).mean_nb=dnD;

data_phys(2).name='delam';
data_phys(2).mean_nb=dnA;
%}
data_phys.div  = dnD/2;
data_phys.delam= dnA;
data_phys.ar   = AreaRatios;


clear dnD dnA AreaRatios TopFolder


%% version vectorized
Gene_vect=[];
AreaRatio_vect=[];
[r,c]=size(data_gene(1).ar);
for f =1:length(data_gene)
    Gene_vect    =[Gene_vect,reshape(data_gene(f).gene,[r*c,1])];
    AreaRatio_vect=[AreaRatio_vect,reshape(data_gene(f).ar,[r*c,1])];
end
list_gene  ={data_gene.name};
%Gene_vect = array2table(Gene_vect,'VariableNames',{data_gene.name});
%AreaRatio_vect = array2table(AreaRatio_vect,'VariableNames',{data_gene.name});
clear f c r 

Phys_vect=[data_phys.div(:),data_phys.delam(:),data_phys.ar(:)];
list_phys={'div','delam','ar'};


% save output
save('/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/GeneAndPhys_Grid_Steph171222')