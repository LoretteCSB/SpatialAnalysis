% ## aim :  find cluster of expression among all the genes
% number of clusters may be helpul to determine regions of interest
% gene expression is normalized (all genes vary between 0 and 1)

clear
close all

%% Import file names
cd('/Users/Lorette/Google Drive/pour_lorette')

k = struct2table(dir('*.tif'));
filenames = k.name(k.bytes>0);
clear k

%% Create matrix with 1 column = 1 gene
nfile= length(filenames);
X=[];

for f =1:nfile
    [mypic,map]=imread(filenames{f});
    mypic=double(mypic);
    mypic=mypic/max(max(mypic)); % normalize the color : all the genes varies between 0 and 1
    [r,c]=size(mypic);
    mypic=reshape(mypic,[r*c,1]);
    X=[X,mypic];
end
clear mypic f


%% remove unwanted 0

ix_row0 = max(X')== 0;
X=double(X(ix_row0==0,:));

%% Clustering
% use K-means, with different numbers of cluster
% Note: Hierarchical classif not adapted ==> too many points

mkdir('ResultCluster')
numClusters = 12;
for cl=6:numClusters
    C = kmeans(X, numClusters);
    picC = zeros(r*c,1);
    picC(ix_row0==0)=C;
    picC=reshape(picC/cl*100,[r,c]);
    imwrite(picC,jet,['ResultCluster/clusters',num2str(cl),' ( ',num2str(nfile),'normalized genes).tif'])
end

