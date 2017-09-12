% ## aim :  find cluster of expression among all the genes
% number of clusters may be helpul to determine regions of interest
% gene expression is normalized (all genes vary between 0 and 1)


% Lorette: changer le code pour pouvoir boucler:
% - avec et sans Bh1
% - avec ou sans normalization des donnees
% - avec ou sans le mask notum
% - look at the different options availables for distance...
% scenario de reference: mask / sans bh1 /
clear
close all

currentFolder = pwd;

UseBh1=1;
UseMask=1;
UseNrmlz=1;

for UseBh1 =0:1
    for UseMask=1
        for UseNrmlz=0:1
            
            % suffix for output filename
            if UseMask
                suffixeMask='MaskNotum';
            else
                suffixeMask='NoMask';
            end
            
            if UseBh1
                suffixeBh1='withBh1';
            else
                suffixeBh1='withoutBh1';
            end
            
            if UseNrmlz
                suffixeNrmlz='Nrmlzd';
            else
                suffixeNrmlz='NotNrmlzd';
            end
            
            %% Mask notum
            
            mydir='/Users/Lorette/Google Drive/YohannsBelaich/Spatial analysis/Maria2Lorette_6-09_test/';
            
            %cd([mydir,'test_2017_sept_6th_masks'])% mask was based on 11 genes :
            cd([mydir,'test_2017_sept_6th_2_mask_without BH1'])% Bh1 was excluded
            
            k = struct2table(dir('*.tif'));
            maskfile = k.name(k.bytes>0);            
            ix=regexpi(maskfile,'Mask_Notum*');
            ix=arrayfun(@(x) length(x{:}),ix);         
            maskfile = maskfile(ix==1);
            [masknotum,map]=imread(maskfile{1});         
            masknotum=double(masknotum);
            %imagesc(masknotum);colorbar;
            %% Import filenames (image with gene expression and phenotype)
            
            %rawfiles
            cd('/Users/Lorette/Google Drive/YohannsBelaich/Spatial analysis/Maria2Lorette_6-09_test/Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data')
            
            k = struct2table(dir('*.tif'));
            filenames = k.name(k.bytes>0);
            
            %split gene expression and phenotype (files starting with z)
            ix =regexp(filenames,'z*');
            ix=arrayfun(@(x) length(x{:}),ix);
            zfilename= filenames(ix==1);%phenotypic image
            filenames= filenames(ix==0);%genetic expression image
            
            % exclude Bh1 from clustering
            if UseBh1==0
                ix =regexp(filenames,'BH1*');
                ix=arrayfun(@(x) length(x{:}),ix);
                filenames= filenames(ix==0);%exlude BH1
            end
            
            clear k ix
            
            %% Create matrix with 1 column = 1 gene
            nfile= length(filenames);
            X=[];
            
            for f =1:nfile
                [mypic,map]=imread(filenames{f});
                mypic=double(mypic);
               
                if UseNrmlz
                    mypic=mypic/max(max(mypic)); % normalize : all the genes varies between 0 and 1
                end
                
                if UseMask
                    %figure
                    %image(mypic)
                    %title([filenames{f},'mask'])
                    %colorbar
                    mypic(masknotum==0)=0;
                    
                    % to see impact mask
                    %figure
                    %image(mypic)
                    %title([filenames{f},'mask'])
                    %colorbar
                end
                [r,c]=size(mypic);
                mypic=reshape(mypic,[r*c,1]);
                X=[X,mypic];
            end
            clear mypic f
            
            
            %% identify 0 (background - regions with no expression)
            ix_bckgrnd = (max(X')== 0)'; %identify pixel where all genes = 0
            
            X=double(X(ix_bckgrnd==0,:));% I keep track of the 0
            
            
            %% Clustering
            % use K-means, with different numbers of cluster
            % Note: Hierarchical classif not adapted ==> too many points
            % use silhouette to assess quality
            cd(currentFolder)
            
            mkdir([mydir,'ResultCluster'])
            numClusters = 10;
            for cl=10:numClusters
                cl
                idx = kmeans(X, cl,'MaxIter',1000,'Replicates',1); % repeat 5 times the kmeans with different initilization
                % silh=silhouette(X,idx);
                picC = zeros(r*c,1);
                picC(ix_bckgrnd==0)=idx;
                picC=reshape(picC,[r,c]);
                % figure
                % imagesc(picC)
                %  saveas(gcf,[mydir,'ResultCluster/',num2str(cl),'clusters',suffixeMask,'_',suffixeNrmlz,'_',suffixeBh1],'tif')
                % saveas(gcf,'Figure','png')
              
                imwrite(picC+ones(size(picC)),[ones(1,3);jet(cl)],...% need to add 1 so that color bar is correct, I don't understand why
                    [mydir,'ResultCluster/',num2str(cl),'clusters',suffixeMask,'_',suffixeNrmlz,'_',suffixeBh1,'.tif'])
            end
            
        end
    end
end
%cd(currentFolder)

