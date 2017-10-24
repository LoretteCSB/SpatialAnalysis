% Gene expression more condensed in nucleus == > bad for segmentation
% interpolation
clear
close all
currentFolder = pwd;
mydir='/Users/Lorette/Google Drive/YohannsBelaich/Spatial analysis/Maria2Lorette_6-09_test/Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';

k = struct2table(dir([mydir,'*.tif']));
filenames = k.name(k.bytes>0);

%split gene expression and phenotype (files starting with z)
ix =regexp(filenames,'z*');
ix=arrayfun(@(x) length(x{:}),ix);
zfilename= filenames(ix==1);%phenotypic image
filenames= filenames(ix==0);%gene

nfile= length(filenames);
X=[];

for f =1:5%:nfile
    ix=regexp(filenames{f},'_');
    [mypic,map]=imread([mydir,filenames{f}]);
    mypic=double(mypic);
    
    figure
    subplot(3,1,1)
    image(mypic)
    title([filenames{f}(1:(ix(1)-1)),' original'])
    
    subplot(3,1,2)
    %figure
    sig =10;
    mypic_filtered=imgaussfilt(mypic,sig);
    image(mypic_filtered)
    title([filenames{f}(1:(ix(1)-1)),' Gaussian ',sig])
    %}
    %{
    figure
    scal =0.08;
    mypic_resize=resizem(mypic,scal,'bicubic');
    image(mypic_resize)
    title([filenames{f}(1:(ix(1)-1)),' Resize ',scal])
    
    %}
    
    subplot(3,1,3)
    %figure
    sig =10;scal =0.08;
    mypic_filtered=imgaussfilt(mypic,sig);
    mypic_filtered_resize=resizem(mypic_filtered,scal,'bicubic');
    image(mypic_filtered_resize)
    title([filenames{f}(1:(ix(1)-1)),'Mix ',sig])
    %}
    
end

figure
    sig =8;scal =0.1;
    mypic_filtered=imgaussfilt(mypic,sig);
    mypic_filtered_resize=resizem(mypic_filtered,scal,'bicubic');
    image(mypic_filtered_resize)
    title([filenames{f}(1:(ix(1)-1)),'Mix ',sig])