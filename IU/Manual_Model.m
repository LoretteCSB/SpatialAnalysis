%% Explore models manually
% Y = phenotype variable (delamination or proliferation)
% X = gene expression

% delamination = + a1*esg  + a2*(esg*h) - a3*sd + a4*(bi*h)
% proliferation
% would be better to do with slider, but can make it work for now

% Waiting for Boris map to calculate r2
% open question: should I work with raw map or filered?
clear all
close all
mydir=  '/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/Maria2Lorette_6-09_test/';
dirmask='test_2017_sept_6th_2_mask_without BH1/';
dirdata='Patt_CellProp_Rescaled_2017_sept_06th_test_raw_data/';
dirphysio='ImageBoris/'


%% Parameters
nrmlz01 =1;% dummy to normalize gene expression

nrow = 4;
ncol = 3;
filename='linear_eq';
plot_yes_no=1;

%% Mask notum
mask =importimage([mydir,dirmask],'Mask_Notum_2.tif',nrmlz01);

%% Y
%{
delam  = importimage([mydir,dirphysio],'Delamination_BIGwt2_nDelMax=4_0015_NEG.png',0);
prolif = importimage([mydir,dirphysio],'Proliferation_BIGwt2_nDivMax=7_0015_NEG.png',0);

delam = double(rgb2gray(delam));
prolif= double(rgb2gray(prolif));
%delam  =importimage([mydir,dirp],'z_Delamination_12h.tif',0);%prolif =importimage([mydir,dirdata],'z_Proliferation_12h.tif',0);

% problem with Boris map. I cannot use the one Maria gave me
% Need Boris to give me a map without marocheat (background), cell border,
% can he recode the background 999999?
% cell border = recoded as background or value
%unique(delam) %  [0;26;64;128;192;255]
%unique(prolif) % [0;26;37;73;110;146;183;219;255]
% apply mask to remove writing
%imfinfo([mydir,dirphysio,'Delamination_BIGwt2_nDelMax=4_0015_NEG.png'])

delam(double(mask)==0,:)=NaN; %delam  = delam/max(max(delam));
prolif(mask==0,:)=NaN;%prolif = prolif/max(max(prolif));

%{
f=figure
hist(delam)
title('Histogram values in delamination image')
saveas(f,'delamination')
f=figure
hist(prolif)
title('Histogram values in proliferation image')
saveas(f,'prolif')
%}
%}
%% X
%[BH1,map]=imread([mydir,dirdata,'BH1_12apf_6.tif']);
esg = importimage([mydir,dirdata],'esg_12apf_4.tif',nrmlz01);
bi  =importimage([mydir,dirdata],'bi_12apf_2.tif',nrmlz01);
sd  =importimage([mydir,dirdata],'sd_12apf_7.tif',nrmlz01);
hth =importimage([mydir,dirdata],'hth_12apf_4.tif',nrmlz01);
DIAP1 =importimage([mydir,dirdata],'DIAP1_12apf_3.tif',nrmlz01);
h =importimage([mydir,dirdata],'h_12apf_5.tif',nrmlz01);
wg =importimage([mydir,dirdata],'wg_12apf_7.tif',nrmlz01);

% missing tsh (Maria is going to redo the experiment, keep for later)
% and h for delamination. DIAP1?
% missing wg, via and fz3 (not expressed at 12apf) for proliferation

% background set to NaN so it can be more easily identified
%{
esg = setBckgrndNaN(esg,mask);
bi  = setBckgrndNaN(bi,mask);
sd  = setBckgrndNaN(sd,mask);
hth = setBckgrndNaN(hth,mask);
DIAP1 = setBckgrndNaN(DIAP1,mask);
h   = setBckgrndNaN(h,mask);
wg  = setBckgrndNaN(wg,mask);
%}

%n1=esg/max(esg(:));
%n2=esg/prctile(esg(:),99);


%% Try various models
IndBckgrnd= ((esg.*h.*sd.*bi)==0)&(mask==0);
ModelEq = @(a1,a2,a3,a4) a1*esg+a2*(esg.*h)-a3*sd-a4*(bi.*h);%delamination

Model=ModelEq(0.5,0,0.5,0);
ModelEq = @(a1,a2,a3,a4) -a1*(wg.*fz3)-a2*(hth.*bi)-a3*(eyg.*bi)-a4*(sr);%prolif
Model(IndBckgrnd==1)=NaN;
%nrmlize image so it will be on a [0-1] scale
Model = (Model-nanmin(Model(:)))/(nanmax(Model(:))-nanmin(Model(:))); % nrmlz so that
Model(IndBckgrnd==1)=NaN;
fig=figure
imagesc(Model)
colorbar
%
p=0;nfig=1;
fig=figure
for a1 =0:0.5:1
    for a2 =0:0.5:1
        for a3 =0:0.5:1
            for a4 =0:0.5:1
                if (a1|a2|a3|a4)
                    Model=ModelEq(a1,a2,a3,a4);
                    Model(IndBckgrnd==1)=NaN;
                    
                    %nrmlize image so it will be on a [0-1] scale
                    Model = (Model-nanmin(Model(:)))/(nanmax(Model(:))-nanmin(Model(:))); % nrmlz so that
                    Model(IndBckgrnd==1)=NaN;
                    
                    if (p==nrow*ncol)
                        tightfig
                        print('-fillpage',strcat(filename,index,num2str(round(nfig,0))),'-dpdf')
                        nfig=nfig+1;%for label pdf
                        p=0;
                        fig=figure
                        if(~plot_yes_no)
                            set(fig,'visible','off');
                        end
                    end
                    p=p+1;
                    subplot(nrow,ncol,p)
                    imagesc(Model)
                    set(gca,'XTick',[])
                    set(gca,'YTick',[])
                  %  title(strcat(num2str(a1),'*esg + ',num2str(a2),'*(esg*hth) - ',num2str(a3),'*sd'),'fontsize',10);
                     title(strcat('a1=',num2str(a1),' a2= ',num2str(a2),' a3=',num2str(a3),' a4=',num2str(a4)),'fontsize',10);
                    
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
                    %-bestfit'-fillpage'
                    %print('-fillpage',strcat(filename,index,num2str(round(nfig,0))),'-dpdf')
                end
            end
        end
    end
end
%

function y = importimage(dirpath,file,nrmlz01)
y=imread([dirpath,file]);
y=double(y);
if nrmlz01
    y = y/max(y(:));
end
end


function y = setBckgrndNaN(y,mask)
ind_bckg=(y==0)&(mask==0);
y(ind_bckg)=NaN;
end
