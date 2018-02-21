clear
close all
% path to package
%Pas esg/exd
% addpath('/Users/Lorette/Documents/MATLAB/libsvm-3.22/matlab');% libsvm
addpath('/Users/Lorette/Documents/MATLAB/liblinear-2.20/matlab');% liblinear


dirdata='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/';
filedata ='GeneAndPhys_Grid_Steph20180201';

load([dirdata,filedata])
clear Folder* dirdata fileout

dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/Model/';


%% Parameter for simulation

list_gen_inc={'eyg','salm'};fileout='Fig_Delam01_BestModel';

list_gen_inc={'PAX', 'exd','esg'};fileout='Fig_Delam01_ModelExpert';

c=1.00E-02;

Gene_bin = MaskSR_vect;
type_model={'L1_lin_svm'};kernel =5; %L1-regul L2-loss SVM
% type_model={'L1_logistic'};kernel =6;%L1 regul logistic


%%
% select gene
ia = ismember(list_gene,list_gen_inc);

% Most restrictive mask : only keep point if area ratio 1
minAR = min(AreaRatio_vect(:,ia),[],2); % min per row
mask = (minAR==1) .* (Phys_vect(:,3)==1);% 1 if area ratio =1 and y has a value

X = Gene_bin(mask==1,ia)+0;
y = double(Phys_vect(mask==1,2)>0);% column for apoptosis / code 0 if no apoptosis 1 otherwise

% encodage des classes sous la forme -1 et +1
classcode = unique(y);
ysvm = y;
ysvm(y==classcode(1))=-1;
ysvm(y==classcode(2))=1;


%% Check that model produce similar result 
% I do that because the mask is changed when you only select some variables

myseed =1;
rand('seed', myseed);
cv = cvpartition(length(y),'holdout',0.15);

% Test set
Xtest = X(test(cv),:);
ytest = ysvm(test(cv),:);

% temp = Training and validation set
Xtemp = X(training(cv),:);
ytemp = ysvm(training(cv),:);

% now split temp into training and validation
myseed2 =991;
rand('seed', myseed2);
cv2 = cvpartition(length(ytemp),'holdout',0.35);

Xvalid = X(test(cv2),:);
yvalid = ytemp(test(cv2),:);

% temp = Training and validation set
Xtrain = X(training(cv2),:);
ytrain = ytemp(training(cv2),:);

clear *temp

myoption=sprintf('-s %g -c %g -B %g -q',kernel,c,1);% s = model type, c = penalty parameter, b = add a bias term w0, -q quiet mode

model = train(ytrain, sparse(Xtrain), myoption);

% Predicted values
[yvalid_pred, acc_valid, ~] = predict(yvalid, sparse(Xvalid), model);
[ytest_pred , acc_test , ~] = predict(ytest,  sparse(Xtest),  model);

% Quality model
cm_valid = confusionmatrix(yvalid,yvalid_pred)
acc_valid = (cm_valid(1,1)+cm_valid(2,2))/sum(cm_valid(:))
cm_test = confusionmatrix(ytest,ytest_pred)
acc_test = (cm_test(1,1)+cm_test(2,2))/sum(cm_test(:))



%% Plot predict versus 
h=figure

subplot(2,2,1)% Overall image
imagesc(double(data_phys.delam))
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Original Map ')

%-----------
subplot(2,2,2)% Overall image
imagesc(double(data_phys.delam>0))
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Original Map dichotomized ')

%-----------
subplot(2,2,3)% Overall image
model = train(y, sparse(X), myoption);

y_pred = NaN(length(Phys_vect),1);
y_pred(mask==1) = predict(y, sparse(X), model);

imagesc(reshape(y_pred,size(data_phys.div))); 
set(gca,'xtick',[]);set(gca,'ytick',[])
title(list_gen_inc)

%-----------
hsub= subplot(2,2,4)% Overall image

cm = confusionmatrix(y,y_pred(mask==1));

acc=(cm(1,1)+cm(2,2))/sum(cm(:));
cm =array2table(cm,'RowNames',{'Obs0','Obs1'},'VariableNames',{'Pred0','Pred1'})


hAxes = gca;
hAxes.XRuler.Axle.LineStyle = 'none';  
axis off
pos = get(hsub, 'position');

%
hTable = uitable('Data',cm{:,:},'ColumnName',cm.Properties.VariableNames,...
    'RowName',cm.Properties.RowNames,'Units', 'Normalized', 'Position',pos.*[1,0.9,1,0.8]);
cellMaxLen = num2cell(60);
set(hTable, 'ColumnWidth', cellMaxLen);

title(['Confusion Matrix, accuracy = ',num2str(round(acc*100,0)),'%'])
saveas(h,[dirout,fileout],'jpeg')