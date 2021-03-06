clear
close all

dirdata='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/';
filedata ='GeneAndPhys_Grid_Steph20180201';

load([dirdata,filedata])
clear Folder* dirdata fileout

dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/Model/';
fileout='ExperimentSummaryDivision.txt';

type_model='lasso';


%% SET UP MODEL PARAMETERS
list_gene_excl = {'Delta','hh','wg'};%sd,wg
meth_lambda    = 'MaxFeature'; %{'1SE';'MaxFeature'}';
max_feature    = 10;
MethGene01     = 'Manual';   %'Manual','Automatic'
Interact01     = 0;


% ---------- Gene Activated or not
switch MethGene01
    case 'Automatic'
        Gene_bin = MaskSR_vect;
    case 'Manual'
        Gene_bin = MaskYB_vect;
    otherwise
        for g=1:length(list_gene)
            cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
            Gene_bin(:,g)=Gene_vect(:,g)>cutoff_gene(g)+0;
        end
end


% EXCLUDE GENE and Apply mask
[list_gene_incl,ia]=setdiff(list_gene,list_gene_excl);

% Most restrictive mask : only keep point if area ratio 1
minAR = min(AreaRatio_vect(:,ia),[],2); % min per row
mask = (minAR==1) .* (Phys_vect(:,3)==1);% 1 if area ratio =1 and y has a value

X = Gene_bin(mask==1,ia)+0;
y = Phys_vect(mask==1,1)+0;

% Create interaction terms
if Interact01==1
    X=x2fx(X,'interaction');
    list_gene_interact=x2fx_names(list_gene_incl,'interaction');
    % remove column corresponding to constant
    X=X(:,2:end);
    list_gene_interact=list_gene_interact(2:end);
else
    list_gene_interact=list_gene_incl;
end

%% split sample training/test
myseed =1;
rand('seed', myseed);
cv = cvpartition(length(y),'holdout',0.15);

% Test set
Xtest = X(test(cv),:);
ytest = y(test(cv),:);

% Training and validation set
Xtrain = X(training(cv),:);
ytrain = y(training(cv),:);


%% Generate models for various prenalty parameter lambda
rand('seed', 19479);

[B,FitInfo] = lasso(Xtrain,ytrain,'CV',10);%,'PredictorNames',list_gene_incl);
[ax,figh]=lassoPlot(B,FitInfo,'PlotType','CV');%  green circle ==> |Lambda| with minimum cross-validation error.

if ismember(meth_lambda,'MaxFeature')
    n_feature_per_model =sum(B~=0);
    idx_max_feature = find(n_feature_per_model<=max_feature);
    [~,temp]=min(FitInfo.MSE(:,idx_max_feature));
    idx_best_lambda = idx_max_feature(temp);
    
else
    idx_best_lambda = FitInfo.(['Index',meth_lambda{:}]);
end

% variable selected
WeightOpti = B(:,idx_best_lambda);% best lambda can be 'MinMSE' or '1SE'
ListVarOpti= list_gene_interact(WeightOpti~=0);
TableParaOpti=table(list_gene_interact',WeightOpti,'VariableNames',{'name','weight'});
TableParaOpti(TableParaOpti.weight==0,:)=[];
n_var_selected = length(ListVarOpti);

% Predicted values
ytrain_pred = Xtrain*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);
ytest_pred  =  Xtest*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);
y_pred      =      X*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);

% R2 and R2adj
[R2train,R2adjtrain] = myR2(ytrain,ytrain_pred,n_var_selected);%train set

pct_err_train          =( ytrain(:)-  ytrain_pred(:))./ytrain(:);
pct_err_test           =( ytest(:)-  ytest_pred(:))./ytest(:);
pct_err_train_within10 =length(find(pct_err_train>=-0.10&pct_err_train<=0.10))/length(pct_err_train);
pct_err_train_within20 =length(find(pct_err_train>=-0.20&pct_err_train<=0.20))/length(pct_err_train);
pct_err_test_within10  =length(find(pct_err_test>=-0.10&pct_err_test<=0.10))/length(pct_err_test);
pct_err_test_within20  =length(find(pct_err_test>=-0.20&pct_err_test<=0.20))/length(pct_err_test);


clear h



%% Diagnostic plot

h=figure

subplot(2,2,1)% histogram error - training set

histogram(ytrain-  ytrain_pred) 
title(['Error - training set' ])


subplot(2,2,2)% obs vs predict  (training) 
scatter(ytrain,ytrain_pred) 
hold on
plot([0,max(ytrain)],[0,max(ytrain)],'LineWidth',1.5)
plot(sort(ytrain),sort(ytrain)*0.8,'k--',sort(ytrain),sort(ytrain)*1.2,'k--','LineWidth',1.5)
hold off
title(['Obs versus Pred (training)(R2 = ',num2str(round(R2train,2)),')'])

subplot(2,2,3)% obs vs predict  (test)
scatter(ytest(:),ytest_pred(:)) 
hold on
plot([0,max(ytest)],[0,max(ytest)],'LineWidth',1.5)
plot(sort(ytest),sort(ytest)*0.8,'k--',sort(ytest),sort(ytest)*1.2,'k--','LineWidth',1.5)

hold off
title('Observed versus Predicted (test)')


subplot(2,2,4)% Overall image
y_pred = NaN(length(Phys_vect),1);
y_pred(mask==1) = X*B(:,idx_best_lambda)+FitInfo.Intercept(idx_best_lambda);

imagesc(reshape(y_pred,size(data_phys.div))); 
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Predicted Map')
colorbar




%%  ****************** Additional functions **************************

% calcul R2
function [R2,R2adj] = myR2(y,ypred,k)
n = length(y);
idx_nan = isnan(y);
SStot = var(y)*n;
SSres = sum( (y(:)-ypred(:)).^2 );
R2    = 1 - SSres/SStot;
R2adj = 1-(1-R2)*(n-1)/(n-k-1);
end

% file exchange
function res=x2fx_names(PredictorNames,model)
%   returns a cell array of "compound" predictor variable names that
%   correpond to the x2fx'ed matrix of "normal" predictor variable names.


%set up a table of primes (don't feel like calling prime() every time)
prime_table = [...
    2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,...
    127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229];  %that's primes(230) for ya. fifty of them.

%check if have enough primes
%can probably do that n/ln(n) thing to figure out how many primes we actually need, but 50 should be enough for anybody (just like 640KB of memory)
if length(prime_table)<length(PredictorNames)
    error('internal error: increase size of prime table (really mixing more than 50 variables???)')
end

%run the x2fx to mix a row of prime numbers
mix = x2fx(prime_table(1:length(PredictorNames)),model);

%factor each number in the mix to figure out which variable is mixed where
for i = 1:length(mix)
    %the following should work -- it checks if x2fx has returned a constant
    if mix(i)==1
        res{i} = 'constant';
        continue;
    end
    %figure out the prime factors
    factors = factor(mix(i));
    
    %figure out unique factor list so that we return powers: instead of
    % the ugly A*A*A*B we return prettier A^3+B
    u_factors = unique(factors);
    
    %rinse and repeat
    for factor_cnt = 1:length(u_factors)
        factor_power = sum(u_factors(factor_cnt)==factors);
        factor_name = PredictorNames{find(u_factors(factor_cnt) == prime_table,1)};
        if (factor_power)>1
            factor_name = [factor_name,'^',num2str(factor_power)];
        end
        if (factor_cnt == 1)
            res{i} = factor_name;
        else
            res{i} = [res{i},'*',factor_name];
        end
    end
end
end
