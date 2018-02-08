%
clear
close all
dirdata='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/';
filedata ='GeneAndPhys_Grid_Steph20180201';

load([dirdata,filedata])%,'data*')
clear Folder* dirdata fileout

dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/Model/';
fileout_txt='ExperimentSummmary.txt';



%{
gene_excl = {'caup','hth','DIAP1','fz3','mirr','exd','sd','wg'};%sd,wg
array_list_gene_excl=[];
for elem=gene_excl
    array_list_gene_excl = [array_list_gene_excl; union({'hh','Delta'},elem)];
end
%}
%
% Delta exclut car pas 12apf ; hh exclut car mauvais ; fz3 : exclut non fourni par St?phane
array_list_gene_excl=[{'hh','Delta','caup','hth','DIAP1','fz3'};{'hh','Delta','mirr','exd','sd','wg'}];
%}
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
%MethGene01=num2str(cutoff_gene);%'ManualImage-Yohanns'%'AutomaticVarianceBased-Stephane'
list_meth_selec_lambda  = {'1SE';'MaxFeature'}';
max_feature=10;
type_model='lasso';

%% Binarize data (rough cutoff selected manually)

%{
MethGene01={num2str(cutoff_gene);'ManualImage-Yohanns';'AutomaticVarianceBased-Stephane'};
Gene_bin = MaskSR_vect;
%}

nsimu=0;
for MethGene01 = {num2str(cutoff_gene),'ManualImage-Yohanns','AutomaticVarianceBased-Stephane'}
    
    
    switch MethGene01{:}
        case 'AutomaticVarianceBased-Stephane'
            Gene_bin = MaskSR_vect;
        case 'ManualImage-Yohanns'
            Gene_bin = MaskYB_vect;
        otherwise
            for g=1:length(list_gene)
                Gene_bin(:,g)=Gene_vect(:,g)>cutoff_gene(g)+0;
            end
    end
    
    %% ------------------------------------------------------------------------
    %% %%%%%%%%%%%%%%%%%%%
    
    for Interact01=0:1
        
        [nl,~] = size(array_list_gene_excl);
        for ix = 1:nl
            list_gene_excl={array_list_gene_excl{ix,:}};
            
            %% EXCLUDE GENE and Apply mask
            [list_gene_incl,ia]=setdiff(list_gene,list_gene_excl);
            
            % Most restrictive mask : only keep point if area ratio 1
            minAR = min(AreaRatio_vect(:,ia),[],2); % min per row
            mask = (minAR==1) .* (Phys_vect(:,3)==1);% 1 if area ratio =1 and y has a value
            
            X = Gene_bin(mask==1,ia)+0;
            y = Phys_vect(mask==1,1)+0;
            
            %% Create interaction terms
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
            % 85% training-validation / 15% test (n=393/163/97)
            myseed =1;
            rand('seed', myseed);
            cv = cvpartition(length(y),'holdout',0.15);
            
            % Test set
            Xtest = X(test(cv),:);
            ytest = y(test(cv),:);
            
            % Training and validation set
            Xtrain = X(training(cv),:);
            ytrain = y(training(cv),:);
            
            %  cv = cvpartition(length(ytrain),'holdout',0.25/0.85);
            
            %% Generate models for various prenalty parameter lambda
            rand('seed', 19479);
            [B,FitInfo] = lasso(Xtrain,ytrain,'CV',10);%,'PredictorNames',list_gene_incl);
            %[ax,figh]=lassoPlot(B,FitInfo,'PlotType','CV');%  green circle ==> |Lambda| with minimum cross-validation error.
            % saveas(figh,[dirout,'Lambda_lasso_Interact',MethGene01,num2str(Interact01),'_',num2str(nsimu)])
            % nsimu=nsimu+1;
            %}
            %% BEST MODEL - weight, R2, nb of variables selected
            for meth =list_meth_selec_lambda
                
                if ismember(meth,'MaxFeature')
                    n_feature_per_model =sum(B~=0);
                    idx_max_feature = find(n_feature_per_model<=max_feature);
                    [~,temp]=min(FitInfo.MSE(:,idx_max_feature));
                    idx_best_lambda = idx_max_feature(temp);
                    
                else
                    idx_best_lambda = FitInfo.(['Index',meth{:}]);
                end
                
                
                
                % variable selected
                WeightOpti = B(:,idx_best_lambda);% best lambda can be 'MinMSE' or '1SE'
                ListVarOpti= list_gene_interact(WeightOpti~=0);
                TableParaOpti=table(list_gene_interact',WeightOpti,'VariableNames',{'name','weight'});
                TableParaOpti(TableParaOpti.weight==0,:)=[];
                n_var_selected = length(ListVarOpti);
                
                % Predicted values
                ytrain_pred = Xtrain*B(:,idx_best_lambda)+FitInfo.Intercept(idx_best_lambda);
                ytest_pred = Xtest*B(:,idx_best_lambda)+FitInfo.Intercept(idx_best_lambda);
                
                % R2 and R2adj
                [R2train,R2adjtrain] = myR2(ytrain,ytrain_pred,n_var_selected);%train set
                [R2test,R2adjtest] = myR2(ytest,ytest_pred,n_var_selected);%test set
                
                median_err_train = nanmedian( ytrain- ytrain_pred );
                median_err_test  = nanmedian( ytest - ytest_pred  )  ;
                
                %stat_err_train = [ min(abs(ytrain-ytrain_pred)),nanmedian((ytrain-ytrain_pred)),max(abs(ytrain-ytrain_pred))];
                %stat_err_test  = [ min(abs(ytest-ytest_pred))  ,nanmedian((ytest-ytest_pred))  ,max(abs(ytest-ytest_pred))];
                pct_err_train          =( ytrain(:)-  ytrain_pred(:))./ytrain(:);
                pct_err_test           =( ytest(:)-  ytest_pred(:))./ytest(:);
                pct_err_train_within10 =length(find(pct_err_train>=-0.10&pct_err_train<=0.10))/length(pct_err_train);
                pct_err_train_within20 =length(find(pct_err_train>=-0.20&pct_err_train<=0.20))/length(pct_err_train);
                pct_err_test_within10  =length(find(pct_err_test>=-0.10&pct_err_test<=0.10))/length(pct_err_test);
                pct_err_test_within20  =length(find(pct_err_test>=-0.20&pct_err_test<=0.20))/length(pct_err_test);
                
                
                %{
    h=figure
    scatter(ytest,ytest_pred);
    hold on
    plot([0,max(ytest)],[0,max(ytest)])
    hold off
    title(['Test set: observed vs predicted (lambda =',meth{:},')'])
    saveas(h,[dirout,'YvsYhat_test_set_MSE'])
                %}
                WriteResultModel(filedata,dirout,fileout_txt,list_gene_interact,list_gene_excl,Interact01,...
                    MethGene01{:},myseed, type_model,meth{:},FitInfo.Lambda(idx_best_lambda),FitInfo.MSE(idx_best_lambda),mat2str([length(ytrain),length(ytest)]),...
                    n_var_selected,TableParaOpti,FitInfo.Intercept(idx_best_lambda),R2train,R2adjtrain,pct_err_train_within10,pct_err_train_within20,R2test,R2adjtest,pct_err_test_within10,pct_err_test_within20);
                %}
                clear h
            end % type of lambda
        end %for ix:list gene
    end % for jx: interaction
end %mask YB SR LN
%% FULL SAMPLE: look for best model (find colorbar)
% comparison original versus predicted
%{
% pred all point
y_pred = NaN(length(Phys_vect),1);
y_pred(mask==1) = X*B(:,idx_best_lambda)+FitInfo.Intercept(idx_best_lambda);

cmin =min(min([y_pred, y_pred,Phys_vect(:,1)]));
cmax =max(max([y_pred, y_pred,Phys_vect(:,1)]));

% reshape for image
y_pred = reshape(y_pred,size(data_phys.div));%training set
%}
%% FULL SAMPLE (train + test sets): graphical inspection
%{
h=figure
subplot(2,2,1)%observed
imagesc(reshape(mask,size(data_phys.div)).*data_phys.div,[cmin,cmax])
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Division (observed)')

subplot(2,2,2)% predicted image
imagesc(y_pred,[cmin,cmax]);
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Division (predicted)')

subplot(2,2,3)% map error
imagesc(data_phys.div-y_pred)
set(gca,'xtick',[]);set(gca,'ytick',[])
title('Error: obs-pred')
colorbar

subplot(2,2,4)% histogram error
histogram(data_phys.div-y_pred)

saveas(h,[dirout,'GraphBestModel_',best_lambda,'.pdf'])
save([dirout,'ModelDivision'])
%}
%% Additional functions

% calcul R2
function [R2,R2adj] = myR2(y,ypred,k)
n = length(y);
idx_nan = isnan(y);
SStot = var(y)*n;
SSres = sum( (y(:)-ypred(:)).^2 );
R2    = 1 - SSres/SStot;
R2adj = 1-(1-R2)*(n-1)/(n-k-1);
end

%% file exchange
function res=x2fx_names(PredictorNames,model)
%function res = x2fx_name(PredictorNames, model)
%   returns a cell array of "compound" predictor variable names that
%   correpond to the x2fx'ed matrix of "normal" predictor variable names.
%

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
