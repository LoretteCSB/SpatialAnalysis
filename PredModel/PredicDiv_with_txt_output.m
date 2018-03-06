clear
close all
dirdata='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/data/preprocessed/';
filedata ='GeneAndPhys_Grid_Steph20180201';

load([dirdata,filedata])
clear Folder* dirdata fileout

dirout='/Users/Lorette/Google Drive/YohannsBellaiche/Spatial analysis/results/PredicPhysio/Model/';
fileout_txt='Comparison_Nrmlz_Gene_Division_20180306.txt';

%{
gene_excl = {'','caup','hth','DIAP1','fz3','mirr','exd','sd','wg'};%sd,wg
array_list_gene_excl=[];
for elem=gene_excl
    array_list_gene_excl = [array_list_gene_excl; union({'hh','Delta'},elem)];
end
%}
%
% Delta exclut car pas 12apf ; hh exclut car mauvais ; fz3 : exclut non fourni par St?phane
array_list_gene_excl=[{'hh','Delta','','','',''};{'hh','Delta','caup','hth','DIAP1','fz3'};{'hh','Delta','mirr','exd','sd','wg'}];
%}
cutoff_gene = [0.035	0.02	0.04	0.075	0.027	0.15	0.05	0.04	0.02	0.07	0.025	0.1	0.07	0.004	0.07	0.16	0.03	0.025	0.055	0.06	0.02	0.1	0.05];
list_meth_selec_lambda  = {'1SE';'MaxFeature10'}';%'MaxFeature5';
max_feature=5;
type_model='lasso';

nsimu=0;
for notum = {'full_notum','half_notum'}
    for iseed=1:3
        myseed = iseed;
        %for max_feature=5:5:10
        
        % ---------- Gene Activated or not
        for MethGene01 = {'Manual_ContinuousExprMinMax','ContinuousExprMinMax','ContinuousExpr1Pct','ContinuousExpr2Pct','ContinuousExpr5Pct'}%'Manual','Automatic'}
            
            switch MethGene01{:}
                case 'Automatic'
                    Gene_bin = MaskSR_vect;
                case 'Manual'
                    Gene_bin = MaskYB_vect;
                case 'Manual_ContinuousExprMinMax'
                    Gene_bin = Nrmlz(Gene_vect,0).*MaskYB_vect;
                case 'ContinuousExprMinMax'
                    Gene_bin = Nrmlz(Gene_vect,0);
                case 'ContinuousExpr1Pct'
                    Gene_bin = Nrmlz(Gene_vect,0.5);
                case 'ContinuousExpr2Pct'
                    Gene_bin = Nrmlz(Gene_vect,1);
                case 'ContinuousExpr5Pct'
                    Gene_bin = Nrmlz(Gene_vect,2.5);
                otherwise
                    for g=1:length(list_gene)
                        Gene_bin(:,g)=Gene_vect(:,g)>cutoff_gene(g)+0;
                    end
            end
            
            % ----------- Include interactions between pair of genes
            for Interact01=0:1
                
                [nl,~] = size(array_list_gene_excl);
                
                %---------------- change list of variables to exclude
                for ix = 1:nl
                    list_gene_excl={array_list_gene_excl{ix,:}};
                    
                    % EXCLUDE GENE and Apply mask
                    [list_gene_incl,ia]=setdiff(list_gene,list_gene_excl);
                    
                    % Most restrictive mask : only keep point if area ratio 1
                    minAR = min(AreaRatio_vect(:,ia),[],2); % min per row
                    mask = (minAR==1) .* (Phys_vect(:,3)==1);% 1 if area ratio =1 and y has a value
                    
                    if ismember(notum,'half_notum')
                        mask_half_notum = zeros(size(data_phys.ar));
                        mask_half_notum(21:end,:)=1;
                        mask = mask.* mask_half_notum(:);
                    end
                    
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
                    rand('seed', myseed);
                    cv = cvpartition(length(y),'holdout',0.20);
                    
                    % Test set
                    Xtest = X(test(cv),:);
                    ytest = y(test(cv),:);
                    
                    % temp = Training and validation set
                    Xtemp = X(training(cv),:);
                    ytemp = y(training(cv),:);
                    
                    % now split temp into training and validation
                    myseed2 = mod(100998*iseed+1,151);
                    rand('seed', myseed2);
                    cv2 = cvpartition(length(ytemp),'holdout',0.25);
                    
                    Xvalid = Xtemp(test(cv2),:);
                    yvalid = ytemp(test(cv2),:);
                    
                    % temp = Training and validation set
                    Xtrain = Xtemp(training(cv2),:);
                    ytrain = ytemp(training(cv2),:);
                    
                    clear *temp
                    
                    
                    %% Generate models for various prenalty parameter lambda
                    rand('seed', myseed2);
                    
                    % lasso from Matlab, does not work like liblinear or libvm
                    [B,FitInfo] = lasso(Xtrain,ytrain,'CV',3);%,'PredictorNames',list_gene_incl);
                    %[ax,figh]=lassoPlot(B,FitInfo,'PlotType','CV');%  green circle ==> |Lambda| with minimum cross-validation error.
                    
                    % BEST MODEL - weight, R2, nb of variables selected
                    for meth =list_meth_selec_lambda
                        
                        switch meth{:}
                            case 'MaxFeature5' % best model = min MSE AND max feature
                                max_feature=5;
                                n_feature_per_model =sum(B~=0);
                                idx_max_feature = find(n_feature_per_model<=max_feature);
                                [~,temp]=min(FitInfo.MSE(:,idx_max_feature));
                                idx_best_lambda = idx_max_feature(temp);
                                %                                meth_select_lambda = {['MaxFeature_',num2str(max_feature)]
                            case 'MaxFeature10'
                                max_feature=10;
                                n_feature_per_model =sum(B~=0);
                                idx_max_feature = find(n_feature_per_model<=max_feature);
                                [~,temp]=min(FitInfo.MSE(:,idx_max_feature));
                                idx_best_lambda = idx_max_feature(temp);
                            otherwise
                                idx_best_lambda = FitInfo.(['Index',meth{:}]);
                        end
                        meth_select_lambda = meth;
                        
                        % variable selected
                        WeightOpti = B(:,idx_best_lambda);% best lambda can be 'MinMSE' or '1SE'
                        ListVarOpti= list_gene_interact(WeightOpti~=0);
                        TableParaOpti=table(list_gene_interact',WeightOpti,'VariableNames',{'name','weight'});
                        TableParaOpti(TableParaOpti.weight==0,:)=[];
                        n_var_selected = length(ListVarOpti);
                        
                        % Predicted values
                        ytrain_pred = Xtrain*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);
                        yvalid_pred = Xvalid*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);
                        ytest_pred  =  Xtest*B(:,idx_best_lambda) + FitInfo.Intercept(idx_best_lambda);
                        
                        % R2 and R2adj (can only be assess on train)
                        [R2train,R2adjtrain] = myR2(ytrain,ytrain_pred,n_var_selected);%train set
                        
                        % Error
                        median_err_valid = nanmedian( yvalid- yvalid_pred );
                        median_err_test  = nanmedian( ytest - ytest_pred  )  ;
                        
                        pct_err_valid          =( yvalid(:)-  yvalid_pred(:))./yvalid(:);
                        pct_err_valid_within10 =length(find(pct_err_valid>=-0.10&pct_err_valid<=0.10))/length(pct_err_valid);
                        pct_err_valid_within20 =length(find(pct_err_valid>=-0.20&pct_err_valid<=0.20))/length(pct_err_valid);
                        pct_err_test           =( ytest(:)-  ytest_pred(:)  )./ytest(:);
                        pct_err_test_within10  =length(find(pct_err_test>=-0.10&pct_err_test<=0.10))/length(pct_err_test);
                        pct_err_test_within20  =length(find(pct_err_test>=-0.20&pct_err_test<=0.20))/length(pct_err_test);
                        
                        
                        WriteResultModel(filedata,dirout,notum{:},fileout_txt,list_gene_interact,list_gene_excl,Interact01,...
                            MethGene01{:},myseed,myseed2, type_model,meth_select_lambda{:},FitInfo.Lambda(idx_best_lambda),FitInfo.MSE(idx_best_lambda),mat2str([length(ytrain),length(ytest)]),...
                            n_var_selected,TableParaOpti,FitInfo.Intercept(idx_best_lambda),R2train,R2adjtrain,pct_err_valid_within10,pct_err_valid_within20,pct_err_valid_within10,pct_err_test_within20);
                        %}
                        clear h
                    end % type of lambda
                end %for ix:list gene
            end % for jx: interaction
        end %mask YB SR LN
    end
    %  end
end


%%  ****************** Additional functions **************************
function [Gene_vect_nrmlzd,limits_nrmlz] = Nrmlz(Gene_vect,pct)

[nr,nc]=size(Gene_vect);
Gene_vect_nrmlzd= NaN(nr,nc);
limits_nrmlz    = NaN(2,nc);

for col=1:nc
    low=prctile(Gene_vect(:,col),pct);
    high=prctile(Gene_vect(:,col),100-pct);
    limits_nrmlz(1,col)=low;
    limits_nrmls(2,col)=high;
    
    Gene_vect_nrmlzd(:,col) =(Gene_vect(:,col) - low)/(high-low);
    Gene_vect_nrmlzd(Gene_vect(:,col)<low,col) =0;
    Gene_vect_nrmlzd(Gene_vect(:,col)>high,col)=1;
    %figure
    %subplot(1,2,1)
    %hist(Gene_vect(:,col),30)
    %subplot(1,2,2)
    %hist(Gene_vect_nrmlzd(:,col),30)
end


end
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
