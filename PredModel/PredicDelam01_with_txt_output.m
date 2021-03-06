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
fileout_txt='ExperimentSummaryDelam01.txt';
fileout_txt='Comparison_Nrmlz_Gene_Delam01_20180306.txt';

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

%list_meth_selec_lambda  = {'1SE';'MaxFeature'}';% here I only use
%MaxFeature
meth='MaxFeature';% for delamination map, do  
type_model={'L1_lin_svm'}%'logistic'};

%%

% ---------- Gene Activated or not

for notum = {'half_notum','full_notum'}
    for iseed=1:3
        myseed = iseed
        for max_feature=5:5:10
            for MethGene01 = {'ContinuousExprMinMax','ContinuousExpr1Pct','ContinuousExpr2Pct','ContinuousExpr5Pct'}%'Manual','Automatic',
                
                
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
                        Gene_bin = MaskYB_vect;% default continuous
                end
                %
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
                        
                        y = double(Phys_vect(mask==1,2)>0);% column for apoptosis / code 0 if no apoptosis 1 otherwise
                        
                        % encodage des classes sous la forme -1 et +1
                        classcode = unique(y);
                        ysvm = y;
                        ysvm(y==classcode(1))=-1;
                        ysvm(y==classcode(2))=1;
                        
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
                        cv = cvpartition(length(y),'holdout',0.15);
                        
                        % Test set
                        Xtest = X(test(cv),:);
                        ytest = ysvm(test(cv),:);
                        
                        % temp = Training and validation set
                        Xtemp = X(training(cv),:);
                        ytemp = ysvm(training(cv),:);
                        
                        % now split temp into training and validation
                        myseed2 =mod(100998*iseed+1,151);
                        rand('seed', myseed2);
                        cv2 = cvpartition(length(ytemp),'holdout',0.35);
                        
                        Xvalid = Xtemp(test(cv2),:);
                        yvalid = ytemp(test(cv2),:);
                        
                        % temp = Training and validation set
                        Xtrain = Xtemp(training(cv2),:);
                        ytrain = ytemp(training(cv2),:);
                        
                        clear *temp
                        
                        %% Generate models for various prenalty parameter lambda
                        for kernel =5:6
                            if kernel==5
                                type_model={'L1_lin_svm'};%L1-regul L2-loss SVM
                            else
                                type_model={'L1_logistic'};%L1 regul logistic
                            end
                            
                            c=logspace(-6,6,13);
                            
                            table_synthese=NaN(length(c),8);
                            
                            for ic=1:length(c)
                                myoption=sprintf('-s %g -c %g -B %g -q',kernel,c(ic),1);% s = model type, c = penalty parameter, b = add a bias term w0, -q quiet mode
                                % fit model on train dataset
                                model = train(ytrain, sparse(Xtrain), myoption);
                                % use validation to choose best lambda
                                [yvalid_pred, accuracy, ~] = predict(yvalid, sparse(Xvalid), model);
                                % what do I want to maximize? Accuracy?
                                cm = confusionmatrix(yvalid,yvalid_pred);
                                table_synthese(ic,:)=[kernel,c(ic),length(find(model.w)),accuracy(1),cm(1,1),cm(1,2),cm(2,1),cm(2,2)];
                            end
                            clear cm model
                            
                            
                            table_synthese=array2table(table_synthese);
                            table_synthese.Properties.VariableNames={'kernel','c','nFeatures','accuracy','TN','FP','FN','TP'};
                           
                            % BEST MODEL satisfies 2 criteria :
                            % choose lambda such as we get the best accuracy AND a max number of features selected             
                            [maxacc,~]=max(table_synthese.accuracy(table_synthese.nFeatures<=max_feature));
                            imax=find(table_synthese.accuracy==maxacc);
                            imax=imax(1); % in case several models with same results
                            myoption=sprintf('-s %g -c %g',kernel,table_synthese.c(imax));
                            
                            % re-run best model
                            model = train(ytrain, sparse(Xtrain), myoption);
                            
                            % variable selected
                            WeightOpti = model.w;% for apoptosis we only keep MaxFeature / 1SE not implemented
                            ListVarOpti= list_gene_interact(WeightOpti~=0);
                            TableParaOpti=table(list_gene_interact',WeightOpti','VariableNames',{'name','weight'});
                            TableParaOpti(TableParaOpti.weight==0,:)=[];
                            n_var_selected = length(ListVarOpti);
                            
                            % Predicted values
                            [yvalid_pred, acc_valid, ~] = predict(yvalid, sparse(Xvalid), model);
                            [ytest_pred , acc_test , ~] = predict(ytest,  sparse(Xtest),  model);
                            
                            % Quality model
                            cm_valid = confusionmatrix(yvalid,yvalid_pred);
                            cm_test = confusionmatrix(ytest,ytest_pred);
                            
                            
                            WriteResultModel01(filedata,dirout,notum{:},fileout_txt,list_gene_interact,list_gene_excl,Interact01,...
                                MethGene01{:},myseed,myseed2, type_model{:},[meth,num2str(max_feature)],mat2str([length(ytrain),length(yvalid),length(ytest)]),table_synthese.c(imax),...
                                n_var_selected,TableParaOpti.name,TableParaOpti.weight,cm_valid,cm_test);
                            
                            clear h
                        end %
                    end
                end
                
            end %mask YB SR LN
        end
    end
    
end



%%*******************************************************************%
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
