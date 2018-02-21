function WriteResultModel(filedata,dirout,fileout_txt,list_gene,list_gene_excl,Interact01,...
 MethGene01,myseed, type_model,meth,best_lambda,MSE,size_sets,...
 n_var_selected,TableParaOpti,Weight_Intercept,R2train,R2adjtrain,pct_err_train_within10,pct_err_train_within20,pct_err_test_within10,pct_err_test_within20)



fileout = [dirout,fileout_txt];
% create

if exist([fileout],'file')==2
    % if file exists, append to the end
    fileID = fopen(fileout,'a');
    
else % if file does not exist, add header
    fileID = fopen(fileout,'w');
    fprintf(fileID,'Date \t DataPath \t Var_Incl \t Var_Excl \t Interact01 \t MethToBinarize \t seed \t Model \t TypeLambda \t lambda \t size(train,test) \tn_var_select \t Var_select \t Weight \t Intercept \t MSE_train \t R2_train \t R2adj_train \t pct_err_train_within10 \t pct_err_train_within20  \t pct_err_test_within10 \t pct_err_test_within20 \n') ;   
end


temp = list_gene';
fprintf(fileID,'%s\t',datetime('now'));
fprintf(fileID,'%s\t',filedata);
fprintf(fileID,'%s,' ,list_gene{:})     ;fprintf(fileID,'\t');
fprintf(fileID,'%s,' ,list_gene_excl{:});fprintf(fileID,'\t');
fprintf(fileID,'%d\t',Interact01); 
fprintf(fileID,'%s\t' ,MethGene01);
fprintf(fileID,'%f\t',myseed); 
fprintf(fileID,'%s\t',type_model); 
fprintf(fileID,'%s\t',meth); 
fprintf(fileID,'%s\t',best_lambda); 
fprintf(fileID,'%s\t',size_sets); 
fprintf(fileID,'%d\t',n_var_selected);
fprintf(fileID,'%s,' ,TableParaOpti.name{:});    fprintf(fileID,'\t');
fprintf(fileID,'%f,' ,TableParaOpti.weight);    fprintf(fileID,'\t');
fprintf(fileID,'%s\t',Weight_Intercept); 


fprintf(fileID,'%s\t',MSE); 

fprintf(fileID,'%.3f\t',R2train);
fprintf(fileID,'%.3f\t',R2adjtrain);
fprintf(fileID,'%.3f\t',pct_err_train_within10);
fprintf(fileID,'%.3f\t',pct_err_train_within20);
fprintf(fileID,'%.3f\t',pct_err_test_within10);
fprintf(fileID,'%.3f\t',pct_err_test_within20);

fprintf(fileID,'\n');
fclose(fileID);


end

