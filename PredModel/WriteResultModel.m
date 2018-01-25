function WriteResultModel(dirdata,dirout,fileout_txt,list_gene,list_gene_excl,Interact01,...
 cutoff_gene,myseed, type_model,best_lambda,MSE,size_sets,...
 n_var_selected,TableParaOpti,R2train,R2adjtrain,pct_err_train_within10,pct_err_train_within20,R2test,R2adjtest,pct_err_test_within10,pct_err_test_within20);



% min error, max error, median
% same thing for test set



fileout = [dirout,fileout_txt];
% create

if exist([fileout],'file')==2
    % if file exists, append to the end
    fileID = fopen(fileout,'a');
    
else % if file does not exist, add header
    fileID = fopen(fileout,'w');
    fprintf(fileID,'DataPath \t Var_Incl \t Var_Excl \t Interact01 \t cutoff \t seed \t Model \t lambda \t size(train,test) \tn_var_select \t Var_select \t Weight \t MSE_train \t R2_train \t R2adj_train \t pct_err_train_within10 \t pct_err_train_within20  \t R2_test \t R2adj_test \t pct_err_test_within10 \t pct_err_test_within20 \n') ;   
    %fprintf(fileID,'DataPath;Var_Incl.;Var_Excl.;Interaction01;threshold;seed;Model;lambda;Model_description;MSE;SizeTest;Weight\n') ;   
end

% write relevant output
%fprintf(fileID,'%s %s %s\n',dirdata,list_gene,list_gene_excl);
%fprintf(fileID,'%s \t %s\t %s,%f,%s,%f \n',dirdata,list_gene{:},list_gene_excl{:},Interact01,cutoff_gene,myseed);
  

temp = list_gene';
%temps = temp{:};
fprintf(fileID,'%s\t',dirdata);
fprintf(fileID,'%s,' ,list_gene{:})     ;fprintf(fileID,'\t');
fprintf(fileID,'%s,' ,list_gene_excl{:});fprintf(fileID,'\t');
fprintf(fileID,'%d\t',Interact01); 
fprintf(fileID,'%.4f,' ,cutoff_gene);   fprintf(fileID,'\t')
fprintf(fileID,'%f\t',myseed); 
fprintf(fileID,'%s\t',type_model); 
fprintf(fileID,'%s\t',best_lambda); 
fprintf(fileID,'%s\t',size_sets); 
fprintf(fileID,'%d\t',n_var_selected);
fprintf(fileID,'%s,' ,TableParaOpti.name{:});    fprintf(fileID,'\t');
fprintf(fileID,'%f,' ,TableParaOpti.weight);    fprintf(fileID,'\t');
fprintf(fileID,'%s\t',MSE); 

fprintf(fileID,'%.3f\t',R2train);
fprintf(fileID,'%.3f\t',R2adjtrain);
fprintf(fileID,'%.3f\t',pct_err_train_within10);
fprintf(fileID,'%.3f\t',pct_err_train_within20);
fprintf(fileID,'%.3f\t',R2adjtest);
fprintf(fileID,'%.3f\t',R2test);
fprintf(fileID,'%.3f\t',pct_err_test_within10);
fprintf(fileID,'%.3f\t',pct_err_test_within20);

fprintf(fileID,'\n');
fclose(fileID);


end

