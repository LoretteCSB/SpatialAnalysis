function WriteResultModel01(filedata,dirout,fileout_txt,list_gene_interact,list_gene_excl,Interact01,...
                    MethGene01,myseed,myseed2, type_model,meth,size_sets,best_lambda,...
                    n_var_selected,ListVarOpti,WeightOpti,cm_valid,cm_test)

                
fileout = [dirout,fileout_txt];
% create

if exist([fileout],'file')==2
    % if file exists, append to the end
    fileID = fopen(fileout,'a');
    
else % if file does not exist, add header
    fileID = fopen(fileout,'w');
    fprintf(fileID,'Date \t DataSource \t Var_Incl \t Var_Excl \t Interact01 \t MethToBinarize \t seed \t seed2 \t Model \t TypeLambda \t lambda (c) \t size(train,valid,test) \t n_var_select \t Var_select \t Weight  \t TN_valid \t FN_valid \t FP_valid \t TP_valid \t Accuracy_valid \t TN_test \t FN_test \t FP_test \t TP_test \tAccuracy_test\n') ;   
end

fprintf(fileID,'%s\t',datetime('now'));
fprintf(fileID,'%s\t',filedata);
fprintf(fileID,'%s,' ,list_gene_interact{:})     ;fprintf(fileID,'\t');
fprintf(fileID,'%s,' ,list_gene_excl{:});fprintf(fileID,'\t');
fprintf(fileID,'%d\t',Interact01); 
fprintf(fileID,'%s\t' ,MethGene01);
fprintf(fileID,'%f\t',myseed); 
fprintf(fileID,'%f\t',myseed2); 
fprintf(fileID,'%s\t',type_model); 
fprintf(fileID,'%s\t',meth); 
fprintf(fileID,'%s\t',best_lambda); 
fprintf(fileID,'%s\t',size_sets); 
fprintf(fileID,'%d\t',n_var_selected);
fprintf(fileID,'%s,' ,ListVarOpti{:});    fprintf(fileID,'\t');
fprintf(fileID,'%f,' ,WeightOpti);    fprintf(fileID,'\t');
fprintf(fileID,'%d\t',cm_valid); 
fprintf(fileID,'%d%%\t',(cm_valid(1,1)+cm_valid(2,2))/sum(cm_valid(:))*100); 

fprintf(fileID,'%d\t',cm_test);
fprintf(fileID,'%d\t',(cm_test(1,1)+cm_test(2,2))/sum(cm_test(:))); 

fprintf(fileID,'\n');
fclose(fileID);


end

