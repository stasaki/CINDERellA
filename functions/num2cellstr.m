function [output]=num2cellstr(gene_id)
output = cell(length(gene_id),1);
for i=1:length(gene_id)
    output{i} = num2str(gene_id(i));
end