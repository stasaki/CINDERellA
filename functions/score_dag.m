function BIC=score_dag(dag,LS,pa_limit,dummy,padim)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


BIC = 0;
for i=1:length(dag)
    paset = find(dag(:,i)==1);
    paset = mysub2ind(padim,[paset',repmat(dummy,1,pa_limit-length(paset'))]);
    BIC = BIC + LS{i}.score(LS{i}.indx==paset);
end

end

