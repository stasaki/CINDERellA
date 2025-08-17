function BIC=score_dag_onfly(dag,expdata)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


BIC = 0;
for i=1:length(dag)
    paset = find(dag(:,i)==1);
    BIC = BIC + CcalcLS_one(expdata,i,paset);
end

end

