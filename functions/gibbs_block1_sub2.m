%% sub code for learn_struct_mcmc_gibbs_block2_con
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License


% select 1 nodes for modification
node1 = nperm(mod(bc1,n)+1);
bc1 = bc1 +1;
pa_node1 = find(dag(:,node1))'; %get current parents for node1

%orphaning node1 and updating path count matrix(cmat)
for i=1:length(pa_node1)
    dag(pa_node1(i),node1)=0;
    cmat = cmat - cmat(:,pa_node1(i))*cmat(node1,:);
end

%retrive descendant of node1
pa_kj1c = find(cmat(node1,:)>0);

%compute parent sets
%find all possible parent combination for node1
n1_indx1=find(~any(LS{node1}.lookup(:,pa_kj1c),2));
%get local score for possible child and parents combination
n1_LS_Spa_kj1_kj2=LS{node1}.score(n1_indx1);


temp2(:)=dummy;
temp2(1:length(pa_node1))=pa_node1;
paset_trt=mysub2ind(padim,temp2);
n1_ls_old = n1_LS_Spa_kj1_kj2(ismembc2(ismembc2(paset_trt, LS{node1}.indx),n1_indx1));


% calc normalizing constant
n1_LS_Spa_kj1_kj2_exp=exp(n1_LS_Spa_kj1_kj2-max(n1_LS_Spa_kj1_kj2));
n1_LS_Spa_kj1_kj2_exp=n1_LS_Spa_kj1_kj2_exp./sum(n1_LS_Spa_kj1_kj2_exp);

%Draw new parent
u=rand;
postPend=0;% store n1_ls for debug use

[n1_loc,~,n1_ls]=drawpa3(n1_LS_Spa_kj1_kj2_exp,n1_LS_Spa_kj1_kj2,postPend,u);
n1_pa=LS{node1}.parent(n1_indx1(n1_loc),:);

%Update dag and cmat
n1_pa=rmzero(n1_pa,dummy);
dag(n1_pa,node1)=1;



for j=1:length(n1_pa)
    cmat = cmat + cmat(:,n1_pa(j))*cmat(node1,:);
end
%Update current LL

LL = LL + n1_ls - n1_ls_old;
