%% sub code for learn_struct_mcmc_gibbs_block12_v2
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License


pa_node1 = find(dag(:,node(1)))'; %get current parents for node1
pa_node2 = find(dag(:,node(2)))'; %get current parents for node2
pa_node3 = find(dag(:,node(3)))'; %get current parents for node2

%%orphaning node1 and node2 and updating path count matrix(cmat)
for i=1:length(pa_node1)
    dag(pa_node1(i),node(1))=0;
    cmat = cmat - cmat(:,pa_node1(i))*cmat(node(1),:);
end
for i=1:length(pa_node2)
    dag(pa_node2(i),node(2))=0;
    cmat = cmat - cmat(:,pa_node2(i))*cmat(node(2),:);
end
for i=1:length(pa_node3)
    dag(pa_node3(i),node(3))=0;
    cmat = cmat - cmat(:,pa_node3(i))*cmat(node(3),:);
end

%%
n1_lookup=LS{node(1)}.lookup;
n2_lookup=LS{node(2)}.lookup;
n3_lookup=LS{node(3)}.lookup;


%retrive parent candidates
pa_kj1c = find(cmat(node(1),:)~=0); %descendant for node1
pa_kj2c = find(cmat(node(2),:)~=0); %descendant for node2
pa_kj3c = find(cmat(node(3),:)~=0); %descendant for node2

%union descendant for node1 and node2 and node3
pa_kj1c_kj2c_kj3c=unionKPM(pa_kj3c,unionKPM(pa_kj1c,pa_kj2c));%back here later

%%compute parent sets
%% find possible parent candidates TYPE0 for node1, node2 and node3
n1_3indx{1}=find(~any(n1_lookup(:,pa_kj1c_kj2c_kj3c),2));
n2_3indx{1}=find(~any(n2_lookup(:,pa_kj1c_kj2c_kj3c),2));
n3_3indx{1}=find(~any(n3_lookup(:,pa_kj1c_kj2c_kj3c),2));

%% find possible parent candidates TYPE1 for node1
%node1<-node2 
n1_3indx{2}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(~any(n1_lookup(:,pa_kj3c),2)));
%node1<-node3 
n1_3indx{3}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj3c),2)&(~any(n1_lookup(:,pa_kj2c),2)));

%find possible parent candidates TYPE1 for node2
%node2<-node1
n2_3indx{2}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(~any(n2_lookup(:,pa_kj3c),2)));
%node2<-node3
n2_3indx{3}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj3c),2)&(~any(n2_lookup(:,pa_kj1c),2)));


%find possible parent candidates TYPE1 for node3
%node3<-node1
n3_3indx{2}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(~any(n3_lookup(:,pa_kj2c),2)));
%node3<-node2
n3_3indx{3}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj2c),2)&(~any(n3_lookup(:,pa_kj1c),2)));


%% find possible parent candidates TYPE2 for node1
%node1<-(node2,node3)
n1_3indx{4}=find((~any(n1_lookup(:,pa_kj1c),2))&any(n1_lookup(:,pa_kj2c),2)&(any(n1_lookup(:,pa_kj3c),2)));


%find possible parent candidates TYPE2 for node2
%node2<-(node1,node3)
n2_3indx{4}=find((~any(n2_lookup(:,pa_kj2c),2))&any(n2_lookup(:,pa_kj1c),2)&(any(n2_lookup(:,pa_kj3c),2)));


%find possible parent candidates TYPE2 for node3
%node3<-(node1,node2)
n3_3indx{4}=find((~any(n3_lookup(:,pa_kj3c),2))&any(n3_lookup(:,pa_kj1c),2)&(any(n3_lookup(:,pa_kj2c),2)));



%% get local score for possible child and parents combination for node1

sctemp=LS{node(1)}.score;
for i=1:4
n1_3LS_indx{i}=sctemp(n1_3indx{i});
end

temp2(:)=dummy;
temp2(1:length(pa_node1))=pa_node1;
paset_trt=mysub2ind(padim,temp2);
n1_ls_old = sctemp(ismembc2(paset_trt, LS{node(1)}.indx));


%% get local score for possible child and parents combination for node2

sctemp=LS{node(2)}.score;
for i=1:4
n2_3LS_indx{i}=sctemp(n2_3indx{i});
end


temp2(:)=dummy;
temp2(1:length(pa_node2))=pa_node2;
paset_trt=mysub2ind(padim,temp2);
n2_ls_old = sctemp(ismembc2(paset_trt, LS{node(2)}.indx));


%% get local score for possible child and parents combination for node3

sctemp=LS{node(3)}.score;
for i=1:4
n3_3LS_indx{i}=sctemp(n3_3indx{i});
end



temp2(:)=dummy;
temp2(1:length(pa_node3))=pa_node3;
paset_trt=mysub2ind(padim,temp2);
n3_ls_old = sctemp(ismembc2(paset_trt, LS{node(3)}.indx));


%% calc normalizing constant


for i=1:12
    pb_3n2{i}=n2_3LS_indx{mar3code2(i,1)}+logsumexp(n3_3LS_indx{mar3code2(i,2)},1);
end


for i=1:25
    pb_3n1{i}=n1_3LS_indx{mar3code1(i,1)}+logsumexp(pb_3n2{mar3code1(i,2)},1);
end

normconst=logsumexp(cat(1,pb_3n1{:}),1);



%Draw new parent
u=rand;
postPend=0;% store n1_ls and n2_ls for debug use

% draw parant for n1
for i=1:25
    [n1_loc,postPend,n1_ls]=drawpa3(exp(pb_3n1{i}-normconst),n1_3LS_indx{mar3code1(i,1)},postPend,u);
    
    if ~isempty(n1_loc)
        sp_loc=i;
        n1_pa = LS{node(1)}.parent(n1_3indx{mar3code1(i,1)}(n1_loc),:);
        break;
    end
end
[n2_loc,postPend,n2_ls]=drawpa3(exp(n1_ls+pb_3n2{mar3code(sp_loc,2)}-normconst),n2_3LS_indx{net3code(sp_loc,2)},postPend,u);
[n3_loc,~,n3_ls]=drawpa3(exp(n1_ls+n2_ls+n3_3LS_indx{mar3code(sp_loc,3)}-normconst),n3_3LS_indx{net3code(sp_loc,3)},postPend,u);
n2_pa = LS{node(2)}.parent(n2_3indx{net3code(sp_loc,2)}(n2_loc),:);
n3_pa = LS{node(3)}.parent(n3_3indx{net3code(sp_loc,3)}(n3_loc),:);


%Update dag and cmat
n1_pa=rmzero(n1_pa,dummy);
n2_pa=rmzero(n2_pa,dummy);
n3_pa=rmzero(n3_pa,dummy);
dag(n1_pa,node(1))=1;
dag(n2_pa,node(2))=1;
dag(n3_pa,node(3))=1;


for j=1:length(n1_pa)
    cmat = cmat + cmat(:,n1_pa(j))*cmat(node(1),:);
end
for j=1:length(n2_pa)
    cmat = cmat + cmat(:,n2_pa(j))*cmat(node(2),:);
end
for j=1:length(n3_pa)
    cmat = cmat + cmat(:,n3_pa(j))*cmat(node(3),:);
end


%Update current LL
    
LL = LL + n1_ls + n2_ls + n3_ls - n1_ls_old - n2_ls_old - n3_ls_old;
