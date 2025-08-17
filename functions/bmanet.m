function [indx,postprob] = bmanet(input_args,LL,net)
%bmanet Summary of this function goes here
% INPUTS:
%	input_args:	sampled graphs
%	LL:     log likelihood 
%   net:   generate consensus network?
% OUTPUTS:
%	indx:		binary edge matirx
%   postprob: post probablity weight
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License

LL = exp(LL-max(LL));
postprob=zeros(size(input_args{1}));
for i=1:length(input_args)
    postprob=postprob+input_args{i}*LL(i);
end


if net ==1
    indx=postprob>postprob';
    indx((sum(LL(1:i))-(postprob+postprob'))>postprob)=0;
else
    indx=postprob>0.5*sum(LL(1:i));
end

end

