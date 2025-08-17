function [ output_args1 ] = countedgefreq( input_args)
%countedgefreq Summary of this function goes here
% INPUTS:
%	input_args:	sampled graphs
% OUTPUTS:
%	output_args1:		edge frequency matrix
%
%   Author: Shinya Tasaki, Ph.D. (stasaki@gmail.com)
%   Code covered by the 3-clause BSD License

output_args1 = input_args{1};
for i=2:length(input_args)
    if isempty(input_args{i})
        break;
    else
    output_args1 = output_args1+input_args{i};
    end
end

output_args1 = output_args1./i;

end

