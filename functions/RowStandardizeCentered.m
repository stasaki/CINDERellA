function [ data ] = RowStandardizeCentered( data )
% sub-function from Matrix eQTL
% Matrix eQTL by Andrey A. Shabalin
% http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
% Build 1.2.0

                div = sqrt(sum(data.^2,2));
				div(div==0) = 1;
				data = bsxfun(@rdivide, data, div);
%                 data = data./(div*ones(size(data,2),1)');
%                 data = data./repmatC(div,1,size(data,2));
end

