function [res] = getCosineSmilarityOfRowVectors(data) 
FZ = data*data'; 
FM = sqrt(sum(   data.^2   ,2));   
res = FZ./(FM+eps)./(FM' +eps) ;
% res( isnan( res )))=0;
% function [res] = cosineSimilarityNaive(data)
% get the dimensions
% [n_row n_col] = size(data);
% % calculate the norm for each row
% %
% norm_r = sqrt(sum(abs(data).^2,2));
% %
% for i = 1:n_row
%     % 
%     for j = i:n_row
%         %
%         res(i,j) = dot(data(i,:), data(j,:)) / (norm_r(i) * norm_r(j));
%         res(j,i) = res(i,j);
%     end
% end
% end
