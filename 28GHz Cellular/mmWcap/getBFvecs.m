% Calculates the optimal(maximizing cov matrices) Tx/Rx BF vectors.

% freq = 28;
% nantBS1 = 8;
% nantUE1 = 4;
% covFn = sprintf('../chanMod/data/covData%d_%dx%dx%dx%d', ...
%     freq, nantBS1, nantBS1, nantUE1, nantUE1);
% load(covFn);
% Qtx: 64x64x1000   Qrx: 16x16x1000
ncov = size(QrxTot,3); 
Wtx_opt1 = zeros(size(QtxTot,1),ncov);
Wrx_opt1 = zeros(size(QrxTot,1),ncov);

for j=1:ncov
    % nb: eigs, at least in my Matlab version returns normalized eigenvecs.
    % However, when calculating e.g., the Rx BF gains, need to normalize 
    % the size of the *Tx* antenna array! 
    [Wtx_opt1(:,j),dtx] = eigs(QtxTot(:,:,j),1);  % Maximum eigen vector
    [Wrx_opt1(:,j),drx] = eigs(QrxTot(:,:,j),1);       
end

save('data/TxBFVecs.mat','Wtx_opt1');
save('data/RxBFVecs.mat','Wrx_opt1');




