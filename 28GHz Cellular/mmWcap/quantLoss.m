% quantLoss:  Estimate the loss due to quantization loss
% ======================================================

WtotMHz = 1e3;      % Total system bandwidth
rateTgtMbps = 10;   % Rate target in Mbps
bwLoss = 0.5*0.8;   % Reduction factor due to half duplex constraint
Wsig = 1;           % sub-signal BW in MHz
Nsig = 4;           % num sub-signals per slot
nbit = 5;           % num bits in the quantizer
nrx = 16;    % antenna sizes
ntx = 64;

% Find SNR target:  SNR = P/(N0*Wtot)
gainMax = nrx*ntx;
snrTgt = (2^(rateTgtMbps/bwLoss/WtotMHz)-1)/gainMax;

% Find PSS SNR = P/(N0*Wsig*Nsig)
snrPSS = snrTgt*WtotMHz/(Wsig*Nsig);

% Find inverse coding gain
[del, alpha, xmean] = unifQuant(nbit);

% Find SNR after quantization
snrPSSQ = (1-alpha)*snrPSS/(1 + alpha*snrPSS);
snrTgtQ = (1-alpha)*snrTgt/(1 + alpha*snrTgt); 

% Decrement
loss = 10*log10( snrPSS / snrPSSQ );
lossTgt = 10*log10( snrTgt / snrTgtQ );

