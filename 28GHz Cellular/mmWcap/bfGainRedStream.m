function [pathLossDL1, pathLossUL1] = bfGainRedStream( pathLossDL, pathLossUL, Icell, bfOpt)
% bfGainSing:  Adjusts the path loss due to restriction to a finite number
% of streams.
%
% Assumes a 2D ULA with lambda/2 antenna spacing.

% Get parameters
nstream = bfOpt.nstream;                                 % num streams
antSp = 0.5;                                                      % normalized antenna spacing

% Get dimensions
[nue,nbs] = size(pathLossDL);
nantH = round( sqrt(bfOpt.nantBS) );                     % num horiz antennas at BS

% Compute the SINR with the max BF gain
p = 10.^(0.1*(repmat(bfOpt.txpow',nue,1) - pathLossDL));  % 3900x390
totPow = sum(p,2);                                                              % 3900x1
sigPowMax = zeros(nue,1);

for iue = 1:nue
    sigPowMax(iue) = p(iue, Icell(iue));                 % power received from the connected BS
end
noisePow = 10^(0.1*bfOpt.noisepow);
intPow = totPow - sigPowMax + noisePow;
% BeamFormed SINR, in capsim.m the input arguments are the PL with BF gain
dlSinrMax = sigPowMax ./ intPow;                        

if 0
    dlSinrMaxdB = 10*log10( dlSinrMax );
    plot(sort(dlSinrMaxdB), (1:nue)/nue);
    grid on;
end

% Loop over serving BSs
gainLoss = zeros(nue,1);
for ibs = 1:nbs
    
    % Find UEs connected to the BS
    Iue1 = find( Icell == ibs );
    nue1 = length(Iue1);
    if (nue1 == 0)
         continue;
    end
    
    % Generate random directions
    costh = cos( pi*rand(1,nue1) );                                          % random nue1 directions
    V = exp(1i*2*pi*antSp*(0:nantH-1)'*costh) / sqrt(nantH); % random Az. directions w.r.t to the BS
    sinrMax = dlSinrMax(Iue1);
    
    % Find the optimal direction
    param.nit = 100;
    param.verbose = false;
    param.nstream = nstream;
    [U, sinr1] = bfSingOpt( V, sinrMax, param );
    gainLoss(Iue1) = 10*log10( sinrMax ./ sinr1 );
    
    
end

% Adjust the path loss
pathLossDL1 = pathLossDL;
pathLossUL1 = pathLossUL;
for iue = 1:nue
    pathLossDL1(iue, Icell(iue)) = pathLossDL1(iue, Icell(iue)) + gainLoss(iue);
    pathLossUL1(Icell(iue),iue) = pathLossUL1(Icell(iue),iue) + gainLoss(iue);
end
end


function [U, sinr] = bfSingOpt( V, sinrMax, param )
% Solves the maximization
%   max_u \sum_i log( rate( (u'*V(:,i))^2*sinrMax ) )

% Get parameters
verbose = param.verbose;
nstream = param.nstream;
nit = param.nit;
alpha = 0.5;

% Initialize to maximal singular vector
[U0,S0,V0] = svd(V);
U = U0(:,1:nstream); % choose largest nstream directions
gain = V'*U;
sinr = sinrMax.*sum(abs(gain).^2,2);
[util, grads] = fneval(sinr);
grad = 2*V* (gain.*repmat(sinrMax.*grads, 1, nstream));
step = 1e-4;
for it = 1:nit
    
    % Try candidate
    U1 = U + step*grad;
    gain1 = abs(V'*U1);
    sinr1 = sinrMax.*sum(gain1.^2,2);
    util1 = fneval(sinr1);
    
    % Check if pass
    dutilEst = sum(sum(real(conj(grad).*(U1-U))));
    pass = (util1 > util + alpha*dutilEst);
    if (pass)
        step = 2*step;
        U = U1 / sqrtm(U1'*U1);
        gain = V'*U;
        sinr = sinrMax.*sum(abs(gain).^2,2);
        [util, grads] = fneval(sinr);
        grad = 2*V* (gain.*repmat(sinrMax.*grads, 1, nstream));        
    else
        step = 0.5*step;
    end
    if (verbose)
        fprintf(1, 'it=%d util=%12.4e step=%12.4e\n', it, util, step);
    end
end


end


function [util, grads, se ] = fneval( sinr )
% Compute the objective function and gradient
% 
% The objective function is util = sum( log( rate )),
% where se = min( log2(1 + sinrLoss*sinr), maxse)
sinrLossdB = 3;
maxse = 4.8;
sinrLoss = 10^(-0.1*sinrLossdB);
se = log2(1+sinrLoss*sinr);
I = (se < maxse);
se = min(se, maxse);
util = sum( log(se) );

% Compute the gradient
if (nargout >= 2)
    dutilrate = 1./se;
    dratesinr = sinrLoss/log(2)./(1+sinrLoss*sinr) .* I;
    grads = dutilrate.*dratesinr;
end
end



