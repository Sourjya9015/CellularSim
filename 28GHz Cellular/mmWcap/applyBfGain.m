function [pathLossDL,pathLossUL,bfGainRxDLdes] = applyBfGain(pathLoss, Icell, bfOpt)

% Get options
intNull = bfOpt.intNull;
covFn = bfOpt.covFn;
noisepow = bfOpt.noisepow;
txpow = bfOpt.txpow;

% Get dimensions
[nbs,nue] = size(pathLoss);

%nb:
bfGainRxDLdes = zeros(nue,1);
%bfGainRxDL2 = zeros(nue,nbs);

% Generate random covariance matrices for each link
% Loads the data from the chan/... data section. contains the Int (interference)
% and Des (desired) beamforming gains. Taken from measurements.
load(covFn);
nantUE = size(QrxTot,1);
nantBS = size(QtxTot,2);
ncov = size(QrxTot,3);
Icov = randi(ncov,nbs,nue);

% Generate random TX gains by selecting a random gain "interfering" gain
% on the interfering links and a random "serving" gain on the serving
% links.
bfGainTxDL = bfGainIntTx(Icov);
bfGainRxDL = bfGainIntRx(Icov);
for iue = 1:nue
    ibs = Icell(iue);
    bfGainTxDL(ibs,iue) = bfGainDesTx(Icov(ibs,iue));
    bfGainRxDL(ibs,iue) = bfGainDesRx(Icov(ibs,iue));
    bfGainRxDLdes(iue,:) = bfGainRxDL(ibs,iue);
end


pathLossUL = pathLoss - bfGainTxDL - bfGainRxDL;
pathLossDL = pathLossUL';


if ~intNull
    return;
end

% Downlink interference nulling
pathLoss2 = pathLoss - bfGainTxDL;

% Compute thermal noise relative to TX power
Pnoise = 10.^(0.1*(noisepow-mean(txpow)));

% Compute the RX gains
% To simplify the computations, we compute the interference covariance
% matrix only on the dominant interferers within a threshold.  Other
% interferes are approximated as white noise.

intThresh = 30;
for iue = 1:nue
    
    % Find the dominant and non-dominant interferers
    ibs = Icell(iue);
    p = (pathLoss2(:,iue) < pathLoss2(ibs,iue) + intThresh);
    Idom = find( p & ((1:nbs)' ~= ibs));
    Inondom = find( ~p & ((1:nbs)' ~= ibs));
    
    % Compute the interference power on the non-dominant interfers
    plint = pathLoss2(Inondom,iue) - bfGainRxDL(Inondom,iue);
    Pint = sum(10.^(-0.1*plint));
    QrxInt = Pint*eye(nantUE);
    
    % Add the interference from the dominant interferers
    nbs1 = length(Idom);
    for ibs1 = 1:nbs1
        ibs2 = Idom(ibs1);
        QrxInt= QrxInt + 10^(-0.1*pathLoss2(ibs2,iue))*QrxTot(:,:,Icov(ibs2,iue))/nantBS;
    end
    
    % Add the thermal noise or a minimum level to keep inverse well-conditioned
    QrxInt = QrxInt + max(trace(QrxInt)/nantUE*1e-4,Pnoise)*eye(nantUE);
    
    % Get the desired RX covariance matrix
    QrxDes = QrxTot(:,:,Icov(ibs,iue));
    
    % Compute the optimal BF vector    
    QrxInt2 = inv(sqrtm(QrxInt));
    Q = QrxInt2'*QrxDes*QrxInt2;
    [wrx,d] = eigs(Q,1);
    wrx = wrx/norm(wrx);
    
    % Compute the BF gain on the desired signal
    bfGainRxDL(ibs,iue) = 10*log10( real(wrx'*QrxDes*wrx/nantBS));

    % Compute the BF gain on the dominant interfering links
    for ibs1 = 1:nbs1
        ibs2 = Idom(ibs1);
        icov = Icov(ibs2,iue);
        bfGainRxDL(ibs,iue) = 10*log10( real(wrx'*QrxTot(:,:,icov)*wrx/nantBS));
    end
    
end
pathLossDL = (pathLoss2 - bfGainRxDL)';


