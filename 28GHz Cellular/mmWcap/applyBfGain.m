function [pathLossDL,pathLossUL,bfGainRxDLdes,IntraCellIntf] = applyBfGain(pathLoss, Icell, bfOpt)

% Get options
intNull = bfOpt.intNull;
covFn = bfOpt.covFn;
noisepow = bfOpt.noisepow;
txpow = bfOpt.txpow;

IntraCellIntf = containers.Map(); % define here, set to 0 if intf nulling is not enabled

% Get dimensions
[nbs,nue] = size(pathLoss);

%nb:
bfGainRxDLdes = zeros(nue,1);
%bfGainRxDL2 = zeros(nue,nbs);

% Generate random covariance matrices for each link
% Loads the data from the chan/... data section. contains the Int (interference)
% and Des (desired) beamforming gains. Taken from measurements.
load(covFn);
nantUE = size(QrxTot,1); % QrxTot -> Rx covariance matrix 16x16x1000
nantBS = size(QtxTot,2); % QtxTot -> Tx covariance matrix 64x64x1000
ncov = size(QrxTot,3);
Icov = randi(ncov,nbs,nue);

bfICIReject = 0; % in dB, How much should this be if any?

% Generate random TX gains by selecting a random gain "interfering" gain
% on the interfering links and a random "serving" gain on the serving
% links.
bfGainTxDL = bfGainIntTx(Icov); % Nbs x Nue
bfGainRxDL = bfGainIntRx(Icov); % Nbs x Nue
for iue = 1:nue
    ibs = Icell(iue);
    bfGainTxDL(ibs,iue) = bfGainDesTx(Icov(ibs,iue));
    bfGainRxDL(ibs,iue) = bfGainDesRx(Icov(ibs,iue));
    bfGainRxDLdes(iue,:) = bfGainRxDL(ibs,iue);
end


pathLossUL = pathLoss - bfGainTxDL - bfGainRxDL;
pathLossDL = pathLossUL';




% The code stops here is intf nulling is not enabled.
% i.e all path losses are reduced by Tx+Rx BF gain.
% No sense of interference cancellation due to directivity.
if ~intNull
    
    bfGainIntraTx = bfGainIntTx(Icov);
    bfGainIntraRx = bfGainIntRx(Icov);
    % Omnidirectional intra-cell interference
    for iue = 1:nue
        ibs = Icell(iue);
        
        pl = pathLoss(ibs,iue);
        BsUei = find(Icell == ibs); % array of index
        BsUei(BsUei == iue) = [];   % remove the current one from the index set
        
        intraTxBfGain = bfGainIntraTx(Icov(ibs, BsUei));
        intraRxBfGain = bfGainIntraRx(Icov(ibs,iue));

        intraOpt.ueList = BsUei; % List of UEs in the current sell
        % Tx BF gain + a Constant Rejection due to receive beamforming
        icellIntf = ones(1,length(BsUei))*pl - intraTxBfGain - intraRxBfGain + bfICIReject; 
        
        intraOpt.pathLoss = icellIntf;
        IntraCellIntf(num2str(iue)) = intraOpt;
    end
    
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
        % Check this!!!! Compare with original code
        bfGainRxDL(ibs2,iue) = 10*log10( real(wrx'*QrxTot(:,:,icov)*wrx/nantBS)); % is this correct? this is just reinitializing the same location!
    end
    
    % Intra cell interference
    pl = pathLoss2(ibs,iue);
    BsUei = find(Icell == ibs); % array of index
    BsUei(BsUei == iue) = []; % remove the current one from the index set
    
    intraOpt.ueList = BsUei; % List of UEs in the current cell
    pathLoss = ones(1,length(BsUei))*pl;
    
    intraBfGainRxDL = zeros(1,length(BsUei));
    
    for iIntraUe = 1:length(BsUei)
        indCov = Icov(ibs,BsUei(iIntraUe));
        intraBfGainRxDL(iIntraUe) = 10*log10( real(wrx'*QrxTot(:,:,indCov)*wrx/nantBS)); 
    end
    
    intraOpt.pathLoss = pathLoss - intraBfGainRxDL; % Pathloss with the BF gains
    
    IntraCellIntf(num2str(iue)) = intraOpt;
end

pathLossDL = (pathLoss2 - bfGainRxDL)'; % Nue x Nbs

end


