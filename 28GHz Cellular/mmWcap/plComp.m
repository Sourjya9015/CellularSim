% Plots path loss vs. distance comparing mmW and 3GPP Urban Micro model
d = linspace(10,1000,100)';  % distances in meters
nd = length(d);

% Path loss indices
MMW28 = 1;
MMW73 = 2;
PLF1 = 3;
PLF2 = 4;
UMi = 5;
FS28 = 6;
FS73 = 7;
LINK_BUD = 8;
nmod = LINK_BUD;
plStr = {'Empirical NYC, fc=28 GHz', 'Empirical NYC, fc=73 GHz', ...
         'PLF1', 'PLF2', '3GPP UMi, fc=2.5 GHz', ...
         'free-space, fc=28 GHz','free-space, fc=73 GHz'};
    
% Path loss models to plot
Iplot = [MMW28 MMW73 UMi FS28 FS73]';
nplot = length(Iplot);

% Path loss models
% ----------------
plMod = zeros(nd,nmod);

% mmW path loss based on Ted's groups data
load ../chanMod/data/chanParam28
plMod(:,MMW28) = plcoeff(1) + 10*plcoeff(2)*log10(d);

load ../chanMod/data/chanParam73
plMod(:,MMW73) = plcoeff(1) + 10*plcoeff(2)*log10(d);


% 3GPP Urban micro path loss model
fc1 = 2.5;
plMod(:,UMi) = 22.7 + 36.7*log10(d) + 26*log10(fc1);
plMod(:,PLF1) = 141.3 + 20*log10(d*1e-3);    % Samsung models
plMod(:,PLF2) = 157.4 + 32*log10(d*1e-3); 

% Free space propagation
c = 3e8;
fc = 28;
plMod(:,FS28) = 20*log10(4*pi*fc*1e9/c) + 20*log10(d);
fc = 73;
plMod(:,FS73) = 20*log10(4*pi*fc*1e9/c) + 20*log10(d);


% Maximum path loss @ SNR=0 dB, 1GHz
bw = 1e9;                   % BW in Hz
txpow = 30;                 % Tx power in dBm
BFgain = 10*log10(64^2)-3;  % BF gain
snrTgt = 0;                 % SNR target (dB)
noiseFig = 7;               % UE noise fig
plMod(:,LINK_BUD) = -snrTgt - noiseFig + BFgain+txpow + 174 - 10*log10(bw);

% Plot the results 
colors = (0:1/(nplot-1):1)';
linetype = {'-', '--', '-.'};
for iplot = 1:nplot
    ind = Iplot(iplot);
    linet = mod(iplot-1,3)+1;
    semilogx(d, plMod(:,ind), linetype{linet}, 'LineWidth', 3, ...
        'color', [1-colors(iplot) colors(iplot) 1]);
    hold on;
end
hold off;
set(gca,'FontSize',16);
xlabel('TX-RX separation (m)');
ylabel('Path loss (dB)');
grid on;
legend(plStr{Iplot},'Location', 'NorthWest');

return
print -depsc plots/plComp