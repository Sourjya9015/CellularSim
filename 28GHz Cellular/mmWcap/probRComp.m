R = [0.001:0.001:0.3];

figure;
prob1 = min(0.018./R,1).*(1-exp(-R/0.063))+exp(-R./0.063);
prob2 = exp(-(R-0.01)/0.2);
prob3 = exp(-(R-0.01)/1.0);
plot(R,[prob1;prob2;prob3],'Linewidth',2);
hold on;

prob1 = min(0.018./R,1).*(1-exp(-R/0.072))+exp(-R./0.072);
prob2 = exp(-(R-0.01)/0.23);
prob3 = exp(-(R-0.01)/1.15);
plot(R,[prob1;prob2;prob3],'--','Linewidth',2);

prob1 = 0.5-min(0.5,5*exp(-0.156./R))+min(0.5,5*exp(-R/0.03));
prob2 = 0.5-min(0.5,3*exp(-0.3./R))+min(0.5,3*exp(-R/0.095));
plot(R,[prob1;prob2],'-.','Linewidth',2);

grid on;
set(gca,'FontSize',16);
xlabel('TX - RX separation, R (km)');
ylabel('probability of RX to be in LOS of TX, p(R)');
title('p(R) @ 2GHz');
xlim([R(1) R(end)]);
legend ('Macro to UE, case1: ISD=500m (urban)', ...
    'Macro to UE, case3: ISD=1732m, suburban', ...
    'Macro to UE, case3: ISD=1732m, rural/suburban',...
    'Macro to relay, case1: ISD=500m (urban)', ...
    'Macro to relay, case3: ISD=1732m, suburban', ...
    'Macro to relay, case3: ISD=1732m, rural/suburban',...
    'Relay to UE, case1: ISD=500m (urban)', ...
    'Relay to UE, case3: ISD=1732m');
