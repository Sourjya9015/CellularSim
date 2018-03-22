% mmWave Uplink Power
%clc; clear all; close all;

%load data/user.mat
sinr_ul = [sinrTestUL{2,:}];
sinr_ul_f = [sinrTestUL{5,:}];
[UL_TDMA, axes] = hist (sinr_ul,50);
[UL_FDMA, xax] = hist (sinr_ul_f,50);
pdf_1 = UL_TDMA./length(sinr_ul);
pdf_2 = UL_FDMA./length(sinr_ul_f);

figure(5);
plot(axes, pdf_1,'-.');
hold on;
plot(xax, pdf_2, '-*r');
cdf_1 = cumsum(pdf_1);
cdf_2 = cumsum(pdf_2);
figure(6);
plot(axes,1-cdf_1,'*b');
hold on;
plot(xax,1-cdf_2,'*r');

% Fitting a Gaussian
mu = mean (sinr_ul);
sig = var(sinr_ul);

xplot = linspace(min(sinr_ul), max(sinr_ul), 1000)';

Fgauss = qfunc ((xplot-mu)./sqrt(sig));

plot (xplot,Fgauss,'--b','linewidth',2);

mu_f = mean (sinr_ul_f);
sig_f = var(sinr_ul_f);

xplot2 = linspace(min(sinr_ul_f), max(sinr_ul_f), 1000)';

Fgauss_f = qfunc ((xplot2-mu_f)./sqrt(sig_f));

plot (xplot2,Fgauss_f,'--r','linewidth',2);

disp (mu);
disp(mu_f);
