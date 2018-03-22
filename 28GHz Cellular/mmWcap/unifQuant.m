function [del, var, xmean] = unifQuant(nbit) 
% unitQuant:  Finds the optimal nbit bit uniform quantizer of
% a uniform Gaussian variable.  The outputs del is the optimal 
% stepsize, var is the variance and xmean is the optimal
% mapping point.  To quantized value should be:
%
%   nlev2 = 2^(nbit-1);
%   iq = floor(x./del);
%   iq = min(nlev2-1, max(-nlev2, iq));
%   q = xmean(nlev2 + iq + 1);

nlev2 = 2^(nbit-1);

if (nbit == 1)
    xmean = sqrt(2/pi);
    xmean = [-xmean xmean]';
    del = 1;
    var = 1 - 2/pi;
    return
end

% Loop over possible stepsizes
delavg = sqrt(12*2^(-2*nbit));
deltest = linspace(0,4*delavg,400)';
ntest = length(deltest);

vartest = zeros(ntest,1);
for it = 1:ntest
    
    % Threshold levels
    del = deltest(it);
    x = [0:nlev2-1]'*del;
    
    % Cummulative moments
    %   fn(i) = int_0^x(i) u^n p(u)du
    % where p(u) is the unit Gaussian pdf
    expx = exp(-x.^2/2);
    a = 1/sqrt(2*pi);
    f0 = 0.5*[erf(x/sqrt(2)); 1]; % Q func.
    f1 = a*[1 - expx; 1];
    f2 = f0 - [a*x.*expx; 0];
    
    % Probability in each level
    prob = diff(f0) + 1e-8;
    
    % Mean in each interval
    xmean = diff(f1)./prob;
    
    % Second moment in each interval
    xvar = diff(f2)./prob - xmean.^2;
    
    % Total variance
    vartest(it) = 2*xvar'*prob;
    
end
% Find minimum variance
[var, imin] = min(vartest);

% Re-compute optimal delta and mean values within interval
del = deltest(imin);
x = [0:nlev2-1]'*del;

% Cummulative moments
%   fn(i) = int_0^x(i) u^n p(u)du
% where p(u) is the unit Gaussian pdf
expx = exp(-x.^2/2);
a = 1/sqrt(2*pi);
f0 = 0.5*[erf(x/sqrt(2)); 1];
f1 = a*[1 - expx; 1];

% Probability in each level
prob = diff(f0);

% Mean in each interval
xmean = diff(f1)./prob;
xmean = [-flipud(xmean); xmean];

return

%plot(deltest./delavg, vartest);
%grid on;