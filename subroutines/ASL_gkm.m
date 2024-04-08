function GKM = ASL_gkm(alpha,T1,Delta_t,tau,dt,t_max,m,f)

%--------------------------------------------------------------------------
% Input parameters:
%   alpha           % inversion efficacy
%   T1              % tissue T1
%   T1b             % perfusive substance T1
%   Delta_t         % transit delay
%   tau             % perfusion duration
%   m               % magnetization
%   f               % CBF ml of perfusive substance per ml tissue per sec
%   lambda          % partition coefficient water/tissue
%--------------------------------------------------------------------------
% vmalis@ucsd.edu
% Vadim Malis

lambda=1;
t = 0:dt:t_max; % Time vector
t(end)=[];
% Define the piecewise function c(t) for the pulsed state
c = zeros(size(t));
c(t >= Delta_t & t < tau + Delta_t) = alpha * exp(-t(t >= Delta_t & t < tau + Delta_t)/T1);

% Define r(t) and m(t)
r  = exp(-lambda*t/tau);

% r(t) and m(t)
rm = r.*m'; % Multiplying

% Convolve c(t) with rm_conv
Delta_M_conv = 2 * f * conv(c, rm, 'full') * dt; % Multiply by dt to approximate the integral

% Truncate the convolution result to match the time vector length
GKM = Delta_M_conv(1:length(t));
GKM(1:ceil(Delta_t/dt)) = 0;

end

