function b = b_cal(grad,t1,t2,dt)

% gradient waveform [G/cm]
% t1, t2 inverse pulse time point [ms]
% dt [ms]
% b [s/mm^2]

N1 = floor(t1/dt)+1;
N2 = floor(t2/dt)+1;
gamma_ba = 1.07e3; % Hz/G

flag = ones(size(grad));
flag(N1:N2) = -1;

k = dt*cumsum(gamma_ba*grad.*flag)*1e-3;
b = (2*pi)^2*sum(k.^2)*dt;
b = b/1e5;