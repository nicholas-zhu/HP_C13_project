function width = trans_est(t, A0, f0, mu, thk, thd)

% adiabatic pulse transaction area calculation
%
% Input:
%   t  : pulse duration (s)
%   A0 : adiabaitc pulse amplitude (G)
%   f0 : adiabatic pulse bandwith (Hz)
%   mu : adiabatic pulse para
%   thk: slice thickness (cm)
%   thd: none zeros pulse (%)
% Output:
%   width: width of transition area
% basic paras
dt = 5e-6;
t_d = 2e-2;
T1 = 20;
T2 = 2;
% gamma_ba = 4.2585e3; %(Hz/G)
gamma_ba = 1.071e3;
p_inv = pulse_HS ( A0, mu, f0, t, dt);
GA_ex = f0 / ( gamma_ba * thk );
g_inv = gradient_design( size( p_inv, 2), GA_ex);

% double spin echo
delay = zeros(1,round(t_d/dt));
gz = cat(2,g_inv,delay,g_inv);
pulse = cat(2,p_inv,delay,p_inv);
g = zeros( 3, length( gz ));
g( 3, :) = gz;

l_N = 10000;
l_z = linspace( -2, 2, l_N);%20mm
l = zeros( 3, length( l_z ));
l( 3, :) = l_z;

[ ~, ~, mz] = bloch( pulse', g', dt, T1, T2, 0, l', 0, 0, 0, 1);
width = sum(mz<thd*exp(-length(pulse)*dt/T1))/l_N*4;
plot(l_z,mz,'LineWidth',1);
xlabel('profile position/cm');
ylabel('Magnetization');
drawnow;
