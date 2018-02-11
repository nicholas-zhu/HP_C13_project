function ratio = trans_calc_inv( A0, f0, mu, threshold)
% adiabatic pulse transaction area calculation
%
% Input:
%   A0 : adiabaitc pulse amplitude
%   f0 : adiabatic pulse bandwith
%   mu : adiabatic pulse para
%   
%   threshold: none zeros pulse

% basic paras
dt = 5e-6;
de_t = 1e-4;
% gamma_ba = 4.2585e3; %(Hz/G)
gamma_ba = 1.071e3;

% adiabatic inverse pulse design
dz_i = 2;
p_inv = pulse_HS( A0, mu, f0, 1e-2, dt);

% diffusion gradient design
GA_diff = 10;
g_diff = gradient_design( round( 5e-3  / dt ), GA_diff);

% adiabatic pulse
GA_inv = f0 / ( gamma_ba * dz_i );
g_inv = gradient_design( size( p_inv, 2), GA_inv);

p_de = time_delay( round( de_t/dt ) );
pulse = [  p_de, zeros(size(g_diff)), p_de, p_inv, p_de,zeros(size(g_diff)),p_de];
gz = [p_de, g_diff, p_de, g_inv, p_de, g_diff, p_de];
g = zeros( 3, length( gz ));
g( 3, :) = gz;

% bloch simulation
T1 = 14;
T2 = 1;
l_N = 2000;
l_z = linspace( -2, 2, l_N);
l = zeros( 3, length( l_z ));
l( 3, :) = l_z;

[ mx, my, mz] = bloch( pulse', g', dt, T1, T2, 0, l', 0, 0, 0, 1);
mz_f = medfilt1(mz,8);
ratio = 1 - mean( mz_f < (-threshold) )/mean(mz_f<threshold); % percentage of excitation
plot(l_z,mz_f);
xlabel('location/cm');
ylabel('Magnetization');