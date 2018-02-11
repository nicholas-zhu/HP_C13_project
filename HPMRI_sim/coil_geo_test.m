% coil geometry and off resonance
% using bloch simulation
close all;
clear;
addpath /Users/xuchengzhu/Documents/Graduate/Research/Projects/HyperpolarizedMRI/HMRI_sim_C
addpath rf_tools/
addpath functions/

%% on resonance testing
dt = 1e-6;
de_t = 5e-3;
diff_t = 1e-3;
gamma_ba = 1.071e3;
% spin info
T1 = 40;
T2 = .2;
l_N = 1000;
l_Nx = 600;
l_z = linspace( -3, 3, l_N);
l_f = 100;
mx = zeros(l_N,1);
my = zeros(l_N,1);
mz = ones(l_N,1);
%% pulse design
% excitation pulse design
FA = 0*pi/180;
dz = .8;
BD_ex = 2e3;
p_ex = pulse_Sinc( FA, dt, BD_ex, 16);

df = 0;
p_ex = p_ex .* exp(-1i*2*pi*df*dt*linspace(-0.5,0.5,length(p_ex))*length(p_ex));

% adiabatic inverse pulse design
dz_i = 1.5;
BD_inv = 9e3;
t_inv = 1.5e-2;

temp = exp(- wthresh(linspace(-5,5,l_Nx),'s',3) .^ 2/9)';
window_inv = (min(temp,1)-min(temp))/(1-min(temp)+eps);
A_inv = window_inv*1.4;
A_inv = padarray(A_inv,(l_N-l_Nx)/2);

% gradient design
% selective pulse design
GA_ex = BD_ex / ( gamma_ba * dz );
g_ex = gradient_design( size( p_ex, 2), GA_ex);
g_ex2 = gradient_design( size( p_ex, 2),  - GA_ex/2);

% adiabatic pulse
GA_inv = BD_inv / ( gamma_ba * dz_i );
g_inv = gradient_design( round(t_inv/dt), GA_inv);

% diffusion gradient design
GA_diff = 4;
g_diff = gradient_design( round( 22e-3  / dt ), GA_diff);


%% simulation
N = 1;
for j = 1:N
    for i = 1:l_N
        % pulse sequence
        p_inv = pulse_HS(A_inv(i), 65, BD_inv, t_inv, dt);
        p_de = time_delay( round( de_t/dt ) );
        p_diff = time_delay( round( diff_t/dt ) );
        pulse = [ p_de, g_ex2*0 , p_ex, g_ex2*0, p_de, zeros(size(g_diff)), p_de, p_inv, p_de, zeros(size(g_diff)), p_diff, zeros(size(g_diff)), p_de, p_inv, p_de, zeros(size(g_diff)), p_de];
        gz = [p_de, g_ex2, g_ex, g_ex2, p_de, g_diff, p_de, g_inv, p_de, g_diff, p_diff, g_diff, p_de, g_inv, p_de, g_diff, p_de];
        g = zeros( 3, length( gz ));
        g( 3, :) = gz;

        [ mx0, my0, mz0] = bloch( pulse', g', dt, T1, T2, df, [0,0,l_z(i)], 0, mx(i), my(i), mz(i));
        mx(i) = mx0;
        my(i) = my0;
        mz(i) = mz0;
    end
    mx = medfilt1(mx,5);
    my = medfilt1(my,5);
    mz = medfilt1(mz,5);
end
plot([abs(mx(1:4:end)+1i*my(1:4:end)),mz(1:4:end)]);
