function [mx,my,mz] = sequence_GRE(FA, ex, mx0, my0, mz0)

T1 = 40;
T2 = 40;
dt = 5e-6;
de_t = 8e-5;
diff_t = 1e-2;
gamma_ba = 1.071e3; %(Hz/G)


% pulse design
% excitation pulse design
dz = ex.dz;
BD_ex = ex.BD_ex;
p_ex = pulse_Gaus( FA, dt, BD_ex );
p_ex = pulse_Sinc( FA, dt, BD_ex, 32);


% gradient design
% selective pulse design
GA_ex = BD_ex / ( gamma_ba * dz );
g_ex = gradient_design( size( p_ex, 2), GA_ex);
g_ex2 = gradient_design( size( p_ex, 2),  - GA_ex / 2);

% pulse sequence
p_de = time_delay( round( de_t/dt ) );
pulse = [ p_de, p_ex, g_ex2 * 0, p_de];
gz = [p_de, g_ex, g_ex2, p_de];
g = zeros( 3, length( gz ));
g( 3, :) = gz;

% bloch simulation
l_N = size( mx0, 1);
l_z = linspace( -1, 1, l_N);
l = zeros( 3, length( l_z ));
l( 3, :) = l_z;

t_l = ones(1,length(g))*dt;
t_l(1,end) = 1;
[ mx, my, mz] = bloch( pulse', g', t_l, T1, T2, 0, l', 0, mx0, my0, mz0);
