function pulse = pulse_Gaus( FA, dt, f0)
% Gaussian pulse 
%
% Inputs:
%       FA : flip angle.(rad)
%       f0 : frequency.(Hz)
%       ta : adiabatic pulse duration(s)(~1e-2)
%       dt : time resolution(s)(4e-6)
%
% Outputs:
%       pulse: pulse shape of Gaussian pulse.
%
% written by Xucheng Zhu, Jan. 2016.



% paras
gamma_ba = 1.071e3; % (G/Hz)
ta = 1 / f0 ;
N_rf = round( ta / dt );
t_N = linspace( -9, 9, N_rf * 9);

% pulse design
pulse = 2 / ( gamma_ba * 2 * pi ) * f0 * exp(- pi * t_N .^ 2) * FA;