function pulse = pulse_Sinc( FA, dt, f0, TBW)
% Sinc pulse 
%
% Inputs:
%       FA : flip angle.(rad)
%       f0 : frequency.(Hz)
%       dt : time resolution(s)(4e-6)
%       TBW : 
%
% Outputs:
%       pulse: pulse shape of Gaussian pulse.
%
% written by Xucheng Zhu, Jan. 2016.



% paras
gamma_ba = 1.071e3; % (G/Hz)
ta = TBW / f0 ;
N_rf = round( ta / dt );
t = linspace(- ta / 2, ta / 2, N_rf);

% pulse design
pulse = FA * f0 / ( gamma_ba * 2 * pi ) * sinc( t * f0 );