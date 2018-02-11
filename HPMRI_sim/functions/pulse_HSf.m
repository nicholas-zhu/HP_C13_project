function pulse = pulse_HSf ( A0, mu, f0, ta, dt)
% hyperbolic adiabatic pulse
% for double spin echo sequence without filtering
%
% Inputs:
%       A0 : B1 Amplitude(G/cm)(0.4~1)
%       mu : sharpness adjustment(>15)
%       f0 : Bandwidth(Hz)(0.65e4)
%       ta : adiabatic pulse duration(s)(~1e-2)
%       dt : time resolution(s)(4e-6)
%
% Outputs:
%       pulse: pulse shape of HS adiabatic pulse.
%
% written by Xucheng Zhu, Jan. 2016.



% paras
N_pulse = round(ta/dt); % sampling point number
t_N = (-N_pulse/2+1:N_pulse/2)*dt; % time point
n = 1;
m = 1;
beta = f0*pi/(mu);
window = 1;

% pulse design
t_phase = abs(t_N).^n/max(abs(t_N)).^(n-1) .* sign(t_N);
t_amplitude = abs(t_N).^m/max(abs(t_N)).^(m-1) .* sign(t_N);
B_phase = 1i*mu*log(sech(beta.*t_phase));
B_amplitude = (A0*sech(beta.*t_amplitude));

B = B_amplitude.*exp(B_phase);
pulse = B.*window;


