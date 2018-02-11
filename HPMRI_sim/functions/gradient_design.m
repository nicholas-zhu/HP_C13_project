function gradient = gradient_design ( N_g, A)
% gradient design
% 
% Inputs:
%   N_g : gradient points
%   A   : gradient intensity(G/cm)

gradient = ones( 1, N_g ) * A ;
