% this function produce CDF according to HG scattering phase function
% written by Cox, Edited by Zahra Vali
% august 2016 vancouver

function [cdf_scatter,petzold_angle_rad] = generate_scatter_HG(g)


%g = 0.924;

gsqr = g^2;
%theta = [[0:0.01:10] [10.1:0.1:180]];
theta =  0:0.01:180;

angle = theta.*pi./180; %radian

intensity = ((1/(4*pi)).*(1 - gsqr)) ./ (1 + gsqr - 2.*g.*cos(angle)).^(3/2);
cdf_scatter = 2*pi*cumtrapz(angle,intensity.*sin(angle));
cdf_scatter = cdf_scatter ./ max(cdf_scatter);
petzold_angle_rad = angle;

end


