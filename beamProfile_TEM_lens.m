function [x,y,ux,uy,uz,type]  = beamProfile_TEM_lens(n,beamWaist,diverg,~) %zv..input:num_photons,beamWidth,beamDiverg,'gaussian'


% Updated 10/31/11 By W. Cox- Start with a Gaussian radius ray bundle distribution
% on the x/y plane. No divergence. Pass this through a THIN LENS ray matrix
% transformation (x1 = x2, theta2 = (-1/f)x1 + theta1) to diverge the beam.
% Choose 1/f to give the desired DIVERGENCE when the ray at the BEAM WAIST
% enters the diverging lens.

% m=0;
% sigma=beamWaist/sqrt(2);
% radius = abs(sigma*randn(n,1) + m);
type=2;
randVals = rand(n,1);  
radius = beamWaist*sqrt(-log(randVals));    
invF = -diverg/beamWaist;
divAng = -invF.*radius;       % Polar divergence angle based on thin lens ray x form matrix                       
phiAng = (2*pi).* rand(n,1);       	% Uniform distribution azimuth angle

cosPhi = cos(phiAng);
sinPhi = sin(phiAng);

uz = cos(divAng);
sin_uz = sqrt(1-uz.^2);
ux = sin_uz.*cosPhi;
uy = sin_uz.*sinPhi;

x = radius.*cosPhi;
y = radius.*sinPhi;

% Normalize the pointing vectors -> ux^2 + uy^2 + uz^2 = 1^2
if abs(1-(ux.^2 + uy.^2 + uz.^2)) > 1e-11
    normLength = sqrt(ux.^2 + uy.^2 + uz.^2);
    ux = ux / normLength;
    uy = uy / normLength;
    uz = uz / normLength;
    %disp('Vector normalization wrong!');
end

 


% %sd = 1/e beam radius, E(b) = (1/e)E(0)
% muD = cos(diverg);
% randVals = rand(n,1);
% ang = (2*pi).* rand(n,1);
% cosPhi = cos(ang);
% sinPhi = sin(ang);
% x = r.*cosPhi;
% y = r.*sinPhi;
% uz = 1 - rand(n,1).*(1 - muD); % randomly sample values from [1, muD], ISOTROPIC over [0,acos(muD)] <- IS THIS RIGHT?
%hist3([x y],[100 100])
%hist3([x y])

end
