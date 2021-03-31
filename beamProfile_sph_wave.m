% This function produces a spherical wave
% no lense
% August 2016 by Zahra 
% Assumption: point source is on the center (0,0,0)
% the area of the cross section plane should be equal to plane and gaussian
% wave
% Assumption: the maximum angle spreading is similar to gaussian wave


function [x,y,ux,uy,uz,type]  = beamProfile_sph_wave(n,beamWaist,diverg)     
type=3;
%% x and y calculation ( the same as plane wave )
R=2.5*beamWaist;
%divAng = pi/2*rand(n,1);
radius=R;
r = sqrt(rand(n,1))*radius;
phiAng = (2*pi).* rand(n,1);       	
cosPhi = cos(phiAng);
sinPhi = sin(phiAng);
x = r.*cosPhi;
y = r.*sinPhi;

%% direction calculation 
%Irr=1/(2*pi*(1-cos(diverg)));
%uz=1-rand(n,1)/(2*pi*Irr);
uz=1-rand(n,1);   %teta between 0 and pi/2
%%uz = cos(divAng);
sin_uz = sqrt(1-uz.^2);
ux = sin_uz.*cosPhi;
uy = sin_uz.*sinPhi;
%teta=acos(uz);
% Normalize the pointing vectors -> ux^2 + uy^2 + uz^2 = 1^2
if abs(1-(ux.^2 + uy.^2 + uz.^2)) > 1e-11
    normLength = sqrt(ux.^2 + uy.^2 + uz.^2);
    ux = ux / normLength;
    uy = uy / normLength;
    uz = uz / normLength;
    %disp('Vector normalization wrong!');
end

end
