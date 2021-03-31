% This function produces a plane wave
% no lense
% the area of cross section plane is equal to the gaussian wave

function [x,y,ux,uy,uz,type]  = beamProfile_plane_wave(n,beamWaist)     
type=1;
R=2.5*beamWaist;
radius=R;
r = sqrt(rand(n,1))*radius;
phiAng = (2*pi).* rand(n,1);       	
cosPhi = cos(phiAng);
sinPhi = sin(phiAng);
x = r.*cosPhi;
y = r.*sinPhi;

% x=0;
% y=0;

ux=0;
uy=0;
uz=1;   % teta z=0



end





