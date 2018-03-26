% initray.m

% Inputs: Number of nodes x2, grid mesh x2, step sizes x2, base flow and derivatives x2, base flow singularity strength, 
% full solution, phase speed and  no-slip fraction.
% Outputs: Initial guess for the singularity strength for the Rayleigh equation.

% This code takes the inputs then evaluates an initial guess for the singularity strength. This involves a 
% number of things which we must calculate, namely the singular solution, full solution and their derivatives, equated 
% at the point nearest the singularity to calculate the strength.



function [ssR0,vi,dzvi,dyvi,viv,dzvi2v,dyvi2v,U0v,dzU0v,dyU0v]=initray(Nz0,Ny0,Z0,Y0,dz0,dy0,U0,dzU0,dyU0,ssB0,E0,c0,delta)



% Calculate polar domain



R=sqrt((Z0-delta).^2+Y0.^2); T=atan2(Y0,Z0-delta); % evaluate polar domain  



% Other solutions and base flow



vi=R.*cos(T) - ( ssB0*(R.^(3/2))/(9*c0) ).*( 2*sin(3*T/2) + 3*(pi-T).*cos(3*T/2) ); % interior singular solution 
[dzvi,dyvi]=gradient(vi,dz0,dy0); % gradients
viv=reshape(vi.',(Nz0-2)*(Ny0-2),1); % vectorise

vi2=vi-R.*cos(T); % interior solution 
[dzvi2,dyvi2]=gradient(vi2,dz0,dy0); % calculate gradients
dzvi2v=reshape(dzvi2.',(Nz0-2)*(Ny0-2),1); dyvi2v=reshape(dyvi2.',(Nz0-2)*(Ny0-2),1); % vectorise

U0v=reshape(U0.',(Nz0-2)*(Ny0-2),1); dzU0v=reshape(dzU0.',(Nz0-2)*(Ny0-2),1); dyU0v=reshape(dyU0.',(Nz0-2)*(Ny0-2),1); % vectorise



% Initial guess at strength 



istar1=1; jstar1=round(delta*(Nz0-2)); % node nearest singularity 

vin=max(max(vi)); pin=max(max(E0)); % normalisation values

ssR0=( E0(istar1,jstar1)/pin - 1 )/( vi(istar1,jstar1) - vin*E0(istar1,jstar1)/pin ); % initial guess at strength 




end