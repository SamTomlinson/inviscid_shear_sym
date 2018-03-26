% resid.m

% Inputs: Finite difference matrices, evalue, wavenumber, rayleigh strength, singular field, base flow and derivatives, 
% singular derivatives, node numbers, z step size, finite difference matrices, domain height, step size in y, base flow 
% singular strength, no-slip fraction
% Outputs: Inhomogeneous rhs vector and combined matrix.

% This code takes the inputs then evaluates the combined finite difference matrix C=A-cB. Then it moves on to evaluating 
% the rhs vector bv, which contains the singularity removal terms.



function [C,bv]=resid(A,B,ev,k,ssR,viv,Uv,dyUv,dzUv,dyvi2,dzvi2,Nz,Ny,dz,a,d,e,alpha,gamma,h,dy,ssB,delta)



% Calculate combined matrix



C=sparse(A-ev*B);



% Inhomogeneous vector



bv=sparse(((k^2)+(k^2)*ssR*viv).*(Uv-ev) + 2*ssR*(dyUv.*dyvi2) + 2*ssR*(dzUv.*dzvi2)); % from pde

zbl=linspace(0,1,Nz-2); ybl=linspace(dy/2,h-dy/2,Ny-2); [Zb,Yb]=meshgrid(zbl,ybl); % boundary points
Rb=sqrt((Zb-delta).^2+Yb.^2); Tb=atan2(Yb,Zb-delta); % boundary polars
vib=Rb.*cos(Tb) - ( ssB.*(Rb.^(3/2))/(9*ev) ).*( 2*sin(3*Tb/2)+3*(pi-Tb).*cos(3*Tb/2) ); % singular solution
[dzvib,~]=gradient(vib,zbl(2)-zbl(1),dy); % gradient

for ii=1:Ny-2 % from left and right bcs in z interior (not y)
    bv((ii-1)*(Nz-2)+1)=bv((ii-1)*(Nz-2)+1)-(e(ii,1)-ev*gamma)*ssR*dz*dzvib(ii,1); % due to z bc
    bv((ii)*(Nz-2))=bv((ii)*(Nz-2))+(d(ii,Nz-2)-ev*gamma)*ssR*dz*dzvib(ii,Nz-2); % due to z bc
end 

zbl=linspace(dz/2,1-dz/2,Nz-2); ybl=linspace(0,h,Ny-2); [Zb,Yb]=meshgrid(zbl,ybl); % boundary points
Rb=sqrt((Zb-delta).^2+Yb.^2); Tb=atan2(Yb,Zb-delta); % boundary polars
vib=Rb.*cos(Tb) - ( ssB*(Rb.^(3/2))/(9*ev) ).*( 2*sin(3*Tb/2)+3.*(pi-Tb).*cos(3*Tb/2) ); % singular solution

bv((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=bv((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))... % due to y free stream bc
    +2*(a(Ny-2,:)-ev*alpha).'.*(1+ssR*vib(Ny-2,:)');



end