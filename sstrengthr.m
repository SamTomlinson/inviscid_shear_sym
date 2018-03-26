% sstrengthr.m

% Inputs: Numer of nodes x2, step size x2, grid vectors x2, grid mesh x2, base flow and derivatives, wavenumbers
% base flow singularity strength, phase speed, no-slip fraction, initial guess and top of domain.
% Outputs: Finite difference matrices D and d and the singularity strength ssR.

% This code takes the inputs then evaluates the singular solution on internal and external nodes. An initial guess
% is imported into the code from solving the full spectrum using eigs. The coefficient is the iterated upon
% by equating the singular and smoothed solutions untill convergence is achieved. % The solution and convergence 
% are then examined.

function [D,n,ssRn,vi,phat,flag] = sstrengthr(Nz,Ny,dz,dy,z,y,Z,Y,U,dzU,dyU,k,ssB,cn,delta,ssR0,h,r)


flag=0;
% Evalute the singular solution



% zv=[0-dz/2,z,1+dz/2]; yv=[0-dy/2,y,h+dy/2]; [Zv,Yv]=meshgrid(zv,yv); % full grid with ghost nodes

% figure() % plot the grid, interior and exterior nodes
% scatter(reshape(Zv,[Nz*Ny,1]),reshape(Yv,[Ny*Nz,1])); hold on; 
% rectangle('Position',[0 0 1 h])
% scatter(reshape(Z,[(Nz-2)*(Ny-2),1]),reshape(Y,[(Ny-2)*(Nz-2),1]))
% legend('Ghost nodes','Internal nodes')

R=sqrt((Z-delta).^2+Y.^2); T=atan2(Y,Z-delta); % convert to polars

vi=(R.*cos(T) - ( ssB*(R.^(3/2))/(9*cn) ).*( 2*sin(3*T/2)+3*(pi-T).*cos(3*T/2) )); % interior solution
vi2=vi-R.*cos(T); % interior solution

[vi2z,vi2y]=gradient(vi2,dz,dy); % calculate gradients for inhomogeneous terms



% Iterate upon the coefficient untill we achieve convergence



[A,B,e,d,a,~,~,gamma,alpha,~]=fdiffs(Nz,Ny,dz,dy,U,dzU,dyU,k); % calculate finite difference matrices 

ssRn=ssR0; tol=1;  svec=ssR0; % set up interation 



% Iterate upon the coefficient untill we achieve convergence


N=0;
while tol>1e-8

    tolo=tol; % tolerance check 

    ssRnm1=ssRn;  % update value of singularity strength

    viv=reshape(vi.',(Nz-2)*(Ny-2),1); viv2z=reshape(vi2z.',(Nz-2)*(Ny-2),1); viv2y=reshape(vi2y.',(Nz-2)*(Ny-2),1); % vectorise
    Uv=reshape(U.',(Nz-2)*(Ny-2),1); dyUv=reshape(dyU.',(Nz-2)*(Ny-2),1); dzUv=reshape(dzU.',(Nz-2)*(Ny-2),1);

    [D,n]=resid(A,B,cn,k,ssRnm1,viv,Uv,dyUv,dzUv,viv2y,viv2z,Nz,Ny,dz,a,d,e,alpha,gamma,h,dy,ssB,delta); % inhom matrices

    phat=D\n; phatm=reshape(phat,Nz-2,Ny-2).'; % invert system for phat and reshape

    phatmn=max(max(phatm)); rhsn=max(max(1+ssRnm1*vi+phatm)); othn=max(max(vi)); othn2=max(max(1+ssRnm1*vi)); % norms

    istar=1; jstar=round(delta*(Nz-2)); % at node nearest singularity
    ssRn=( (1+ssRnm1*vi(istar,jstar)+phatm(istar,jstar))/rhsn- 1 )/...
        ( vi(istar,jstar) - othn*(1+ssRnm1*vi(istar,jstar)+phatm(istar,jstar))/rhsn ); % iterate upon strength using fact solution equal to singularity 

    tol=norm(ssRn-ssRnm1); % calculate new tol

    if N==100000
         flag=1;
         fprintf('Passed max number of iterations')
         break
    end

    svec=[svec,ssRn]; % calculate tolerance and store answers
    
    N=N+1;

end

sf=svec(end); % store final strength 



% Examine converged solution 



% figure() % plot part of solution without singularity 
% plot(z,real(phatm(1,:)))
% xlabel('x')
% ylabel('ubar(x,y<<1)')

%fracom=4;

% figure() % plot real part of solution without singularity near wall 
% contourf(z,y(1:round(end/fracom)),real(phatm(1:round(end/fracom),:)),20)
% xlabel('z')
% ylabel('y')
% colorbar
% figure() % plot imag of solution without singularity near wall 
% contourf(z,y(1:round(end/fracom)),imag(phatm(1:round(end/fracom),:)),20)
% xlabel('z')
% ylabel('y')
% colorbar
% 
% figure() % plot real convergence of A
% plot(1:length(svec),real(svec))
% xlabel('Iterates')
% ylabel('Real convergence')
% figure() % plot imag convergence of A
% plot(1:length(svec),imag(svec))
% xlabel('Iterates')
% ylabel('Imaginary convergence')



end