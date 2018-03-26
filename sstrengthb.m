% sstrengthb.m

% Inputs: Number of nodes x2, step size x2, grid vectors x2, grid mesh x2, base flow, no-slip fraction and top of 
% domain.
% Outputs: Singularity strength for the base flow.

% This code takes the inputs then evaluates the singular solution on internal and external nodes. An initial guess
% is then evaluated by equating the full solution with local near the singularity. The coefficient is the iterated 
% upon by equating these two solutions until convergence is achieved. The solution and convergence are then 
% examined, and the error is calculated.



function [ssB] = sstrengthb(Nz,Ny,dz,dy,z,y,Z,Y,U,delta,h)



% Evalute the singular solution 



zv=[0-dz/2,z,1+dz/2]; yv=[0-dy/2,y,h+dy/2]; [Zv,Yv]=meshgrid(zv,yv); % full grid with ghost nodes

% figure() % plot the grid, interior and exterior nodes
% scatter(reshape(Zv,[Nz*Ny,1]),reshape(Yv,[Ny*Nz,1])); hold on;              
% scatter(reshape(Z,[(Nz-2)*(Ny-2),1]),reshape(Y,[(Ny-2)*(Nz-2),1]))
% rectangle('Position',[0 0 1 h])
% legend('Ghost nodes','Internal nodes')

Rv=sqrt((Zv-delta).^2+Yv.^2); R=sqrt((Z-delta).^2+Y.^2); % convert to polars
Tv=atan2(Yv,Zv-delta); T=atan2(Y,Z-delta);

v=(Rv.^(1/2)).*sin(Tv./2); % full singular solution
vi=(R.^(1/2)).*sin(T./2); % interior singular solution 

% figure() % plot  singular solution 
% contourf(z,y,vi,20)
% colorbar
% xlabel('x')
% ylabel('y')



% Calculate the intitial guess 



alpha=1/(dy^2); beta=-2/(dy^2)-2/(dz^2); gamma=1/(dz^2); % finite difference coefficients 

alpham=diag(alpha*ones(Nz-2,1),0); % diagonal matrix of alphas

S=diag(beta*ones(Nz-2,1),0)+diag(gamma*ones(Nz-3,1),-1)+diag(gamma*ones(Nz-3,1),1); % construct constant matrix 
S(1,1)=S(1,1)+gamma; S(Nz-2,Nz-2)=S(Nz-2,Nz-2)+gamma; % update due to bcs

S1=S; % first matrix with mixed bcs
S1(1:round(delta*(Nz-2))-1,1:round(delta*(Nz-2))-1)... % slip condition matrix
    =S(1:round(delta*(Nz-2))-1,1:round(delta*(Nz-2))-1)+diag(alpha*ones(round(delta*(Nz-2))-1,1),0); 
S1(round(delta*(Nz-2)):end,round(delta*(Nz-2)):end)... % no-slip condition matrix
    =S(round(delta*(Nz-2)):end,round(delta*(Nz-2)):end)-diag(alpha*ones(Nz-2-round(delta*(Nz-2))+1,1),0); 

SNm1=S+diag(alpha*ones(Nz-2,1),0); % end matrix due to bcs

B=speye((Ny-2)*(Nz-2),(Ny-2)*(Nz-2)); % enter block matrices down the three diagonals
for ii=2:Ny-2 
    B((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=S;
    B((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-2)*(Nz-2)+1:(ii-1)*(Nz-2))=alpham;
    B((ii-2)*(Nz-2)+1:(ii-1)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=alpham;
end          
B(1:Nz-2,1:Nz-2)=S1; B((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2),(Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=SNm1; % enter 1st and last mat

b=speye((Ny-2)*(Nz-2),1); b(1,1)=0; % inhomogeneous rhs vec
b((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=-alpha*dy*ones(length((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2)),1); 

ubar=B\b; ubarmat=reshape(ubar,Nz-2,Ny-2).'; % invert system for ubar and arrange to matrix form 
 
% figure() % plot ubar
% contourf(z,y,full(ubarmat),20)
% colorbar
% xlabel('z')
% ylabel('y')
% 
% figure() % plot solution nearest the wall 
% plot(z,ubarmat(1,:))
% xlabel('z')
% ylabel('u(x,y<<1)')

istar=1; jstar=round(delta*(Nz-2)); % choose point near singularity  

A0=ubarmat(istar,jstar)/vi(istar,jstar); % evaluate initial guess

Avec=A0; A1=A0; tol=1; % set up for iteration                                            



% Iterate upon the coefficient until we achieve convergence 



while tol>1e-10
    
    A0=A1; % update value of singularity strength    

    b=speye((Ny-2)*(Nz-2),1); b(1,1)=0; % enter varying block vectors components down vector
    for ii=1:Ny-2  
        b((ii-1)*(Nz-2)+1)=-gamma*A0*(v(ii+1,2)-v(ii+1,1)); % due to z bc
        b((ii)*(Nz-2))=gamma*A0*(v(ii+1,end)-v(ii+1,end-1)); % due to z bc
    end 
    b((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=b((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))... % due to y bc
        -dy*alpha + alpha*A0*(v(end,2:end-1).'-v(end-1,2:end-1).');

    ubar=B\b; ubarmat=reshape(ubar,Nz-2,Ny-2).'; % invert system for ubar and reshape 

    istar=1; jstar=round(delta*(Nz-2)); % at node nearest singularity 
    
    A1=A0+(ubarmat(istar,jstar)*vi(istar,jstar)...  % iterate upon strength using fact solution equal to singularity
        +ubarmat(istar,jstar+1)*vi(istar,jstar+1)...
        +ubarmat(istar,jstar-1)*vi(istar,jstar-1))/(vi(istar,jstar)...
        +vi(istar,jstar+1)+vi(istar,jstar-1)); 

    tol=norm(A1-A0); Avec=[Avec,A1]; % calculate tolerance and store answers

end

ssB=Avec(end); % save singularity strength



% Examine converged solution 



% figure() % plot part of solution without singularity 
% plot(z,A1*vi(1,:)+ubarmat(1,:))
% xlabel('x')
% ylabel('ubar(x,y<<1)')
% 
% figure() % plot part of solution without singularity near wall 
% contourf(z,y,A1*vi+ubarmat,20)
% xlabel('z')
% ylabel('y')
% colorbar
% 
% figure() % plot convergence of A
% plot(1:length(Avec),Avec)
% xlabel('Iterates')
% ylabel('Convergence')

err=norm(U-(A1*vi+ubarmat)); % calculate error as check (should be small)



end