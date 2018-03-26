% fdiffs.m

% Inputs: Number of nodes x2, step sizes x2, base flow, base flow derivatives x2 and wave number.
% Outputs: Finite difference matrices A and B, and smaller parts a, b, c, d, e, gamma, alpha and beta.

% This code takes the base flow and constructs the main finite difference matrices a, b, c, d, e. Then it 
% forms the blocks for the boundary conditions and main diagonals. Lastly it puts all these together into 
% two larger matrices. These matrices are only valid for the symmetric mode.



function [A,B,e,d,a,b,c,gamma,alpha,beta] = fdiffs(Nz,Ny,dz,dy,U,dzU,dyU,k)



% Main matrices



a=zeros(Ny-2,Nz-2); b=a; c=a; d=a; e=a; % calculate finite difference coefficients
for ii=1:Ny-2
    for jj=1:Nz-2                                                    
        a(ii,jj)=U(ii,jj)/(dy^2)-dyU(ii,jj)/dy;
        b(ii,jj)=-2*U(ii,jj)/(dy^2)-2*U(ii,jj)/(dz^2)-(k^2)*U(ii,jj);
        c(ii,jj)=U(ii,jj)/(dy^2)+dyU(ii,jj)/dy;
        d(ii,jj)=U(ii,jj)/(dz^2)-dzU(ii,jj)/dz;
        e(ii,jj)=U(ii,jj)/(dz^2)+dzU(ii,jj)/dz;
    end
end
alpha=1/(dy^2); beta=-2/(dy^2)-2/(dz^2)-k^2; gamma=1/(dz^2); 



% Boundary conditions and constant matrices



T1=diag(d(1,1:Nz-3),1)+diag(b(1,1:Nz-2)+c(1,1:Nz-2),0)+diag(e(1,2:Nz-2),-1); % construct T1                                
T1(1,1)=T1(1,1)+e(1,1); T1(Nz-2,Nz-2)=T1(Nz-2,Nz-2)+d(1,Nz-2); % z bcs

TNm1=diag(d(Ny-2,1:Nz-3),1)+diag(b(Ny-2,1:Nz-2)-a(Ny-2,1:Nz-2),0)+diag(e(Ny-2,2:Nz-2),-1); % construct Tnm1 
TNm1(1,1)=TNm1(1,1)+e(Ny-2,1); TNm1(Nz-2,Nz-2)=TNm1(Nz-2,Nz-2)+d(Ny-2,Nz-2); % z bcs

alpham=diag(alpha.*ones(Nz-2,1),0); % construct alpha

S=diag(beta.*ones(Nz-2,1),0)+diag(gamma.*ones(Nz-3,1),-1) +diag(gamma.*ones(Nz-3,1),1); % construct S            
S(1,1)=S(1,1)+gamma; S(Nz-2,Nz-2)=S(Nz-2,Nz-2)+gamma; % z bcs

S1=S+diag(alpha.*ones(Nz-2,1),0); % construct S1

SNm1=S-diag(alpha.*ones(Nz-2,1),0); % construct SNm1



% Varying matrices



Tk=zeros(Nz-2,Nz-2,Ny-2); Ak=zeros(Nz-2,Nz-2,Ny-2); Ck=zeros(Nz-2,Nz-2,Ny-2); % construct interior blocks 
for ii=2:Ny-3                                                         
    Tk(:,:,ii)=diag(b(ii,1:Nz-2),0)+diag(e(ii,2:Nz-2),-1)+diag(d(ii,1:Nz-3),1);                  
    Tk(1,1,ii)=Tk(1,1,ii)+e(ii,1); Tk(Nz-2,Nz-2,ii)=Tk(Nz-2,Nz-2,ii)+d(ii,Nz-2);
end
for ii=1:Ny-2                                              
    Ak(:,:,ii)=diag(a(ii,1:Nz-2),0); Ck(:,:,ii)=diag(c(ii,1:Nz-2),0);     
end



% Main matrices out of smaller blocks 



A=speye((Ny-2)*(Nz-2),(Ny-2)*(Nz-2)); B=speye((Ny-2)*(Nz-2),(Ny-2)*(Nz-2)); % assemble A and B out of blocks
for ii=2:Ny-2                                                          
    A((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=Tk(:,:,ii);
    A((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-2)*(Nz-2)+1:(ii-1)*(Nz-2))=Ck(:,:,ii);
    A((ii-2)*(Nz-2)+1:(ii-1)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=Ak(:,:,ii-1);
    B((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=S;
    B((ii-1)*(Nz-2)+1:(ii)*(Nz-2),(ii-2)*(Nz-2)+1:(ii-1)*(Nz-2))=alpham;
    B((ii-2)*(Nz-2)+1:(ii-1)*(Nz-2),(ii-1)*(Nz-2)+1:(ii)*(Nz-2))=alpham;
end
A(1:Nz-2,1:Nz-2)=T1; A((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2),(Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=TNm1;            
B(1:Nz-2,1:Nz-2)=S1; B((Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2),(Ny-3)*(Nz-2)+1:(Ny-2)*(Nz-2))=SNm1; 



end