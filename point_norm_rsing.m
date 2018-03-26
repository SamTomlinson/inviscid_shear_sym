% point_norm_rsing.m

% Inputs: Domain height, no-slip fraction, stream wise wavenumber, initial grid sizes 0a, 0b & 0c, final grid size, tolerances for 
% most unstable mode, fraction of most unstable modes. 
% Outputs: Most unstable mode (eigenvalue and eigenvector).

% This code takes the inputs and solves the singular eigenvalue problem at three resolutions to get initial guesses for the most 
% unstable mode. From this it calculates initial guesses for the r-singularity strength and residuals for the Muller method. 
% It then cranks up the resolution and iterates on the eigenvalue and r-singularity strength until they have both converged.
  


%clear all; close all; clc



%% Input parameters of the problem 



h=5; hi=h; % top of the computational domain 
delta=0.5; % no-slip fraction defined by delta=a/b
k=0.1; % steam wise wavenumber 
fracom=2; % fraction of the contour plot to be viewed



%% Evaluation of the numerical grid 0a 



Nz0a=15+2; Ny0a=75+2; % initial number of nodes in z and y (+2 for ghost nodes)
dy0a=h/(Ny0a-2); dz0a=1/(Nz0a-2); % step sizes in z and y
z0a=linspace(dz0a/2,1-dz0a/2,Nz0a-2); y0a=linspace(dy0a/2,h-dy0a/2,Ny0a-2); [Z0a,Y0a]=meshgrid(z0a,y0a); % generate grid

zv0a=[0-dz0a/2,z0a,1+dz0a/2]; yv0a=[0-dy0a/2,y0a,h+dy0a/2]; [Zv0a,Yv0a]=meshgrid(zv0a,yv0a); % full grid with ghost nodes



%% Initial estimate 0a



U0a=imag( (2/pi)*acos( cos(pi*(Z0a+sqrt(-1)*Y0a)/2)/cos(pi*delta/2) ) ); % evaluate singular base flow 
[dzU0a,dyU0a]=gradient(U0a,dz0a,dy0a); % calculate gradients of base flow

% figure() % plot base flow
% contourf(z0a,y0a(1:round(end/fracom)),U0a(1:round(end/fracom),:),10);
% xlabel('z'); ylabel('y');
% colorbar

[A0a,B0a,e0a,d0a,a0a,b0a,c0a,gamma0a,alpha0a,beta0a]=fdiffs(Nz0a,Ny0a,dz0a,dy0a,U0a,dzU0a,dyU0a,k); % calculate fd matrices 

[V0a,D0a]=eig(full(A0a),full(B0a),'qz'); % solve full eigenvalue problem 
ev0a=diag(D0a); evr0a=real(ev0a); evi0a=imag(ev0a); % extract eigenvalues

% figure() % plot initial spectrum
% scatter(evr0a,evi0a); hold on;
% xlabel('Re(lambda)')
% ylabel('Im(lambda)')

[~,ii0a]=max(imag(ev0a)); ev0a=ev0a(ii0a); % extract eigenvalue with maximum imaginary part 
E0av=V0a(:,ii0a); E0a=reshape(E0av,Nz0a-2,Ny0a-2).'; % corresponding eigenmode
E0an=E0a/max(max(E0a)); E0arn=real(E0an); E0ain=imag(E0an); % normalised eigenmode

% figure() % plot real eigenvector                                                           
% contourf(z0a,y0a(1:round(end/fracom)),E0arn(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')
% 
% figure() % plot imaginary eigenvector                                                
% contourf(z0a,y0a(1:round(end/fracom)),E0ain(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')



%% Evaluation of the numerical grid 0b 



Nz0b=17+2; Ny0b=85+2; % initial number of nodes in z and y (+2 for ghost nodes)
dy0b=h/(Ny0b-2); dz0b=1/(Nz0b-2); % step sizes in z and y
z0b=linspace(dz0b/2,1-dz0b/2,Nz0b-2); y0b=linspace(dy0b/2,h-dy0b/2,Ny0b-2); [Z0b,Y0b]=meshgrid(z0b,y0b); % generate grid

zv0b=[0-dz0b/2,z0b,1+dz0b/2]; yv0b=[0-dy0b/2,y0b,h+dy0b/2]; [Zv0b,Yv0b]=meshgrid(zv0b,yv0b); % full grid with ghost nodes

% figure() % plot the grid
% scatter(reshape(Zv0b,[(Nz0b)*(Ny0b),1]),reshape(Yv0b,[(Ny0b)*(Nz0b),1])); hold on;
% scatter(reshape(Z0b,[(Nz0b-2)*(Ny0b-2),1]),reshape(Y0b,[(Ny0b-2)*(Nz0b-2),1]))
% xlabel('z'); ylabel('y');
% rectangle('Position',[0 0 1 h])
% legend('Ghost nodes','Interior nodes')
% grid on



%% Initial estimate 0b



U0b=imag( (2/pi)*acos( cos(pi*(Z0b+sqrt(-1)*Y0b)/2)/cos(pi*delta/2) ) ); % evaluate base flow 
[dzU0b,dyU0b]=gradient(U0b,dz0b,dy0b); % calculate gradients of base flow 

% figure() % plot base flow
% contourf(z0b,y0b(1:round(end/fracom)),U0b(1:round(end/fracom),:),20);
% xlabel('z'); ylabel('y');
% colorbar

[A0b,B0b,e0b,d0b,a0b,b0b,c0b,gamma0b,alpha0b,beta0b]=fdiffs(Nz0b,Ny0b,dz0b,dy0b,U0b,dzU0b,dyU0b,k); % calculate fd matrices 

[V0b,D0b]=eig(full(A0b),full(B0b)); % solve full eigenvalue problem 
ev0b=diag(D0b); evr0b=real(ev0b); evi0b=imag(ev0b); % extract eigenvalues

% figure() % plot initial spectrum
% scatter(evr0b,evi0b); hold on;
% xlabel('Re(lambda)')
% ylabel('Im(lambda)')

[~,ii0b]=max(imag(ev0b)); ev0b=ev0b(ii0b); % extract eigenvalue with maximum imaginary part
E0bv=V0b(:,ii0b); E0b=reshape(E0bv,Nz0b-2,Ny0b-2).'; % corresponding eigenmode
E0bn=E0b/max(max(E0b)); E0brn=real(E0bn); E0bin=imag(E0bn); % normalised eigenmode

% figure() % plot real eigenvector
% contourf(z0b,y0b(1:round(end/fracom)),E0brn(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')
% 
% figure() % plot imaginary eigenvector 
% contourf(z0b,y0b(1:round(end/fracom)),E0bin(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')



%% Evaluation of the numerical grid 0c 



Nz0c=19+2; Ny0c=95+2; % initial number of nodes in z and y (+2 for ghost nodes)
dy0c=h/(Ny0c-2); dz0c=1/(Nz0c-2); % step sizes in z and y
z0c=linspace(dz0c/2,1-dz0c/2,Nz0c-2); y0c=linspace(dy0c/2,h-dy0c/2,Ny0c-2); [Z0c,Y0c]=meshgrid(z0c,y0c); % generate grid

zv0c=[0-dz0c/2,z0c,1+dz0c/2]; yv0c=[0-dy0c/2,y0c,h+dy0c/2]; [Zv0c,Yv0c]=meshgrid(zv0c,yv0c); % full grid with ghost nodes

% figure() % plot the grid
% scatter(reshape(Zv0c,[(Nz0c)*(Ny0c),1]),reshape(Yv0c,[(Ny0c)*(Nz0c),1])); hold on;
% scatter(reshape(Z0c,[(Nz0c-2)*(Ny0c-2),1]),reshape(Y0c,[(Ny0c-2)*(Nz0c-2),1]))
% xlabel('z'); ylabel('y');
% rectangle('Position',[0 0 1 h])
% legend('Ghost nodes','Interior nodes')
% grid on



%% Initial estimate 0c



U0c=imag( (2/pi)*acos( cos(pi*(Z0c+sqrt(-1)*Y0c)/2)/cos(pi*delta/2) ) ); % evaluate base flow 
[dzU0c,dyU0c]=gradient(U0c,dz0c,dy0c); % calculate gradients of base flow 

% figure() % plot base flow
% contourf(z0c,y0c(1:round(end/fracom)),U0c(1:round(end/fracom),:),20);
% xlabel('z'); ylabel('y');
% colorbar

[A0c,B0c,e0c,d0c,a0c,b0c,c0c,gamma0c,alpha0c,beta0c]=fdiffs(Nz0c,Ny0c,dz0c,dy0c,U0c,dzU0c,dyU0c,k); % calculate fd matrices 

[V0c,D0c]=eig(full(A0c),full(B0c)); % solve full eigenvalue problem 
ev0c=diag(D0c); evr0c=real(ev0c); evi0c=imag(ev0c); % extract eigenvalues

% figure() % plot initial spectrum
% scatter(evr0c,evi0c); hold on;
% xlabel('Re(lambda)')
% ylabel('Im(lambda)')

[~,ii0c]=max(imag(ev0c)); ev0c=ev0c(ii0c); % extract eigenvalue with maximum imaginary part
E0cv=V0c(:,ii0c); E0c=reshape(E0cv,Nz0c-2,Ny0c-2).'; % corresponding eigenmode
E0cn=E0c/max(max(E0c)); E0crn=real(E0cn); E0cin=imag(E0cn); % normalised eigenmode

% figure() % plot real eigenvector
% contourf(z0c,y0c(1:round(end/fracom)),E0crn(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')
% 
% figure() % plot imaginary eigenvector 
% contourf(z0c,y0c(1:round(end/fracom)),E0cin(1:round(end/fracom),:),20)
% colorbar
% xlabel('z')
% ylabel('y')



%% Calculate intitial guess for Rayleigh coefficient



ssB0a=sstrengthb(Nz0a,Ny0a,dz0a,dy0a,z0a,y0a,Z0a,Y0a,U0a,delta,h); % base flow singularity strength 
ssB0b=sstrengthb(Nz0b,Ny0b,dz0b,dy0b,z0b,y0b,Z0b,Y0b,U0b,delta,h);
ssB0c=sstrengthb(Nz0c,Ny0c,dz0c,dy0c,z0c,y0c,Z0c,Y0c,U0c,delta,h);

[ssR0a,via,dzvia,dyvia,viva,dzvi2va,dyvi2va,U0va,dzU0va,dyU0va]=initray(Nz0a,Ny0a,Z0a,Y0a,dz0a,dy0a,U0a,dzU0a,dyU0a,...
    ssB0a,E0a,ev0a,delta); % initial rayleigh singularity strength
[ssR0b,vib,dzvib,dyvib,vivb,dzvi2vb,dyvi2vb,U0vb,dzU0vb,dyU0vb]=initray(Nz0b,Ny0b,Z0b,Y0b,dz0b,dy0b,U0b,dzU0b,dyU0b,...
    ssB0b,E0b,ev0b,delta);
[ssR0c,vic,dzvic,dyvic,vivc,dzvi2vc,dyvi2vc,U0vc,dzU0vc,dyU0vc]=initray(Nz0c,Ny0c,Z0c,Y0c,dz0c,dy0c,U0c,dzU0c,dyU0c,...
    ssB0c,E0c,ev0c,delta);

% q=1+ssR0b*vib; q=q/max(max(q)); % examine (normalised) singular solution 
%
% figure() % plot v
% contourf(z0b,y0b(1:round(end/fracom)),real(q(1:round(end/fracom),:)),20)
% colorbar
% xlabel('z')
% ylabel('y')
% 
% figure() % plot v
% contourf(z0b,y0b(1:round(end/fracom)),imag(q(1:round(end/fracom),:)),20)
% colorbar
% xlabel('z')
% ylabel('y')



%% Evaluate residuals of intitial conditions



[C0a,bv0a]=resid(A0a,B0a,ev0a,k,ssR0a,viva,U0va,dyU0va,dzU0va,dyvi2va,dzvi2va,Nz0a,Ny0a,dz0a,a0a,d0a,e0a,alpha0a,gamma0a,...
    h,dy0a,ssB0a,delta); % calculate inhom matrices

ij0a=round(delta*(Nz0a-2)); % normalisation point

bv0an=bv0a; bv0an(ij0a)=1; % rhs vector
C0an=C0a; C0an(ij0a,:)=zeros(1,(Nz0a-2)*(Ny0a-2)); C0an(ij0a,ij0a)=1; % set point equal to one loosing one equation

p0an=C0an\bv0an; ff0a=C0a*p0an-bv0a; ff0ai=ff0a(ij0a); % calculate residual

[C0b,bv0b]=resid(A0b,B0b,ev0b,k,ssR0b,vivb,U0vb,dyU0vb,dzU0vb,dyvi2vb,dzvi2vb,Nz0b,Ny0b,dz0b,a0b,d0b,e0b,alpha0b,gamma0b,...
    h,dy0b,ssB0b,delta);  % calculate inhom matrices

ij0b=round(delta*(Nz0b-2)); % normalisation point

bv0bn=bv0b; bv0bn(ij0b)=1; % rhs vectors
C0bn=C0b; C0bn(ij0b,:)=zeros(1,(Nz0b-2)*(Ny0b-2)); C0bn(ij0b,ij0b)=1; % set point equal to one loosing one equation

p0bn=C0bn\bv0bn; ff0b=C0b*p0bn-bv0b; ff0bi=ff0b(ij0b); % calculate residual

[C0c,bv0c]=resid(A0c,B0c,ev0c,k,ssR0c,vivc,U0vc,dyU0vc,dzU0vc,dyvi2vc,dzvi2vc,Nz0c,Ny0c,dz0c,a0c,d0c,e0c,alpha0c,gamma0c,...
    h,dy0c,ssB0c,delta);  % calculate inhom matrices

ij0c=round(delta*(Nz0c-2)); % normalisation point

bv0cn=bv0c; bv0cn(ij0c)=1; % rhs vectors
C0cn=C0c; C0cn(ij0c,:)=zeros(1,(Nz0c-2)*(Ny0c-2)); C0cn(ij0c,ij0c)=1; % set point equal to one loosing one equation

p0cn=C0cn\bv0cn; ff0c=C0c*p0cn-bv0c; ff0ci=ff0c(ij0c); % calculate residual

cn=[ev0a,ev0b,ev0c]; fn=[ff0ai,ff0bi,ff0ci]; % set up arrays with initial guesses

% figure() % plot residual
% contourf(z0a,y0a,imag(reshape(ff0a,Nz0a-2,Ny0a-2).'),20)
% colorbar
% xlabel('z')
% ylabel('y')



%% Iterate until converged using point normalisation for increased resolution



Nz=151+2; Ny=755+2; % initial number of nodes in z and y (+2 for ghost nodes)
dy=h/(Ny-2); dz=1/(Nz-2); % step sizes in z and y
z=linspace(dz/2,1-dz/2,Nz-2); y=linspace(dy/2,h-dy/2,Ny-2); [Z,Y]=meshgrid(z,y); % generate grid

U=imag((2/pi)*acos(cos(pi*(Z+sqrt(-1)*Y)/2)/cos(pi*delta/2))); % evaluate base flow 
[dzU,dyU]=gradient(U,dz,dy); % calculate gradients of base flow 

ssB=sstrengthb(Nz,Ny,dz,dy,z,y,Z,Y,U,delta,h); % base flow singularity strength 

ii=4; tol=1e-6; ssRn=ssR0b; r=1; svec=[ssR0b]; % set counter and tolerance and strength store 

ij=round(delta*(Nz-2)); % normalisation node

while norm(cn(ii-1)-cn(ii-2))>tol % whilst difference of estimates is less than tolerance interate
    
    ssRnm1=ssRn; % update singularity strength 
    
    %cn=[cn,cn(ii-1)-0.1*fn(ii-1)*(cn(ii-1)-cn(ii-2))/(fn(ii-1)-fn(ii-2))]; % calculate new estimate using secant method

    w=fn(ii-1)/(cn(ii-1)-cn(ii-2))+fn(ii-2)/(cn(ii-2)-cn(ii-1))+fn(ii-1)/(cn(ii-1)-cn(ii-3))...
        +fn(ii-3)/(cn(ii-3)-cn(ii-1))-fn(ii-2)/(cn(ii-2)-cn(ii-3))-fn(ii-3)/(cn(ii-3)-cn(ii-2));
    fdd=fn(ii-1)/((cn(ii-1)-cn(ii-2))*(cn(ii-1)-cn(ii-3)))+fn(ii-2)/((cn(ii-2)-cn(ii-1))*(cn(ii-2)-cn(ii-3)))...
        +fn(ii-3)/((cn(ii-3)-cn(ii-1))*(cn(ii-3)-cn(ii-2)));
    if abs(w+sqrt(w^2-4*fn(ii-1)*fdd))>abs(w-sqrt(w^2-4*fn(ii-1)*fdd))
        cnnew=cn(ii-1)-0.1*2*fn(ii-1)/(w+sqrt(w^2-4*fn(ii-1)*fdd));
    else
        cnnew=cn(ii-1)-0.1*2*fn(ii-1)/(w-sqrt(w^2-4*fn(ii-1)*fdd));
    end
    cn=[cn,cnnew];
    
    
    [Dn,dn,ssRn,vi,phat,flag]=sstrengthr(Nz,Ny,dz,dy,z,y,Z,Y,U,dzU,dyU,k,ssB,cn(end),delta,ssRnm1,h,r); % new strength 
    
    test=interp2(Z,Y,ones(Ny-2,Nz-2)+ssRn*vi+reshape(phat,Nz-2,Ny-2).',Z0b,Y0b,'spline'); % calculate full solution
    testn=test/(max(max(test))); % normalise 
    
    err=norm(E0b/max(max(E0b))-testn); % calculate error estimes
    
    if flag==1 % exit flag to break loop
        break
    end
    
    bvn=dn; bvn(ij)=1; % normalise rhs vector
    Dnn=Dn; Dnn(ij,:)=zeros(1,(Nz-2)*(Ny-2)); Dnn(ij,ij)=1; % set point equal to one loosing one equation
    pn=Dnn\bvn; ff=Dn*pn-dn; ffi=ff(ij); % calculate residual 
    
    fn=[fn,ffi]; svec=[svec,ssRn]; % store new functional value and strength 
    
    ii=ii+1; % proceed to new iteration 
    
end


% fracom=2; 
% 
% figure() % plot convergence 
% plot(3:ii-1,real(cn(3:end))/real(cn(end))); hold on;
% plot(3:ii-1,imag(cn(3:end))/imag(cn(end)))
% xlabel('Iterations'); ylabel('Growth rate')  
% legend('Real','Imaginary')
% 
% figure() % plot convergence of coefficient 
% plot(3:ii-2,real(svec(3:end))/real(svec(end))); hold on;
% plot(3:ii-2,imag(svec(3:end))/imag(svec(end)))
% xlabel('Iterations'); ylabel('Rayleigh singular coefficient')  
% legend('Real','Imaginary')
% 
% figure() % plot constructed solution
% contourf(z0b,y0b(1:round(length(y0b)/fracom)),real(testn(1:round(length(y0b)/fracom),:)),20)
% xlabel('z')
% ylabel('y')
% colorbar
% 
% figure() % eigenmode
% contourf(z0b,y0b(1:round(length(y0b)/fracom)),real(E0b(1:round(length(y0b)/fracom),:)/max(max(E0b))),20)
% colorbar
% 
% figure() % plot constructed solution
% contourf(z0b,y0b(1:round(length(y0b)/fracom)),imag(testn(1:round(length(y0b)/fracom),:)),20)
% xlabel('z')
% ylabel('y')
% colorbar
% 
% figure() % eigenmode
% contourf(z0b,y0b(1:round(length(y0b)/fracom)),imag(E0b(1:round(length(y0b)/fracom),:)/max(max(E0b))),20)
% colorbar


% %% other modes plotting 
% 
% fracom=round(3);
% 
% phatm=abs(reshape(phat,Nz-2,Ny-2).')/max(abs(phat));
% [dzphat,dyphat]=gradient(phatm,dz,dy);
% vhat=-dyphat./(sqrt(-1)*k*(U-cn(end))); vhatn=abs(vhat)/max(max(abs(vhat)));
% what=-dzphat./(sqrt(-1)*k*(U-cn(end))); whatn=abs(what)/max(max(abs(what)));
% vhati=interp2(Z,Y,vhatn,Z0a,Y0a); whati=interp2(Z,Y,whatn,Z0a,Y0a); 
% [~,dyvhat]=gradient(vhatn,dz,dy); [dzwhat,~]=gradient(whatn,dz,dy);
% uhat=(ones(Ny-2,Nz-2)./(sqrt(-1)*k*(U-cn(end)))).*(-sqrt(-1)*k*phatm - vhat.*dyU - what.*dzU); uhatn=abs(uhat)/max(max(abs(uhat)));
% 
% %% contour plot of couette flow
% 
% figure() % plot constructed solution
% [c,h]=contourf(z,y(1:round(length(y)/fracom)),phatm(1:round(length(y)/fracom),1:end),5,'LineWidth',2); hold on; %clabel(c,h, 'labelspacing', 250,'fontsize',20);
% colormap(summer)
% quiver(Z0a(1:round(length(y0a)/fracom),1:end),Y0a(1:round(length(y0a)/fracom),1:end),vhati(1:round(length(y0a)/fracom),1:end),whati(1:round(length(y0a)/fracom),1:end),'k'); hold on;
% plot(linspace(delta,1,100),zeros(length(linspace(delta,1,100))),'LineWidth',5,'Color','k')
% set(gca,'fontsize',20)
% xlabel('$z$','FontSize',30,'Interpreter','latex'); 
% ylabel('$y$','FontSize',30,'Interpreter','latex');
% l1=legend('$ \ \hat{p}(z,y)$','$ \ \hat{v}/\hat{w}(z,y)$');
% set(l1,'FontSize',30,'Interpreter','latex')
% set(gcf, 'Position', [0, 0, 600, 900])
% ylim([dy/2,hi/fracom-dy])
% xlim([dz/2,1-dz/2])
% colorbar
% figure() % plot constructed solution
% [c,h]=contourf(z,y(1:round(length(y)/fracom)),uhatn(1:round(length(y)/fracom),1:end),5,'LineWidth',2); hold on; %clabel(c,h, 'labelspacing', 250,'fontsize',20);
% colormap(summer)
% quiver(Z0a(1:round(length(y0a)/fracom),1:end),Y0a(1:round(length(y0a)/fracom),1:end),vhati(1:round(length(y0a)/fracom),1:end),whati(1:round(length(y0a)/fracom),1:end),'k'); hold on;
% plot(linspace(delta,1,100),zeros(length(linspace(delta,1,100))),'LineWidth',5,'Color','k')
% set(gca,'fontsize',20)
% xlabel('$z$','FontSize',30,'Interpreter','latex'); 
% ylabel('$y$','FontSize',30,'Interpreter','latex'); 
% l1=legend('$ \ \hat{u}(z,y)$','$ \ \hat{v}/\hat{w}(z,y)$');
% set(l1,'FontSize',30,'Interpreter','latex')
% set(gcf, 'Position', [0, 0, 600, 900])
% ylim([dy/2,hi/fracom-dy])
% xlim([dz/2,1-dz/2])
% colorbar
% figure() % plot constructed solution
% [c,h]=contourf(z,y(1:round(length(y)/fracom)),vhatn(1:round(length(y)/fracom),1:end),5,'LineWidth',2);hold on; %clabel(c,h, 'labelspacing', 250,'fontsize',20); hold on;
% plot(linspace(delta,1,100),zeros(length(linspace(delta,1,100))),'LineWidth',5,'Color','k')
% set(gca,'fontsize',20)
% xlabel('$z$','FontSize',30,'Interpreter','latex'); 
% ylabel('$y$','FontSize',30,'Interpreter','latex'); 
% set(gcf, 'Position', [0, 0, 600, 900])
% figure() % plot constructed solution
% [c,h]=contourf(z,y(1:round(length(y)/fracom)),whatn(1:round(length(y)/fracom),1:end),5,'LineWidth',2); hold on; %clabel(c,h, 'labelspacing', 250,'fontsize',20); hold on;
% plot(linspace(delta,1,100),zeros(length(linspace(delta,1,100))),'LineWidth',5,'Color','k')
% set(gca,'fontsize',20)
% xlabel('$z$','FontSize',30,'Interpreter','latex'); 
% ylabel('$y$','FontSize',30,'Interpreter','latex'); 
% set(gcf, 'Position', [0, 0, 600, 900])






