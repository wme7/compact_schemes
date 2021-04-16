%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Solving 3-D Euler Equations with 
%                   Conservative compact schemes
%
%           coded by Manuel A. Diaz, manuel.ade'at'gmail.com 
%                 Institut PPRIME, ENSMA, 2020.12.20
% 
% Last modif: 15.04.2012.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refs: 
% [1] Brady, Peter T., and Daniel Livescu. "High-order, stable, and
%     conservative boundary schemes for central and compact finite
%     differences." Computers & Fluids 183 (2019): 84-101. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; %close all;
global gamma Dx Dy Dz idL idR idD idU idB idF IC method

%% Parameters
      IC = 01;	% 1:Inviscid Density Pulse, 2:Inviscid vortex (z-periodic)
   gamma = 1.4;	% Ratio of specific heats for ideal di-atomic gas;
RKscheme = 04;	% 3:ERK3, 4:ERK4
FDscheme ='lele643'; % lele643, pade43, T4, T6, T8, E4, E6, E8.
  solver = 01;  % % 1:EE, 2:NS, 3:EE_charac, 4:NS_charac.
  method = 01;  % 1:Direct, 2:Skew-symmetric
plot_fig = true;

switch IC
    case 1, L=5;  H=5;  W=5;     nx=41; ny=41; nz=41; % Gaussian test
            periodic_x=0; periodic_y=0; periodic_z=0; 
    case 2, L=20; H=10; W=10;    nx=41; ny=41; nz=41; % Taylor-Vortex test
            periodic_x=0; periodic_y=0; periodic_z=1;
    otherwise, error('ERROR: IC not available');
end

% Discretize spatial domain
switch periodic_x
    case 0, xa=0; xb=L; lx=xb-xa; dx=(xb-xa)/(nx-1);
    case 1, xa=0; xb=L; lx=xb-xa; dx=(xb-xa)/(nx);
end
switch periodic_y
    case 0, ya=0; yb=H; ly=yb-ya; dy=(yb-ya)/(ny-1);
    case 1, ya=0; yb=H; ly=yb-ya; dy=(yb-ya)/(ny);
end
switch periodic_z
    case 0, za=0; zb=W; lz=zb-za; dz=(zb-za)/(nz-1);
    case 1, za=0; zb=W; lz=zb-za; dz=(zb-za)/(nz);
end
[x,y,z]=meshgrid(linspace(xa,xb,nx),linspace(ya,yb,ny),linspace(za,zb,nz));

% Set IC
[r0,u0,v0,w0,p0,E0,tEnd,CFL] = EE3d_IC(IC,x,y,z); %tEnd=3.0; %-override
a0 = sqrt(gamma*p0./r0);               % Speed of sound
q0 = [r0,r0.*u0,r0.*v0,r0.*w0,r0.*E0]; % array of conserved properties

% Discretize time domain
lambda0=max(abs(u0+a0)); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

% Select Solver
switch solver
    case 1, dF = @EE3d_FD_RHS; BC = @EE3d_BCs; disp('EE - simulation');
    %case 2, dF = @NS3d_FD_RHS; BC = @NS3d_BCs; disp('NS - simulation');
    %case 3, dF = @EE3d_FD_characteristic_RHS; BC = @EE3d_BCs; disp('EE - simulation');
    %case 4, dF = @NS3d_FD_characteristic_RHS; BC = @NS3d_BCs; disp('NS - simulation');
    otherwise, error('ERROR: Solver not available');
end

%% Build numerical scheme
% Build scheme
FD = compactSchemes(FDscheme,[nx,ny,nz],[periodic_x,periodic_y,periodic_z]);

% Diff-operators
Dx = FD.Dx/dx;
Dy = FD.Dy/dy;
Dz = FD.Dz/dz;

% Quadrature weights
W = FD.w;
A = FD.A;

% Boundary masks
idL = FD.index_L; idR = FD.index_R;
idD = FD.index_D; idU = FD.index_U;
idB = FD.index_B; idF = FD.index_F;

%% Solver Loop

% If skew-symmetric
if method==2, FDscheme = [FDscheme,'_SS']; end

% Load initial condition
q=q0; it=0; dt=dt0; t=0; if plot_fig~=false, figure(5); end

% Integral of q from 0>x>L
int_q_dxdydz=[]; time=[]; 

while t<tEnd
    % Iteration local time
    if t+dt>tEnd; dt=tEnd-t; end; t=t+dt;
 
    switch RKscheme
         case 3 % SSP-RK33
            qo= q;
            q = qo+dt*dF(q,t);                    q=BC(q,t);
            q = 0.75*qo+0.25*(q+dt*dF(q,t+dt/2)); q=BC(q,t+dt/2);
            q = (qo+2*(q+dt*dF(q,t+dt)))/3;       q=BC(q,t+dt);
        case 4 % ERK4
            qo= q;                            L1=dF(q,t);
            q = qo+dt/2*L1; q=BC(q,t+0.5*dt); L2=dF(q,t+0.5*dt);
            q = qo+dt/2*L2; q=BC(q,t+0.5*dt); L3=dF(q,t+0.5*dt);
            q = qo+dt*L3;   q=BC(q,t+dt);     L4=dF(q,t+dt);
            q = qo+dt*(L1+2*(L2+L3)+L4)/6;     q=BC(q,t+dt);
        otherwise, error('ERROR: RK scheme not set :P');
    end
   
    % compute primary properties
    u=q(:,2)./q(:,1); v=q(:,3)./q(:,1); w=q(:,4)./q(:,1); rE=reshape(q(:,5),[ny,nx,nz]);
    p=(gamma-1)*(q(:,5)-0.5*(q(:,2).^2+q(:,3).^2+q(:,4).^2)./q(:,1)); P = reshape(p,[ny,nx,nz]);
    if min(p)<0; error('negative pressure found!'); end; a=sqrt(gamma*p./q(:,1));
    
    % Update time step, dt
    lambda=max(sqrt(u.^2+v.^2+w.^2)+a); dt=CFL*dx/lambda; 
    
    % Update iteration counter
	it=it+1;
    
    % Verify conservation
    int_q_dxdydz(it,:)=W'*(A*q);  time(it)=t; %#ok<SAGROW>
    
    % Plot figure
    if plot_fig && rem(it,50) == 0
        switch IC
            case 1, mesh(x(:,:,31),y(:,:,31),rE(:,:,31)); axis square; view(35,25);
            case 2, mesh(x(:,:,31),y(:,:,31),P (:,:,31)); axis equal;  view(50,50);
        end, axis tight; drawnow
    end
end

%% Post-process
% compute flow properties
r=q(:,1); u=q(:,2)./r; v=q(:,3)./r; w=q(:,4)./r; E=q(:,5)./r; 
p=(gamma-1)*(q(:,5)-0.5*(q(:,2).^2+q(:,3).^2+q(:,4).^2));

%% Final plot
if ~exist('./figures','dir'), mkdir('./figures'); end
fig=figure(5); 
switch IC
    case 1, mesh(x(:,:,31),y(:,:,31),rE(:,:,31)); axis square; view(35,25);
    case 2, mesh(x(:,:,31),y(:,:,31),P (:,:,31)); axis equal;  view(50,50);
end%, axis tight
xlabel('$x$','interpreter','latex','fontsize',20); 
ylabel('$y$','interpreter','latex','fontsize',20); 
switch IC
    case { 1 }, zlabel('$\rho E$','interpreter','latex','fontsize',20);
    case { 2 }, zlabel(  '$p$'   ,'interpreter','latex','fontsize',20);
end, model = {'EE','NS','EE-charac','NS-charac'};
title([FDscheme,'-',model{solver},', time ',num2str(t)],'interpreter','latex','fontsize',20);
switch solver
    case {1,3}, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE3d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'x',num2str(nz)],'-dpng');
    case {2,4}, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_NS3d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'x',num2str(nz)],'-dpng');
end
fig=figure(6);
int_q_dxdydz = int_q_dxdydz./int_q_dxdydz(1,:); % Normalize integration values
subplot(2,3,[1,1.2]); plot(time,int_q_dxdydz(:,1),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho\,dx$','interpreter','latex');
subplot(2,3,4); plot(time,int_q_dxdydz(:,2),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho u\,dx$','interpreter','latex');
subplot(2,3,5); plot(time,int_q_dxdydz(:,3),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho v\,dx$','interpreter','latex');
subplot(2,3,6); plot(time,int_q_dxdydz(:,4),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho w\,dx$','interpreter','latex');
subplot(2,3,[2.8,3]); plot(time,int_q_dxdydz(:,5),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho E\,dx$','interpreter','latex');
switch solver
    case {1,3}, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE3d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'x',num2str(nz),'_conservation'],'-dpng');
    case {2,4}, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_NS3d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'x',num2str(nz),'_conservation'],'-dpng');
end