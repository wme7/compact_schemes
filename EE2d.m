%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Solving 2-D Euler Equations with 
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
global gamma Dx Dy idL idR idB idU IC method

%% Parameters
      IC = 01;	% 1:Inviscid Density Pulse, 2:Inviscid vortex
   gamma = 1.4;	% Ratio of specific heats for ideal di-atomic gas;
RKscheme = 04;	% 3:ERK3, 4:ERK4
FDscheme ='lele643'; % lele643, pade43, T4, T6, T8, E4, E6, E8.
  solver = 01;  % % 1:EE, 2:NS, 3:EE_charac, 4:NS_charac.
  method = 01;  % 1:Direct, 2:Skew-symmetric
plot_fig = true;

switch IC
    case 1, L=5.0; H=2.5; nx=061; ny=061; periodic_x=0; BCx=2; periodic_y=0; BCy=2; % Inviscid Density Pulse
    case 2, L=20.; H=10.; nx=101; ny=051; periodic_x=0; BCx=2; periodic_y=1; BCy=0; % Inviscid Vortex in/out
    otherwise, error('ERROR: IC not available');
end

% Discretize spatial domain
switch periodic_x
    case 0, xa=0; xb=L; lx=xb-xa; dx=(xb-xa)/(nx-1);
    case 1, xa=0; xb=L; lx=xb-xa; dx=(xb-xa)/(nx);
end
switch periodic_y
    case 0, ya=-H; yb=H; ly=yb-ya; dy=(yb-ya)/(ny-1);
    case 1, ya=-H; yb=H; ly=yb-ya; dy=(yb-ya)/(ny);
end
[x,y] = meshgrid(linspace(xa,xb,nx),linspace(ya,yb,ny));

% Set IC
[r0,u0,v0,p0,E0,tEnd,CFL] = EE2d_IC(IC,x,y); %tEnd=3.0; %-override
a0 = sqrt(gamma*p0./r0);        % Speed of sound
q0 = [r0,r0.*u0,r0.*v0,r0.*E0]; % array of conserved properties

% Discretize time domain
lambda0=max(abs(u0+a0)); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

% Select Solver
switch solver
    case 1, dF = @EE2d_FD_RHS; BC = @EE2d_BCs; disp('EE - simulation');
    %case 2, dF = @NS2d_FD_RHS; BC = @NS2d_BCs; disp('NS - simulation');
    %case 3, dF = @EE2d_FD_characteristics_RHS; BC = @EE2d_BCs; disp('EE - simulation');
    %case 4, dF = @NS2d_FD_characteristics_RHS; BC = @NS2d_BCs; disp('NS - simulation');
    otherwise, error('ERROR: Solver not available');
end

%% Build numerical scheme
% Build scheme
FD = compactSchemes(FDscheme,[nx,ny],[periodic_x,periodic_y]);

% Diff-operators
Dx = FD.Dx/dx;
Dy = FD.Dy/dy;

% Quadrature weights
W = FD.w;
A = FD.A;

% Boundary masks
idL = FD.index_L; idR = FD.index_R;
idB = FD.index_D; idU = FD.index_U;

%% Solver Loop

% If skew-symmetric
if method==2, FDscheme = [FDscheme,'_SS']; end

% Load initial condition
q=q0; it=0; dt=dt0; t=0; if plot_fig~=false, figure(3); end

% Integral of q from 0>x>L
int_q_dxdy=[]; time=[];

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
    r=q(:,1); u=q(:,2)./r; v=q(:,3)./r; E=q(:,4)./r; p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));
    if min(p)<0; error('negative pressure found!'); end; a=sqrt(gamma*p./r);
    
    % Update time step, dt
    lambda=max(sqrt(u.^2+v.^2)+a); dt=CFL*dx/lambda; 
    
    % Update iteration counter
	it=it+1;
    
    % Verify conservation
    int_q_dxdy(it,:)=W'*(A*q);  time(it)=t; %#ok<SAGROW>
    
    % Plot figure
    if plot_fig && rem(it,20) == 0
        switch IC
            case 1, mesh(x,y,reshape(q(:,4),[ny,nx])); axis square; view(35,25);
            case 2, mesh(x,y,reshape(  p   ,[ny,nx])); axis equal;  view(50,50);
        end, axis tight; drawnow
    end
end

%% Post-process
% compute flow properties
r=q(:,1); u=q(:,2)./r; v=q(:,3)./r; E=q(:,4)./r; p=(gamma-1)*r.*(E-0.5*(u.^2+v.^2));

%% Final plot
if ~exist('./figures','dir'), mkdir('./figures'); end
fig=figure(3); 
switch IC
    case 1, mesh(x,y,reshape(q(:,4),[ny,nx])); axis square; view(35,25); zlim([2.2,2.8]);
    case 2, mesh(x,y,reshape(  p   ,[ny,nx])); axis equal;  view(50,50); 
end, axis tight
xlabel('$x$','interpreter','latex','fontsize',20); 
ylabel('$y$','interpreter','latex','fontsize',20); 
zlabel('$\rho E$','interpreter','latex','fontsize',20);  
title([FDscheme,', time ',num2str(t)],'interpreter','latex','fontsize',20);
switch solver
    case 1, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE2d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny)],'-dpng');
    case 2, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_NS2d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny)],'-dpng');
end
fig=figure(4);
int_q_dxdy = int_q_dxdy./int_q_dxdy(1,:); % Normalize integration values
subplot(2,2,1); plot(time,int_q_dxdy(:,1),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho\,dx$','interpreter','latex');
subplot(2,2,3); plot(time,int_q_dxdy(:,2),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho u\,dx$','interpreter','latex');
subplot(2,2,4); plot(time,int_q_dxdy(:,3),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho v\,dx$','interpreter','latex');
subplot(2,2,2); plot(time,int_q_dxdy(:,4),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho E\,dx$','interpreter','latex');
switch solver
    case 1, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE2d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'_conservation'],'-dpng');
    case 2, print(fig,['./figures/',FDscheme,'_RK',num2str(RKscheme),'_NS2d_IC',num2str(IC),'_',num2str(nx),'x',num2str(ny),'_conservation'],'-dpng');
end