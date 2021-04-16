%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                 Solving 1-D Euler Equations with 
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
global gamma Dx idL idR IC BCs method

%% Parameters
      IC = 01;	% 1:Inviscid Density Pulse, 2:Inflow/Outflow
   gamma = 1.4;	% Ratio of specific heats for ideal di-atomic gas;
RKscheme = 04;	% 3:ERK3, 4:ERK4
FDscheme ='lele643'; % lele643, pade43, T4, T6, T8, E4, E6, E8.
  solver = 01;  % % 1:EE, 2:NS, 3:EE_charac, 4:NS_charac.
  method = 01;  % 1:Direct, 2:Skew-symmetric
plot_fig = true;

switch IC
    case 1, lx=5; nx=101; periodic_x=false; % 1:Inviscid Density Pulse
    case 2, lx=5; nx=101; periodic_x=false; % 2:Inviscid vortex
        BCs=01;  % 1:wall-wall, 2:wall/in-out, 3:out-wall/in, 4:out-out.
end

% Discretize spatial domain
switch IC
    case 1, xa=0; xb=lx; dx=(xb-xa)/(nx-1);   % periodic 
    case 2, xa=0; xb=lx; dx=(xb-xa)/(nx); % non-periodic
end, x=linspace(xa,xb,nx)';

% Set IC
[r0,u0,p0,E0,tFinal,CFL] = EE1d_IC(x,IC);
a0 = sqrt(gamma*p0./r0);     % Speed of sound
q0 = [r0, r0.*u0, r0.*E0];   % vec. of conserved properties

% Discretize time domain
lambda0=max(abs(u0)+a0); dt0=CFL*dx/lambda0;  % using the system's largest eigenvalue

% Select Solver
switch solver
    case 1, dF = @EE1d_FD_RHS; BC = @EE1d_BCs; % The component-wise solver
    %case 2, dF = @NS1d_FD_RHS; BC = @NS1d_BCs; % The component-wise solver
    %case 3, dF = @EE1d_FD_characteristics_RHS; BC = @EE1d_BCs;
	%case 4, dF = @NS1d_FD_characteristics_RHS; BC = @NS1d_BCs;
	otherwise, error('Solver now available.');
end

%% Build numerical scheme
% Build scheme
FD = compactSchemes(FDscheme,nx,periodic_x);

% Diff-operators
Dx = FD.Dx/dx;

% Quadrature weights
W = FD.w;
A = FD.A;

% Boundary masks
idL = FD.index_L; idR = FD.index_R;

%% Solver Loop

% If skew-symmetric
if method==2, FDscheme = [FDscheme,'_SS']; end

% Load initial condition
q=q0; it=0; dt=dt0; t=0; if plot_fig~=false, figure(1); end

% Integral of q from 0>x>L
int_q_dx=[]; time=[]; 

while t<tFinal
    % Iteration local time
    if t+dt>tFinal; dt=tFinal-t; end; t=t+dt;
 
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
    r=q(:,1); u=q(:,2)./r; E=q(:,3)./r; p=(gamma-1)*r.*(E-0.5*u.^2);
    if min(p)<0; error('negative pressure found!'); end; a=sqrt(gamma*p./r);
    
    % Update time step, dt
    lambda=max(abs(u)+a); dt=CFL*dx/lambda; 
    
    % Update iteration counter
	it=it+1;
    
    % Verify conservation
    int_q_dx(it,:)=W'*(A*q);  time(it)=t; %#ok<SAGROW>
    
    % Plot figure
    if plot_fig && rem(it,20) == 0
        subplot(221); plot(x,q(:,1),'.c'); axis([xa,xb, 0.8,1.4]);
        subplot(222); plot(x,q(:,2),'.b'); axis([xa,xb,-0.2,0.2]);
        subplot(223); plot(x, p    ,'.k'); axis([xa,xb, 0.8,1.6]);
        subplot(224); plot(x,q(:,3),'.m'); axis([xa,xb, 2.0,3.6]);
        drawnow
    end
end

%% Post-process
% compute flow properties
r=q(:,1); u=q(:,2)./r; E=q(:,3)./r; p=(gamma-1)*r.*(E-0.5*u.^2);

% Calculation of flow parameters
a = sqrt(gamma*p./r); M = u./a; % Mach number [-]
p_ref = 101325;           % Reference air pressure (N/m^2)
r_ref = 1.225;            % Reference air density (kg/m^3)
s_ref = 1/(gamma-1)*(log(p/p_ref)+gamma*log(r_ref./r)); 
                          % Entropy w.r.t reference condition
s = log(p./r.^gamma);     % Dimensionless Entropy
Q = r.*u;                 % Mass Flow rate per unit area
e = p./((gamma-1)*r);     % internal Energy

%% Final plots
if ~exist('./figures','dir'), mkdir('./figures'); end
figure(1);
subplot(221); plot(x, r  ,'.-r'); axis([xa,xb, 0.8,1.4]); xlabel('$x$','interpreter','latex','fontsize',20); ylabel('$\rho$','interpreter','latex','fontsize',20); 
subplot(222); plot(x, u  ,'.-r'); axis([xa,xb,-0.2,0.2]); xlabel('$x$','interpreter','latex','fontsize',20); ylabel('$u$','interpreter','latex','fontsize',20); 
subplot(223); plot(x, p  ,'.-r'); axis([xa,xb, 0.8,1.6]); xlabel('$x$','interpreter','latex','fontsize',20); ylabel('$p$','interpreter','latex','fontsize',20); 
subplot(224); plot(x, E  ,'.-r'); axis([xa,xb, 2.0,3.6]); xlabel('$x$','interpreter','latex','fontsize',20); ylabel('$E$','interpreter','latex','fontsize',20); 
switch solver
    case 1, print(['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE1d_IC',num2str(IC),'_Nx',num2str(nx)],'-dpng');
    case 2, print(['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE1d_IC',num2str(IC),'_Nx',num2str(nx)],'-dpng');
end
figure(2);
int_q_dx = int_q_dx./int_q_dx(1,:); % Normalize integration values
subplot(221); plot(time,int_q_dx(:,1),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho\,dx$','interpreter','latex');
subplot(223); plot(time,int_q_dx(:,2),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho u\,dx$','interpreter','latex');
subplot(222); plot(time,int_q_dx(:,3),'.-r'); xlabel('$t$','interpreter','latex'); ylabel('$\int_{0}^{L}\rho E\,dx$','interpreter','latex');
switch solver
    case 1, print(['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE1d_IC',num2str(IC),'_Nx',num2str(nx),'_conservation'],'-dpng');
    case 2, print(['./figures/',FDscheme,'_RK',num2str(RKscheme),'_EE1d_IC',num2str(IC),'_Nx',num2str(nx),'_conservation'],'-dpng');
end