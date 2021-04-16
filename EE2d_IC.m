function [r0,u0,v0,p0,E0,tEnd,cfl] = EE2d_IC(IC,x,y)
% Load a smooth gaussian profile to the 1-d Euler Equations.
%
%       Coded by Manuel A. Diaz, ENSMA, 2021.02.26.
% 
% Notation:
%  r = Density
%  u = Velocity in x-direction
%  p = pressure
%  E = Enerty
global gamma p_inf M_inf

switch IC
    case 1
        fprintf('Setting a Gaussian density wave\n');
        % Parameters
        mu_x=5; mu_y=2.5; sigma=2;
        % Selected primitives
        r0 = 1 + exp(-((x-mu_x).^2+(y-mu_y).^2)/(2*sigma))/sqrt(2*pi*sigma);
        u0 = zeros(size(x));
        v0 = zeros(size(x));
        E0 = r0.^(gamma-1)/(gamma-1);
        p0 = (gamma-1)*r0.*(E0-0.5*(u0.^2+v0.^2));
        % Evolution parameters
        tEnd=27.0; cfl=0.50;
    case 2
        fprintf('Setting an inviscid vortex on a uniform flow\n');
        % Parameter
        M_inf=2.0; epsilon=0.1; K=1.0; x0=10; y0=5; 
        % Vortex function
        psi = (epsilon/(2*pi)) * exp(0.5*(1-(K^2)*((x-x0).^2+(y-y0).^2)));
        % Selected primitives
        r0 = (1 - 0.5*(gamma-1) * psi.^2 ).^(1/(gamma-1));
        u0 = M_inf + K*y.*psi;
        v0 = -K*x.*psi;
        p0 = (r0).^gamma;
        E0 = (p0./r0)/(gamma-1) + 0.5*(u0.^2+v0.^2);
        % Evolution parameters
        tEnd=1000.0; cfl=0.10;
end
% If and outflow BC is requested
p_inf = min(p0);

% Reshape as a vector
r0=r0(:); u0=u0(:); v0=v0(:); E0=E0(:); p0=p0(:);

end % set IC