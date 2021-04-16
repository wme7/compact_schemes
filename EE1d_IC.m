function [r0,u0,p0,E0,tEnd,cfl] = EE1d_IC(x,IC)
% Load a smooth gaussian profile to the 1-d Euler Equations.
%
%       Coded by Manuel A. Diaz, ENSMA, 2021.02.26.
% 
% Notation:
%  r = Density
%  u = Velocity in x-direction
%  p = pressure
%  E = Enerty
global gamma p_inf
switch IC
    case 1
        fprintf('Density Gaussian Wave\n');
        % Parameters
        mu=5; sigma=2;
        % Selected primitives
        r0 = 1 + exp(-(x-mu).^2/(2*sigma))/sqrt(2*pi*sigma);
        u0 = zeros(size(x));
        E0 = r0.^(gamma-1)/(gamma-1);
        p0 = (gamma-1)*r0.*(E0-0.5*u0.^2);
        % Evolution parameters
        tEnd=10.5; cfl=0.50;
    case 2
        fprintf('Domain at rest\n');
        % Selected primitives
        r0 = ones(size(x));
        u0 = zeros(size(x));
        E0 = r0.^(gamma-1)/(gamma-1);
        p0 = (gamma-1)*r0.*(E0-0.5*u0.^2);
        % Evolution parameters
        tEnd=10.5; cfl=0.50;
end
% If and outflow BC is requested
p_inf = min(p0);

end % set IC