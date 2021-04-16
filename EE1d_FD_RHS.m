function RHS = EE1d_FD_RHS(q,~)
global gamma Dx method

switch method
    case 1
        % Get primites
        u = q(:,2)./q(:,1);
        p = (gamma-1)*(q(:,3)-0.5*(q(:,2).^2)./q(:,1));
        H = (q(:,3)+p)./q(:,1);

        % Compute RHS using conservative variables
        RHS = -Dx*[q(:,2),q(:,2).*u+p,q(:,2).*H];
    case 2
        % Get primites variables
        r = q(:,1);
        u = q(:,2)./q(:,1);
        p = (gamma-1)*(q(:,3)-0.5*(q(:,2).^2)./q(:,1));
        H = (q(:,3)+p)./q(:,1);
        
        % Compute RHS using skewsymmetric form of df/dx
        RHS = -[...
            0.5*( r  ).*(Dx*u)+0.5*u.*(Dx*( r  ))+0.5*Dx*( r.*u  ), ...
            0.5*(r.*u).*(Dx*u)+0.5*u.*(Dx*(r.*u))+0.5*Dx*(r.*u.*u) + Dx*p, ...
            0.5*(r.*u).*(Dx*H)+0.5*H.*(Dx*(r.*u))+0.5*Dx*(r.*u.*H)];
end