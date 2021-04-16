function RHS = EE2d_FD_RHS(q,~)
global gamma Dx Dy method

switch method
    case 1
        % Get primites
        u = q(:,2)./q(:,1);
        v = q(:,3)./q(:,1);
        p = (gamma-1)*(q(:,4)-0.5*(q(:,2).^2+q(:,3).^2)./q(:,1));
        H = (q(:,4)+p)./q(:,1);

        % Compute RHS using conservative variables
        RHS = - Dx*[q(:,2),q(:,2).*u+p,q(:,2).*v  ,q(:,2).*H] ...
              - Dy*[q(:,3),q(:,3).*u  ,q(:,3).*v+p,q(:,3).*H];
    case 2
        % Get primitive variables
        r = q(:,1);
        u = q(:,2)./q(:,1);
        v = q(:,3)./q(:,1);
        p = (gamma-1)*(q(:,4)-0.5*(q(:,2).^2+q(:,3).^2)./q(:,1));
        H = (q(:,4)+p)./q(:,1);
        
        % Compute RHS using skew-symmetric form of df/dx
        RHS = - [...
            0.5*( r  ).*(Dx*u)+0.5*u.*(Dx*( r  ))+0.5*Dx*( r.*u  ), ...
            0.5*(r.*u).*(Dx*u)+0.5*u.*(Dx*(r.*u))+0.5*Dx*(r.*u.*u) + Dx*p, ...
            0.5*(r.*u).*(Dx*v)+0.5*v.*(Dx*(r.*u))+0.5*Dx*(r.*u.*v), ...
            0.5*(r.*u).*(Dx*H)+0.5*H.*(Dx*(r.*u))+0.5*Dx*(r.*u.*H)] ...
            - [...
            0.5*( r  ).*(Dy*v)+0.5*v.*(Dy*( r  ))+0.5*Dy*( r.*v  ), ...
            0.5*(r.*v).*(Dy*u)+0.5*u.*(Dy*(r.*v))+0.5*Dy*(r.*v.*u), ...
            0.5*(r.*v).*(Dy*v)+0.5*v.*(Dy*(r.*v))+0.5*Dy*(r.*v.*v) + Dy*p, ...
            0.5*(r.*v).*(Dy*H)+0.5*H.*(Dy*(r.*v))+0.5*Dy*(r.*v.*H)];
end