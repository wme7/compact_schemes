function RHS = EE3d_FD_RHS(q,~)
global gamma Dx Dy Dz method

switch method
    case 1
        % Get primites
        u = q(:,2)./q(:,1);
        v = q(:,3)./q(:,1);
        w = q(:,4)./q(:,1);
        p = (gamma-1)*(q(:,5)-0.5*(q(:,2).^2+q(:,3).^2+q(:,4).^2)./q(:,1));
        H = (q(:,5)+p)./q(:,1);

        % Compute RHS using conservative variables
        RHS =-Dx*[q(:,2),q(:,2).*u+p,q(:,2).*v  ,q(:,2).*w  ,q(:,2).*H] ...
             -Dy*[q(:,3),q(:,3).*u  ,q(:,3).*v+p,q(:,3).*w  ,q(:,3).*H] ...
             -Dz*[q(:,4),q(:,4).*u  ,q(:,4).*v  ,q(:,4).*w+p,q(:,4).*H];
    case 2
        % Get primitive variables
        r = q(:,1);
        u = q(:,2)./q(:,1);
        v = q(:,3)./q(:,1);
        w = q(:,4)./q(:,1);
        p = (gamma-1)*(q(:,5)-0.5*(q(:,2).^2+q(:,3).^2+q(:,4).^2)./q(:,1));
        H = (q(:,5)+p)./q(:,1);

        % Set boundary conditions on primitives
        u = reshape(u,[ny,nx,nz]);  u(:,1,:)=0;  u(:,nx,:)=0;  u=u(:);
        v = reshape(v,[ny,nx,nz]);  v(1,:,:)=0;  v(ny,:,:)=0;  v=v(:);
        w = reshape(w,[ny,nx,nz]);  w(:,:,1)=0;  w(:,:,nz)=0;  w=w(:);
        
        % Compute RHS using skewsymmetric form of df/dx
        RHS =-[...
            0.5*( r  ).*(Dx*u)+0.5*u.*(Dx*( r  ))+0.5*Dx*( r.*u  ), ...
            0.5*(r.*u).*(Dx*u)+0.5*u.*(Dx*(r.*u))+0.5*Dx*(r.*u.*u) + Dx*p, ...
            0.5*(r.*u).*(Dx*v)+0.5*v.*(Dx*(r.*u))+0.5*Dx*(r.*u.*v), ...
            0.5*(r.*u).*(Dx*w)+0.5*w.*(Dx*(r.*u))+0.5*Dx*(r.*u.*w), ...
            0.5*(r.*u).*(Dx*H)+0.5*H.*(Dx*(r.*u))+0.5*Dx*(r.*u.*H)] ...
            - [...
            0.5*( r  ).*(Dy*v)+0.5*v.*(Dy*( r  ))+0.5*Dy*( r.*v  ), ...
            0.5*(r.*v).*(Dy*u)+0.5*u.*(Dy*(r.*v))+0.5*Dy*(r.*v.*u), ...
            0.5*(r.*v).*(Dy*v)+0.5*v.*(Dy*(r.*v))+0.5*Dy*(r.*v.*v) + Dy*p, ...
            0.5*(r.*v).*(Dy*w)+0.5*w.*(Dy*(r.*v))+0.5*Dy*(r.*v.*w), ...
            0.5*(r.*v).*(Dy*H)+0.5*H.*(Dy*(r.*v))+0.5*Dy*(r.*v.*H)] ...
            - [...
            0.5*( r  ).*(Dz*w)+0.5*w.*(Dz*( r  ))+0.5*Dz*( r.*w  ), ...
            0.5*(r.*w).*(Dz*u)+0.5*u.*(Dz*(r.*w))+0.5*Dz*(r.*w.*u), ...
            0.5*(r.*w).*(Dz*v)+0.5*v.*(Dz*(r.*w))+0.5*Dz*(r.*w.*v), ...
            0.5*(r.*w).*(Dz*w)+0.5*w.*(Dz*(r.*w))+0.5*Dz*(r.*w.*w) + Dz*p, ...
            0.5*(r.*w).*(Dz*H)+0.5*H.*(Dz*(r.*w))+0.5*Dz*(r.*w.*H)];
end