function q = EE1d_BCs(q,t)
    global IC BCs idL idR

    % Get primites variables
    r = q(:,1);
    u = q(:,2)./q(:,1);
    E = q(:,3)./q(:,1);

    % Set Dirichlet boundary conditions on primitives
    switch IC
        case 1 % Gaussian Pulse
        switch BCs
            case 1 % Rigid walls at x=0 and x=L.
                u(idL) = 0;
                u(idR) = 0;
            case 2 % Rigid walls at x=0.
                u(idL) = 0;
            case 3 % Rigid walls at x=L.
                u(idR) = 0;
            case 4 % Outflows at both ends.
                % Impose nothing !
        end
        case 2 % Test inflow BC
        switch BCs
            case 1 % Rigid walls at x=0 and x=L.
                u(idL) = 0;
                u(idL) = 0;
            case 2 % Inflow at x=0.
                %mu=0.10; sigma=2; u( 1) = mu*0.5*(1+tanh((t-2.5)*sigma));
                mu=0.05; sigma=2; u(idR) = mu*0.5*(1+tanh((t-2.5)*sigma))*cos(pi*t);
            case 3 % Inflow at x=L.
                %mu=0.10; sigma=2; u(nx) =-mu*0.5*(1+tanh((t-2.5)*sigma));
                mu=0.05; sigma=2; u(idL) =-mu*0.5*(1+tanh((t-2.5)*sigma))*cos(pi*t);
            case 4 % Outflows at both ends.
                % Impose nothing !
        end
    end
    
    % Get conservative variables
    q = [r, r.*u, r.*E];
end  