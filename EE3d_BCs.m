function q = EE3d_BCs(q,~)
    global IC idL idR idD idU idB idF

    % Get primites variables
    r = q(:,1);
    u = q(:,2)./q(:,1);
    v = q(:,3)./q(:,1);
    w = q(:,4)./q(:,1);
    E = q(:,5)./q(:,1);

    % Set Dirichlet boundary conditions on primitives
    switch IC
        case 1
            u(idL) = 0;	% slip-wall
            u(idR) = 0;	% slip-wall
            v(idD) = 0;	% slip-wall
            v(idU) = 0;	% slip-wall
            w(idB) = 0;	% slip-wall
            w(idF) = 0;	% slip-wall
        case 2
            u(idL) = M_inf;	%
            % outflow with m_inf at idR
    end
    % Get conservative variables
    q = [r, r.*u, r.*v, r.*w, r.*E]; 
end