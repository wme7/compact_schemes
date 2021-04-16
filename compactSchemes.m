classdef compactSchemes
% Coded by: Manuel A. Diaz @ Pprime-ENSMA, 2021.
properties
    % Diff-Operators:
    dim         % Scheme dimensions
    scheme      % Name of the scheme
    nbs_id      % Number of scheme in table (only for E4,E6,E8,T4,T6,T8)
    periodic_x  % Perodicity in x-direction (true:1 if requested)
    Nx          % Size of nodes in x-direction
    Dx          % Derivative operator in x-direction
    periodic_y  % ... idem
    Ny
    Dy
    periodic_z  % ... idem
    Nz
    Dz
    % Quadrature:
    A           % Associated A matrix of a compact scheme.
    w           % Associated weights of the finite difference scheme.
    % Boundaries:
    index_L     % Boundary indexes at x=1  (left)
    index_R     % Boundary indexes at x=Nx (right)
    index_U     % Boundary indexes at y=1  (up)
    index_D     % Boundary indexes at y=Ny (down)
    index_B     % Boundary indexes at z=1  (front)
    index_F     % Boundary indexes at z=Nz (back)
end

methods
    
    function obj = compactSchemes(varargin)
        % 0. Set default parameters
        try
            obj.scheme = varargin{1};
            domainSize = varargin{2};
            obj.periodic_x = false;
            obj.periodic_y = false;
            obj.periodic_z = false;
            switch obj.scheme
                case 'E4', obj.nbs_id = 101; % : 4/3
                case 'E6', obj.nbs_id = 011; % : 6/5
                case 'E8', obj.nbs_id = 003; % : 8/7
                case 'T4', obj.nbs_id = 101; % : 4/3
                case 'T6', obj.nbs_id = 016; % : 6/5
                case 'T8', obj.nbs_id = 005; % : 8/7
                otherwise, obj.nbs_id = [];
            end
        catch
            disp("Usage: ");
            disp("  Object = compactSchemes(scheme,[Nx,Ny,Nz])");
            disp("  Object = compactSchemes(scheme,[Nx,Ny,Nz],[perodic_x,periodic_y,periodic_z])");
            disp("  Object = compactSchemes(scheme,[Nx,Ny,Nz],[perodic_x,periodic_y,periodic_z],nb_id)");
            disp("where:")
            disp("  periodic = true | false .")
            disp("  Schemes available: 'pade43','lele643','E4','E6','E8','T4','T6','T8'.");
            return
        end
        % 1. Identify the input parameters
        switch nargin
            case 2, obj.dim = numel(domainSize);
            case 3, obj.dim = numel(domainSize);
                    try obj.periodic_x = varargin{3}(1); catch, obj.periodic_x=[]; end
                    try obj.periodic_y = varargin{3}(2); catch, obj.periodic_y=[]; end
                    try obj.periodic_z = varargin{3}(3); catch, obj.periodic_z=[]; end
            case 4, obj.dim = numel(domainSize);
                    try obj.periodic_x = varargin{3}(1); catch, obj.periodic_x=[]; end
                    try obj.periodic_y = varargin{3}(2); catch, obj.periodic_y=[]; end
                    try obj.periodic_z = varargin{3}(3); catch, obj.periodic_z=[]; end
                    obj.nbs_id = varargin{4};
        end
        % 2 Identify the domain size
        switch obj.dim
            case 1, obj.dim= '1D'; obj.Nx = domainSize(1); nx=obj.Nx;
            case 2, obj.dim= '2D'; obj.Nx = domainSize(1); nx=obj.Nx;
                                   obj.Ny = domainSize(2); ny=obj.Ny;
            case 3, obj.dim= '3D'; obj.Nx = domainSize(1); nx=obj.Nx; 
                                   obj.Ny = domainSize(2); ny=obj.Ny; 
                                   obj.Nz = domainSize(3); nz=obj.Nz;
            otherwise, error('Cannot determine the dimensions of the domain.');
        end
        % 3. Set scheme's controling parameters
        switch obj.scheme 
            case 'pade43',  s=1; p=1; % : 4/3,  -> Ref. Moin (2007)
            case 'lele643', s=1; p=2; % : 6/4/3 -> Ref. Lele (2001)
            case 'E4',           p=2; % : 4/3, |
            case 'E6',           p=3; % : 6/5, |
            case 'E8',           p=4; % : 8/7, |-> Ref. Brady & Livescu (2019)
            case 'T4',      s=1; p=1; % : 4/3, |
            case 'T6',      s=1; p=2; % : 6/5, |
            case 'T8',      s=1; p=3; % : 8/7, |
            otherwise, error('Case not available.');
        end
        % 4 Build the desired scheme
        switch obj.dim
            
            case '1D'

                switch obj.scheme
                
                %---------------------
                case 'pade43' % 1D
                %---------------------
                
                % Discrete derivative based on a tri-diagonal scheme
                [beta_L,alpha_L] = TaylorTableAlgorithm(1,[  1 ],[ 0:2]); %#ok<NBRAK> % 3rd-order scheme
                [beta_C,alpha_C] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK> % 4th-order scheme
                [beta_R,alpha_R] = TaylorTableAlgorithm(1,[ -1 ],[-2:0]); %#ok<NBRAK> % 3rd-order scheme
                
                % Build A and B matrices
                diags = repmat([beta_C(1),1,beta_C(2)],[nx,1]);
                Ax = spdiags(diags,[-1:1],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Ax(1,1 + 1)    = beta_L; % Left  boundary
                        Ax(nx,nx - 1)  = beta_R; % Right boundary
                    case 1 % periodic
                        Ax( 1 ,nx ) = beta_C(2); % Left  boundary
                        Ax( nx, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Bx(1,1+( 0:2))   = alpha_L; % Left  boundary
                        Bx(nx,nx+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        Bx( 1 ,nx ) = alpha_C(1); % Left  boundary
                        Bx( nx, 1 ) = alpha_C(3); % Right boundary
                end
                %disp(full(Ax)); disp(full(Bx)); % Debug
                
                % Build Dx-operator
                obj.Dx = (Ax\Bx);

                % quadrature weights
                wx = ones(nx,1);
                
                % Integration weights
                obj.w = wx;
                obj.A = Ax;
                
                %---------------------
                case 'lele643' % 1D
                %---------------------
                    
                % Discrete derivative based on a tri-diagonal scheme
                [beta_L,alpha_L] = TaylorTableAlgorithm(1,[  1 ],[ 0:2]); %#ok<NBRAK> % 3rd-order scheme
                [beta_P,alpha_P] = TaylorTableAlgorithm(1,[-1,1],[-1:1]); %#ok<NBRAK> % 4rd-order scheme
                [beta_C,alpha_C] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK> % 6th-order scheme
                [beta_R,alpha_R] = TaylorTableAlgorithm(1,[ -1 ],[-2:0]); %#ok<NBRAK> % 3rd-order scheme

                % Discrete derivative using tri-diagonal scheme
                diags = repmat([beta_C(1),1,beta_C(2)],[nx,1]);
                Ax = spdiags(diags,[-1:1],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Ax( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Ax( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Ax(nx-1,[nx-2,nx]) = beta_P;
                        Ax( nx , nx-1 )    = beta_R;  % Right boundary c
                    case 1 % periodic
                        Ax( 1 ,nx ) = beta_C(2); % Left  boundary
                        Ax( nx, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Bx( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        Bx( 2 ,1+( 0:3))   = [alpha_P;0];
                        Bx(nx-1,nx+(-3:0)) = [0;alpha_P];
                        Bx( nx ,nx+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        Bx(1,nx-1:nx) = alpha_C(1:2); % Left  boundary
                        Bx(2,   nx  ) = alpha_C( 1 );
                        Bx(nx-1, 1  ) = alpha_C( 5 );
                        Bx( nx ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = (Ax\Bx);

                % quadrature weights
                wx = ones(nx,1);
                
                % Integration weights
                obj.w = wx;
                obj.A = Ax;
                
                %--------------------------
                case {'E4','E6','E8'} % 1D
                %--------------------------
                    
                % Discrete derivative using tri-diagonal scheme
                alpha = fornbergAlgorithm(1,0,[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                Ax = speye(nx);
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,[],alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = (Ax\Bx);
                
                % Integration weights
                obj.w = wx;
                obj.A = Ax;
                
                %--------------------------
                case {'T4','T6','T8'} % 1D
                %--------------------------
                    
                % Discrete derivative using tri-diagonal scheme
                [beta,alpha] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[nx,1]);
                Ax = spdiags(diags,[-s:s],nx,nx); %#ok<NBRAK>
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,beta,alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = (Ax\Bx);
                
                % Integration weights
                obj.w = wx;
                obj.A = Ax;
                    
                end % 1D schemes
                
                % Build boundary mask
                obj.index_L=1;  obj.index_R=nx;

            case '2D'

                switch obj.scheme
                
                %----------------------
                case 'lele643' % 2D
                %----------------------
                    
                % Discrete derivative based on a tri-diagonal scheme
                [beta_L,alpha_L] = TaylorTableAlgorithm(1,[  1 ],[ 0:2]); %#ok<NBRAK> % 3rd-order scheme
                [beta_P,alpha_P] = TaylorTableAlgorithm(1,[-1,1],[-1:1]); %#ok<NBRAK> % 4rd-order scheme
                [beta_C,alpha_C] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK> % 6th-order scheme
                [beta_R,alpha_R] = TaylorTableAlgorithm(1,[ -1 ],[-2:0]); %#ok<NBRAK> % 3rd-order scheme

                % Discrete derivative using tri-diagonal scheme
                diags = repmat([beta_C(1),1,beta_C(2)],[nx,1]);
                Ax = spdiags(diags,[-1:1],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Ax( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Ax( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Ax(nx-1,[nx-2,nx]) = beta_P;
                        Ax( nx , nx-1 )    = beta_R;  % Right boundary c
                    case 1 % periodic
                        Ax( 1 ,nx ) = beta_C(2); % Left  boundary
                        Ax( nx, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Bx( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        Bx( 2 ,1+( 0:3))   = [alpha_P;0];
                        Bx(nx-1,nx+(-3:0)) = [0;alpha_P];
                        Bx( nx ,nx+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        Bx(1,nx-1:nx) = alpha_C(1:2); % Left  boundary
                        Bx(2,   nx  ) = alpha_C( 1 );
                        Bx(nx-1, 1  ) = alpha_C( 5 );
                        Bx( nx ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = kron((Ax\Bx),speye(ny));

                % quadrature weights
                wx = ones(nx,1);

                % Build the A and B Matrices:
                diags = repmat([beta_C(1),1,beta_C(2)],[ny,1]);
                Ay = spdiags(diags,[-1:1],ny,ny); %#ok<NBRAK>
                switch obj.periodic_y
                    case 0 % non-periodic
                        Ay( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Ay( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Ay(ny-1,[ny-2,ny]) = beta_P;
                        Ay( ny , ny-1 )    = beta_R;  % Right boundary
                    case 1 % periodic
                        Ay( 1 ,ny ) = beta_C(2); % Left  boundary
                        Ay( ny, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                switch obj.periodic_y
                    case 0 % non-periodic
                        By( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        By( 2 ,1+( 0:3))   = [alpha_P;0];
                        By(ny-1,ny+(-3:0)) = [0;alpha_P];
                        By( ny ,ny+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        By(1,ny-1:ny) = alpha_C(1:2); % Left  boundary
                        By(2,   ny  ) = alpha_C( 1 );
                        By(ny-1, 1  ) = alpha_C( 5 );
                        By( ny ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Ay)); disp(full(By)); % Debug

                % Build Dx-operator
                obj.Dy = kron(speye(nx),(Ay\By));

                % quadrature weights
                wy = ones(ny,1);

                % Integration weights
                obj.w = kron(wx,wy);
                obj.A = kron(Ax,Ay);
                
                %--------------------------
                case {'E4','E6','E8'} % 2D
                %--------------------------

                % Discrete derivative using tri-diagonal scheme
                alpha = fornbergAlgorithm(1,0,[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                Ax = speye(nx);
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,[],alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = kron((Ax\Bx),speye(ny));

                % Build the A and B Matrices:
                Ay = speye(ny);
                diags = repmat(alpha',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_y
                    case 0, [Ay,By,wy] = conservBoundarySchemes(Ay,By,obj.scheme,obj.nbs_id);
                    case 1, [Ay,By,wy] = set_periodicBCs(Ay,By,[],alpha);
                end
                % disp(full(Ay)); disp(full(By)); % Debug

                % Build Dx-operator
                obj.Dy = kron(speye(nx),(Ay\By));

                % Integration weights
                obj.w = kron(wx,wy);
                obj.A = kron(Ax,Ay);
                
                %--------------------------
                case {'T4','T6','T8'} % 2D
                %--------------------------

                % Discrete derivative using tri-diagonal scheme
                [beta,alpha] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[nx,1]);
                Ax = spdiags(diags,[-s:s],nx,nx); %#ok<NBRAK>
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,beta,alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = kron((Ax\Bx),speye(ny));

                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[ny,1]);
                Ay = spdiags(diags,[-s:s],ny,ny); %#ok<NBRAK>
                diags = repmat(alpha',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_y
                    case 0, [Ay,By,wy] = conservBoundarySchemes(Ay,By,obj.scheme,obj.nbs_id);
                    case 1, [Ay,By,wy] = set_periodicBCs(Ay,By,beta,alpha);
                end
                % disp(full(Ay)); disp(full(By)); % Debug

                % Build Dx-operator
                obj.Dy = kron(speye(nx),(Ay\By));

                % Integration weights
                obj.w = kron(wx,wy);
                obj.A = kron(Ax,Ay);
                
                end % 2D schemes
                
                % Build boundary mask
                idx = reshape(1:nx*ny,[ny,nx]);
                obj.index_L = idx(:,1); obj.index_R = idx(:,nx);
                obj.index_D = idx(1,:); obj.index_U = idx(ny,:);

            case '3D'

                switch obj.scheme
                
                %----------------------
                case 'lele643' % 3D
                %----------------------
                    
                % Discrete derivative based on a tri-diagonal scheme
                [beta_L,alpha_L] = TaylorTableAlgorithm(1,[  1 ],[ 0:2]); %#ok<NBRAK> % 3rd-order scheme
                [beta_P,alpha_P] = TaylorTableAlgorithm(1,[-1,1],[-1:1]); %#ok<NBRAK> % 4rd-order scheme
                [beta_C,alpha_C] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK> % 6th-order scheme
                [beta_R,alpha_R] = TaylorTableAlgorithm(1,[ -1 ],[-2:0]); %#ok<NBRAK> % 3rd-order scheme

                % Discrete derivative using tri-diagonal scheme
                diags = repmat([beta_C(1),1,beta_C(2)],[nx,1]);
                Ax = spdiags(diags,[-1:1],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Ax( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Ax( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Ax(nx-1,[nx-2,nx]) = beta_P;
                        Ax( nx , nx-1 )    = beta_R;  % Right boundary c
                    case 1 % periodic
                        Ax( 1 ,nx ) = beta_C(2); % Left  boundary
                        Ax( nx, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                switch obj.periodic_x
                    case 0 % non-periodic
                        Bx( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        Bx( 2 ,1+( 0:3))   = [alpha_P;0];
                        Bx(nx-1,nx+(-3:0)) = [0;alpha_P];
                        Bx( nx ,nx+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        Bx(1,nx-1:nx) = alpha_C(1:2); % Left  boundary
                        Bx(2,   nx  ) = alpha_C( 1 );
                        Bx(nx-1, 1  ) = alpha_C( 5 );
                        Bx( nx ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = 0.5*kron(speye(nz),kron((Ax\Bx),speye(ny)));

                % quadrature weights
                wx = ones(nx,1);
                
                % Build the A and B Matrices:
                diags = repmat([beta_C(1),1,beta_C(2)],[ny,1]);
                Ay = spdiags(diags,[-1:1],ny,ny); %#ok<NBRAK>
                switch obj.periodic_y
                    case 0 % non-periodic
                        Ay( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Ay( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Ay(ny-1,[ny-2,ny]) = beta_P;
                        Ay( ny , ny-1 )    = beta_R;  % Right boundary
                    case 1 % periodic
                        Ay( 1 ,ny ) = beta_C(2); % Left  boundary
                        Ay( ny, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                switch obj.periodic_y
                    case 0 % non-periodic
                        By( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        By( 2 ,1+( 0:3))   = [alpha_P;0];
                        By(ny-1,ny+(-3:0)) = [0;alpha_P];
                        By( ny ,ny+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        By(1,ny-1:ny) = alpha_C(1:2); % Left  boundary
                        By(2,   ny  ) = alpha_C( 1 );
                        By(ny-1, 1  ) = alpha_C( 5 );
                        By( ny ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Ay)); disp(full(By)); % Debug
                
                % Build Dy-operator
                obj.Dy = 0.5*kron(speye(nz),kron(speye(nx),(Ay\By)));
                
                % quadrature weights
                wy = ones(ny,1);
                
                % Build the A and B Matrices:
                diags = repmat([beta_C(1),1,beta_C(2)],[nz,1]);
                Az = spdiags(diags,[-1:1],nz,nz); %#ok<NBRAK>
                switch obj.periodic_z
                    case 0 % non-periodic
                        Az( 1 , 1+1 )      = beta_L;  % Left  boundary
                        Az( 2 ,[1,3])      = beta_P;  %#ok<*SPRIX>
                        Az(nz-1,[nz-2,nz]) = beta_P;
                        Az( nz , nz-1 )    = beta_R;  % Right boundary
                    case 1 % periodic
                        Az( 1 ,nz ) = beta_C(2); % Left  boundary
                        Az( nz, 1 ) = beta_C(1); % Right boundary
                end
                diags = repmat(alpha_C',[nz,1]);
                Bz = spdiags(diags,[-p:p],nz,nz); %#ok<NBRAK>
                switch obj.periodic_z
                    case 0 % non-periodic
                        Bz( 1 ,1+( 0:2))   = alpha_L; % Left  boundary
                        Bz( 2 ,1+( 0:3))   = [alpha_P;0];
                        Bz(nz-1,nz+(-3:0)) = [0;alpha_P];
                        Bz( nz ,nz+(-2:0)) = alpha_R; % Right boundary
                    case 1 % periodic
                        Bz(1,nz-1:nz) = alpha_C(1:2); % Left  boundary
                        Bz(2,   nz  ) = alpha_C( 1 );
                        Bz(nz-1, 1  ) = alpha_C( 5 );
                        Bz( nz ,1:2 ) = alpha_C(4:5); % Right boundary
                end
                % disp(full(Az)); disp(full(Bz)); % Debug
                
                % Build Dz-operator
                obj.Dz = 0.5*kron((Az\Bz),kron(speye(nx),speye(ny)));
                
                % quadrature weights
                wz = ones(nz,1);
                
                % Integration weights
                obj.w = kron(wz,kron(wx,wy));
                obj.A = kron(Az,kron(Ax,Ay));
                
                %--------------------------
                case {'E4','E6','E8'} % 3D
                %--------------------------
                    
                % Discrete derivative using tri-diagonal scheme
                alpha = fornbergAlgorithm(1,0,[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                Ax = speye(nx);
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,[],alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = 0.5*kron(speye(nz),kron((Ax\Bx),speye(ny)));

                % Build the A and B Matrices:
                Ay = speye(ny);
                diags = repmat(alpha',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_y
                    case 0, [Ay,By,wy] = conservBoundarySchemes(Ay,By,obj.scheme,obj.nbs_id);
                    case 1, [Ay,By,wy] = set_periodicBCs(Ay,By,[],alpha);
                end
                % disp(full(Ay)); disp(full(By)); % Debug

                % Build Dy-operator
                obj.Dy = 0.5*kron(speye(nz),kron(speye(nx),(Ay\By)));

                % Build the A and B Matrices:
                Az = speye(nz);
                diags = repmat(alpha',[nz,1]);
                Bz = spdiags(diags,[-p:p],nz,nz); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_z
                    case 0, [Az,Bz,wz] = conservBoundarySchemes(Az,Bz,obj.scheme,obj.nbs_id);
                    case 1, [Az,Bz,wz] = set_periodicBCs(Az,Bz,[],alpha);
                end
                % disp(full(Az)); disp(full(Bz)); % Debug

                % Build Dz-operator
                obj.Dz = 0.5*kron((Az\Bz),kron(speye(nx),speye(ny)));
                
                % Integration weights
                obj.w = kron(wz,kron(wx,wy));
                obj.A = kron(Az,kron(Ax,Ay));
                    
                %--------------------------
                case {'T4','T6','T8'} % 3D
                %--------------------------
                    
                % Discrete derivative using tri-diagonal scheme
                [beta,alpha] = TaylorTableAlgorithm(1,[-s,s],[-p:p]); %#ok<NBRAK>
                
                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[nx,1]);
                Ax = spdiags(diags,[-s:s],nx,nx); %#ok<NBRAK>
                diags = repmat(alpha',[nx,1]);
                Bx = spdiags(diags,[-p:p],nx,nx); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_x
                    case 0, [Ax,Bx,wx] = conservBoundarySchemes(Ax,Bx,obj.scheme,obj.nbs_id);
                    case 1, [Ax,Bx,wx] = set_periodicBCs(Ax,Bx,beta,alpha);
                end
                % disp(full(Ax)); disp(full(Bx)); % Debug

                % Build Dx-operator
                obj.Dx = 0.5*kron(speye(nz),kron((Ax\Bx),speye(ny)));

                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[ny,1]);
                Ay = spdiags(diags,[-s:s],ny,ny); %#ok<NBRAK>
                diags = repmat(alpha',[ny,1]);
                By = spdiags(diags,[-p:p],ny,ny); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_y
                    case 0, [Ay,By,wy] = conservBoundarySchemes(Ay,By,obj.scheme,obj.nbs_id);
                    case 1, [Ay,By,wy] = set_periodicBCs(Ay,By,beta,alpha);
                end
                % disp(full(Ay)); disp(full(By)); % Debug

                % Build Dy-operator
                obj.Dy = 0.5*kron(speye(nz),kron(speye(nx),(Ay\By)));

                % Build the A and B Matrices:
                diags = repmat([beta(1),1,beta(2)],[nz,1]);
                Az = spdiags(diags,[-s:s],nz,nz); %#ok<NBRAK>
                diags = repmat(alpha',[nz,1]);
                Bz = spdiags(diags,[-p:p],nz,nz); %#ok<NBRAK>
                % Introduce the conservative boundary schemes of Brady & Livescu [1]
                switch obj.periodic_z
                    case 0, [Az,Bz,wz] = conservBoundarySchemes(Az,Bz,obj.scheme,obj.nbs_id);
                    case 1, [Az,Bz,wz] = set_periodicBCs(Az,Bz,beta,alpha);
                end
                % disp(full(Az)); disp(full(Bz)); % Debug

                % Build Dz-operator
                obj.Dz = 0.5*kron((Az\Bz),kron(speye(nx),speye(ny)));
                
                % Integration weights
                obj.w = kron(wz,kron(wx,wy));
                obj.A = kron(Az,kron(Ax,Ay));
                    
                end % 3D schemes
                
                % Build boundary mask
                idx = reshape(1:nx*ny*nz,[ny,nx,nz]);
                obj.index_L = idx(:,1,:); obj.index_R = idx(:,nx,:);
                obj.index_D = idx(1,:,:); obj.index_U = idx(ny,:,:);
                obj.index_B = idx(:,:,1); obj.index_F = idx(:,:,nz);
        end
    end
    
end % Public Methods
end % compactSchemes object ;)

function c = fornbergAlgorithm(k,xbar,x)
% Compute coefficients for finite difference approximation for the
% derivative of order k at xbar based on grid values at points in x.
%
% This function returns a row vector c of dimension 1 by n, where n=length(x),
% containing coefficients to approximate u^{(k)}(xbar),
% the k'th derivative of u evaluated at xbar,  based on n values
% of u at x(1), x(2), ... x(n).
%
% If U is a column vector containing u(x) at these n points, then
% c*U will give the approximation to u^{(k)}(xbar).
%
% Note for k=0 this can be used to evaluate the interpolating polynomial
% itself.
%
% Requires length(x) > k.
% Usually the elements x(i) are monotonically increasing
% and x(1) <= xbar <= x(n), but neither condition is required.
% The x values need not be equally spaced but must be distinct.
%
% This program should give the same results as fdcoeffV.m, but for large
% values of n is much more stable numerically.
%
% Based on the program "weights" in
%   B. Fornberg, "Calculation of weights in finite difference formulas",
%   SIAM Review 40 (1998), pp. 685-691.
%
% Note: Forberg's algorithm can be used to simultaneously compute the
% coefficients for derivatives of order 0, 1, ..., m where m <= n-1.
% This gives a coefficient matrix C(1:n,1:m) whose k'th column gives
% the coefficients for the k'th derivative.
%
% In this version we set m=k and only compute the coefficients for
% derivatives of order up to order k, and then return only the k'th column
% of the resulting C matrix (converted to a row vector).
% This routine is then compatible with fdcoeffV.
% It can be easily modified to return the whole array if desired.
n = length(x);
if k >= n
    error('*** length(x) must be larger than k')
end

m = k;   % change to m=n-1 if you want to compute coefficients for all
% possible derivatives.  Then modify to output all of C.
c1 = 1;
c4 = x(1) - xbar;
C = zeros(n-1,m+1);
C(1,1) = 1;
for i=1:n-1
    i1 = i+1;
    mn = min(i,m);
    c2 = 1;
    c5 = c4;
    c4 = x(i1) - xbar;
    for j=0:i-1
        j1 = j+1;
        c3 = x(i1) - x(j1);
        c2 = c2*c3;
        if j==i-1
            for s=mn:-1:1
                s1 = s+1;
                C(i1,s1) = c1*(s*C(i1-1,s1-1) - c5*C(i1-1,s1))/c2;
            end
            C(i1,1) = -c1*c5*C(i1-1,1)/c2;
        end
        for s=mn:-1:1
            s1 = s+1;
            C(j1,s1) = (c4*C(j1,s1) - s*C(j1,s1-1))/c3;
        end
        C(j1,1) = c4*C(j1,1)/c3;
    end
    c1 = c2;
end
c = C(:,end); % last column of c gives desired row vector
end % Fornberg Scheme

function [beta,alpha] = TaylorTableAlgorithm(k,ax,bx)
% Number of coefs
p = numel(bx);
n = numel(ax) + p;
% Order-of-accuracy
q = n-1;
% Build Matrix A and vector b
A_ij = zeros(n);  A_ij(:,1:numel(bx)) = 1./factorial((0:q)').*power(bx,(0:q)');
if k==1, A_ij(2:n,(numel(bx)+1):n) = 1./factorial((0:q-1)').*power(ax,(0:q-1)'); end
if k==2, A_ij(3:n,(numel(bx)+1):n) = 1./factorial((0:q-2)').*power(ax,(0:q-2)'); end
if k==3, A_ij(4:n,(numel(bx)+1):n) = 1./factorial((0:q-3)').*power(ax,(0:q-3)'); end
b_j = zeros(n,1); b_j(k+1) = 1; % disp(A); disp(b); % Debug
% Solve system A*w = b
coefs = A_ij\b_j; alpha=coefs(1:p); beta=-coefs(p+1:end);
end % Taylor Table

function [A,B,w] = conservBoundarySchemes(A,B,scheme,nbs_id)
% Controling Parameters: scheme, nbs_id (set-ID)
% for scheme = 'E4'; nbs_id = 1:101  sets;
% for scheme = 'E6'; nbs_id = 1:16   sets;
% for scheme = 'E8'; nbs_id = 1:3    sets;
% for scheme = 'T4'; nbs_id = 1:1079 sets;
% for scheme = 'T6'; nbs_id = 1:16   sets;
% for scheme = 'T8'; nbs_id = 1:25   sets;

% Add current path to the search path
if count(py.sys.path,'') == 0
    insert(py.sys.path,int32(0),'');
end

% Import my python3 module:
mod = py.importlib.import_module('compactSchemes'); disp(mod);

% Reload my module (for sanity: In case we do any modif in python)
mod = py.importlib.reload(mod);

% Evaluate my function:
dict = mod.get_alphas(['./BradyLivescu2019/',scheme,'/',scheme,'.db'],nbs_id); clear mod;

% Convert python dictionary to struct
alphas = struct(dict); disp(alphas); clear dict;

% From struct to single variables
v = fieldnames(alphas);
for i=1:numel(v)
    eval(sprintf('%s = alphas.%s;',v{i},v{i}));
end; clear alphas i v;

% Evaluate the variables
run(['./BradyLivescu2019/',scheme,'/',scheme,'.m'])

% capture the size of A and B
[N,M] = size(A); if N~=M, error('A should be a square mat.'); end
[N,M] = size(B); if N~=M, error('B should be a square mat.'); end

% Build weights
w=ones(N,1);

% Introduce Boundary schemes into mat A:
explicit=false;
try A( 1 , 2 )=beta_0_p1;  %#ok<ALIGN>
    A( 2 , 1 )=beta_1_m1;  A( 2 , 3 )=beta_1_p1;
    A( 3 , 2 )=beta_2_m1;  A( 3 , 4 )=beta_2_p1; catch, explicit=true; end
try A( 4 , 3 )=beta_3_m1;  A( 4 , 5 )=beta_3_p1; catch, disp('not T6'); end
try A( 5 , 4 )=beta_4_m1;  A( 5 , 6 )=beta_4_p1;
    A( 6 , 5 )=beta_5_m1;  A( 6 , 7 )=beta_5_p1; catch, disp('not T8'); end

try A( N ,N-1)=beta_0_p1;  %#ok<ALIGN>
    A(N-1,N-2)=beta_1_p1;  A(N-1, N )=beta_1_m1;
    A(N-2,N-3)=beta_2_p1;  A(N-2,N-1)=beta_2_m1; catch, disp('not T4'); end
try A(N-3,N-4)=beta_3_p1;  A(N-3,N-2)=beta_3_m1; catch, disp('not T6'); end
try A(N-4,N-5)=beta_4_p1;  A(N-4,N-3)=beta_4_m1;
    A(N-5,N-6)=beta_5_p1;  A(N-5,N-4)=beta_5_m1; catch, disp('not T8'); end

% Introduce Boundary schemes into mat B:
if explicit
    try
        B( 1 , 1 )= alpha_0_0;  B( N , N )=-alpha_0_0;
        B( 1 , 2 )= alpha_0_1;  B( N ,N-1)=-alpha_0_1;
        B( 1 , 3 )= alpha_0_2;  B( N ,N-2)=-alpha_0_2;
        B( 1 , 4 )= alpha_0_3;  B( N ,N-3)=-alpha_0_3;
        B( 1 , 5 )= alpha_0_4;  B( N ,N-4)=-alpha_0_4;
        B( 2 , 1 )= alpha_1_0;  B(N-1, N )=-alpha_1_0;
        B( 2 , 2 )= alpha_1_1;  B(N-1,N-1)=-alpha_1_1;
        B( 2 , 3 )= alpha_1_2;  B(N-1,N-2)=-alpha_1_2;
        B( 2 , 4 )= alpha_1_3;  B(N-1,N-3)=-alpha_1_3;
        B( 2 , 5 )= alpha_1_4;  B(N-1,N-4)=-alpha_1_4;
        B( 3 , 1 )= alpha_2_0;  B(N-2, N )=-alpha_2_0;
        B( 3 , 2 )= alpha_2_1;  B(N-2,N-1)=-alpha_2_1;
        B( 3 , 3 )= alpha_2_2;  B(N-2,N-2)=-alpha_2_2;
        B( 3 , 4 )= alpha_2_3;  B(N-2,N-3)=-alpha_2_3;
        B( 3 , 5 )= alpha_2_4;  B(N-2,N-4)=-alpha_2_4;
        w(1) = w0;              w( N ) = w0;
        w(2) = w1;              w(N-1) = w1;
        w(3) = w2;              w(N-2) = w2;
    catch, disp('not E4');
    end
    try
        B( 1 , 6 )= alpha_0_5;  B( N ,N-5)=-alpha_0_5;
        B( 1 , 7 )= alpha_0_6;  B( N ,N-6)=-alpha_0_6;
        B( 2 , 6 )= alpha_1_5;  B(N-1,N-5)=-alpha_1_5;
        B( 2 , 7 )= alpha_1_6;  B(N-1,N-6)=-alpha_1_6;
        B( 3 , 6 )= alpha_2_5;  B(N-2,N-5)=-alpha_2_5;
        B( 3 , 7 )= alpha_2_6;  B(N-2,N-6)=-alpha_2_6;
        B( 4 , 1 )= alpha_3_0;  B(N-3, N )=-alpha_3_0;
        B( 4 , 2 )= alpha_3_1;  B(N-3,N-1)=-alpha_3_1;
        B( 4 , 3 )= alpha_3_2;  B(N-3,N-2)=-alpha_3_2;
        B( 4 , 4 )= alpha_3_3;  B(N-3,N-3)=-alpha_3_3;
        B( 4 , 5 )= alpha_3_4;  B(N-3,N-4)=-alpha_3_4;
        B( 4 , 6 )= alpha_3_5;  B(N-3,N-5)=-alpha_3_5;
        B( 4 , 7 )= alpha_3_6;  B(N-3,N-6)=-alpha_3_6;
        B( 4 , 8 )= alpha_3_7;  B(N-3,N-7)=-alpha_3_7;
        B( 5 , 1 )= alpha_4_0;  B(N-4, N )=-alpha_4_0;
        B( 5 , 2 )= alpha_4_1;  B(N-4,N-1)=-alpha_4_1;
        B( 5 , 3 )= alpha_4_2;  B(N-4,N-2)=-alpha_4_2;
        B( 5 , 4 )= alpha_4_3;  B(N-4,N-3)=-alpha_4_3;
        B( 5 , 5 )= alpha_4_4;  B(N-4,N-4)=-alpha_4_4;
        B( 5 , 6 )= alpha_4_5;  B(N-4,N-5)=-alpha_4_5;
        B( 5 , 7 )= alpha_4_6;  B(N-4,N-6)=-alpha_4_6;
        B( 5 , 8 )= alpha_4_7;  B(N-4,N-7)=-alpha_4_7;
        w(4) = w3;              w(N-3) = w3;
        w(5) = w4;              w(N-4) = w4;
    catch, disp('not E6');
    end
    try
        B( 1 , 8 )= alpha_0_7;  B( N ,N-7)=-alpha_0_7;
        B( 1 , 9 )= alpha_0_8;  B( N ,N-8)=-alpha_0_8;
        B( 2 , 8 )= alpha_1_7;  B(N-1,N-7)=-alpha_1_7;
        B( 2 , 9 )= alpha_1_8;  B(N-1,N-8)=-alpha_1_8;
        B( 3 , 8 )= alpha_2_7;  B(N-2,N-7)=-alpha_2_7;
        B( 3 , 9 )= alpha_2_8;  B(N-2,N-8)=-alpha_2_8;
        B( 4 , 9 )= alpha_3_8;  B(N-3,N-8)=-alpha_3_8;
        B( 5 , 9 )= alpha_4_8;  B(N-4,N-8)=-alpha_4_8;
        B( 6 , 1 )= alpha_5_0;  B(N-5, N )=-alpha_5_0;
        B( 6 , 2 )= alpha_5_1;  B(N-5,N-1)=-alpha_5_1;
        B( 6 , 3 )= alpha_5_2;  B(N-5,N-2)=-alpha_5_2;
        B( 6 , 4 )= alpha_5_3;  B(N-5,N-3)=-alpha_5_3;
        B( 6 , 5 )= alpha_5_4;  B(N-5,N-4)=-alpha_5_4;
        B( 6 , 6 )= alpha_5_5;  B(N-5,N-5)=-alpha_5_5;
        B( 6 , 7 )= alpha_5_6;  B(N-5,N-6)=-alpha_5_6;
        B( 6 , 8 )= alpha_5_7;  B(N-5,N-7)=-alpha_5_7;
        B( 6 , 9 )= alpha_5_8;  B(N-5,N-8)=-alpha_5_8;
        B( 6 , 10)= alpha_5_9;  B(N-5,N-9)=-alpha_5_9;
        B( 7 , 1 )= alpha_6_0;  B(N-6, N )=-alpha_6_0;
        B( 7 , 2 )= alpha_6_1;  B(N-6,N-1)=-alpha_6_1;
        B( 7 , 3 )= alpha_6_2;  B(N-6,N-2)=-alpha_6_2;
        B( 7 , 4 )= alpha_6_3;  B(N-6,N-3)=-alpha_6_3;
        B( 7 , 5 )= alpha_6_4;  B(N-6,N-4)=-alpha_6_4;
        B( 7 , 6 )= alpha_6_5;  B(N-6,N-5)=-alpha_6_5;
        B( 7 , 7 )= alpha_6_6;  B(N-6,N-6)=-alpha_6_6;
        B( 7 , 8 )= alpha_6_7;  B(N-6,N-7)=-alpha_6_7;
        B( 7 , 9 )= alpha_6_8;  B(N-6,N-8)=-alpha_6_8;
        B( 7 , 10)= alpha_6_9;  B(N-6,N-9)=-alpha_6_9;
        B( 7 , 11)= alpha_6_10; B(N-6,N-10)=-alpha_6_10;
        w(6) = w5;              w(N-5) = w5;
        w(7) = w6;              w(N-6) = w6;
    catch, disp('not E8');
    end
else
    try
        B( 1 , 1 )= alpha_0_0;  B( N , N )=-alpha_0_0;
        B( 1 , 2 )= alpha_0_1;  B( N ,N-1)=-alpha_0_1;
        B( 1 , 3 )= alpha_0_2;  B( N ,N-2)=-alpha_0_2;
        B( 1 , 4 )= alpha_0_3;  B( N ,N-3)=-alpha_0_3;
        B( 2 , 1 )= alpha_1_0;  B(N-1, N )=-alpha_1_0;
        B( 2 , 3 )= alpha_1_2;  B(N-1,N-2)=-alpha_1_2;
        B( 2 , 4 )= alpha_1_3;  B(N-1,N-3)=-alpha_1_3;
        B( 3 , 2 )= alpha_2_1;  B(N-2,N-1)=-alpha_2_1;
        B( 3 , 4 )= alpha_2_3;  B(N-2,N-3)=-alpha_2_3;
        w(1) = w0;              w( N ) = w0;
        w(2) = w1;              w(N-1) = w1;
        w(3) = w2;              w(N-2) = w2;
    catch, disp('not T4');
    end
    try
        B( 1 , 5 )= alpha_0_4;  B( N ,N-4)=-alpha_0_4;
        B( 1 , 6 )= alpha_0_5;  B( N ,N-5)=-alpha_0_5;
        B( 2 , 5 )= alpha_1_4;  B(N-1,N-4)=-alpha_1_4;
        B( 2 , 6 )= alpha_1_5;  B(N-1,N-5)=-alpha_1_5;
        B( 3 , 1 )= alpha_2_0;  B(N-2, N )=-alpha_2_0;
        B( 3 , 5 )= alpha_2_4;  B(N-2,N-4)=-alpha_2_4;
        B( 3 , 6 )= alpha_2_5;  B(N-2,N-5)=-alpha_2_5;
        B( 4 , 1 )= alpha_3_0;  B(N-3, N )=-alpha_3_0;
        B( 4 , 2 )= alpha_3_1;  B(N-3,N-1)=-alpha_3_1;
        B( 4 , 3 )= alpha_3_2;  B(N-3,N-2)=-alpha_3_2;
        B( 4 , 5 )= alpha_3_4;  B(N-3,N-4)=-alpha_3_4;
        B( 4 , 6 )= alpha_3_5;  B(N-3,N-5)=-alpha_3_5;
        w(4) = w3;              w(N-3) = w3;
    catch, disp('not T6');
    end
    try
        B( 1 , 7 )= alpha_0_6;  B( N ,N-6)=-alpha_0_6;
        B( 1 , 8 )= alpha_0_7;  B( N ,N-7)=-alpha_0_7;
        B( 2 , 7 )= alpha_1_6;  B(N-1,N-6)=-alpha_1_6;
        B( 2 , 8 )= alpha_1_7;  B(N-1,N-7)=-alpha_1_7;
        B( 3 , 7 )= alpha_2_6;  B(N-2,N-6)=-alpha_2_6;
        B( 3 , 8 )= alpha_2_7;  B(N-2,N-7)=-alpha_2_7;
        B( 4 , 7 )= alpha_3_6;  B(N-3,N-6)=-alpha_3_6;
        B( 4 , 8 )= alpha_3_7;  B(N-3,N-7)=-alpha_3_7;
        B( 5 , 1 )= alpha_4_0;  B(N-4, N )=-alpha_4_0;
        B( 5 , 2 )= alpha_4_1;  B(N-4,N-1)=-alpha_4_1;
        B( 5 , 3 )= alpha_4_2;  B(N-4,N-2)=-alpha_4_2;
        B( 5 , 4 )= alpha_4_3;  B(N-4,N-3)=-alpha_4_3;
        B( 5 , 6 )= alpha_4_5;  B(N-4,N-5)=-alpha_4_5;
        B( 5 , 7 )= alpha_4_6;  B(N-4,N-6)=-alpha_4_6;
        B( 5 , 8 )= alpha_4_7;  B(N-4,N-7)=-alpha_4_7;
        B( 5 , 9 )= alpha_4_8;  B(N-4,N-8)=-alpha_4_8;
        B( 6 , 1 )= alpha_5_0;  B(N-5, N )=-alpha_5_0;
        B( 6 , 2 )= alpha_5_1;  B(N-5,N-1)=-alpha_5_1;
        B( 6 , 3 )= alpha_5_2;  B(N-5,N-2)=-alpha_5_2;
        B( 6 , 4 )= alpha_5_3;  B(N-5,N-3)=-alpha_5_3;
        B( 6 , 5 )= alpha_5_4;  B(N-5,N-4)=-alpha_5_4;
        B( 6 , 7 )= alpha_5_6;  B(N-5,N-6)=-alpha_5_6;
        B( 6 , 8 )= alpha_5_7;  B(N-5,N-7)=-alpha_5_7;
        B( 6 , 9 )= alpha_5_8;  B(N-5,N-8)=-alpha_5_8;
        w(5) = w4;              w(N-4) = w4;
        w(6) = w5;              w(N-5) = w5;
    catch, disp('not T8');
    end
end
end % conservative boundaries schemes

function [A,B,w] = set_periodicBCs(A,B,beta_C,alpha_C)
% Impose boundary conditions

% Check inputs is a tri-diagonal | explicit central scheme
[ny,nx] = size(A); if nx~=ny; error('matrix A is not square!'); end
if numel(beta_C)~=2, disp('central scheme is not tridiagonal'); end
if rem(numel(alpha_C),2)~=1, error('scheme is uneven'); end
N = numel(alpha_C); n = (N-1)/2;
% Set A coefs:
if numel(beta_C)==2
    A( 1 ,nx ) = beta_C(2);
    A( ny, 1 ) = beta_C(1);
end
% Set B coefs:
switch n
    case 1
        B(1,nx) = alpha_C( 1 );
        B(ny,1) = alpha_C( N );
    case 2
        B(1,nx-1:nx) = alpha_C( 1:2 );
        B(2,   nx  ) = alpha_C(  1  );
        B(ny-1, 1  ) = alpha_C(  N  );
        B( ny ,1:2 ) = alpha_C(N-1:N);
    case 3
        B(1,nx-2:nx) = alpha_C( 1:3 );
        B(2,nx-1:nx) = alpha_C( 1:2 );
        B(3,   nx  ) = alpha_C(  1  );
        B(ny-2, 1  ) = alpha_C(  N  );
        B(ny-1,1:2 ) = alpha_C(N-1:N);
        B( ny ,1:3 ) = alpha_C(N-2:N);
    case 4
        B(1,nx-3:nx) = alpha_C( 1:4 );
        B(2,nx-2:nx) = alpha_C( 1:3 );
        B(3,nx-1:nx) = alpha_C( 1:2 );
        B(4,   nx  ) = alpha_C(  1  );
        B(ny-3, 1  ) = alpha_C(  N  );
        B(ny-2,1:2 ) = alpha_C(N-1:N);
        B(ny-1,1:3 ) = alpha_C(N-2:N);
        B( ny ,1:4 ) = alpha_C(N-3:N);
    otherwise, error('periodic BCs not set for this case!');
end
% Weights
w = ones(nx,1);
end