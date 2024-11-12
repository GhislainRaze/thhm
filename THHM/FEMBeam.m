% 
% Euler-Bernoulli finite element model of a beam
%
% Input:
%   * lx,ly,lz: beam dimensions
%   * E: Young modulus
%   * rho: density
%   * Ne: number of elements
%   * xe: nodal coordinates (optional)
%
% Output:
%   * M,K: mass and stiffness matrix
%   * xe: nodal coordinates
%
function [M,K,xe] = FEMBeam(lx,ly,lz,E,rho,Ne,xe)

  
  if nargin < 7
    xe = linspace(0,lx,Ne+1);
    l = lx/Ne*ones(Ne,1);
  else
    l = diff(xe);
  end
  
  
  J = zeros(1,3);
  J(1,:) = ly*[lz,0,2*(lz/2)^3/3];
  
  % Structural mass and stiffness matrices
  ndofs = 3+3*Ne;
  K = zeros(ndofs,ndofs);
  M = zeros(ndofs,ndofs);
  
  % Inertia/stiffness factors of the beam
  Ib = zeros(size(J));
  Ab = zeros(size(J));
  Ib(1,:) = rho*J(1,:);
  Ab(1,:) = E*J(1,:);  
  
  % Assembly
  for nn = 1:Ne
    % Element global degrees of freedom
    edofs = 3*nn-2:3*nn+3;

    % Element matrices
    [Ke0,Ke1,Ke2,Me0,Me1,Me2] = elementMatrices(l(nn));
    
    Ae = Ab(1,:);
    Ie = Ib(1,:);
    
    % Beam and piezoelectric patch masses and stiffnesses
    K(edofs,edofs) = K(edofs,edofs)+ Ae(1)*Ke0 + Ae(2)*Ke1 + Ae(3)*Ke2;
    M(edofs,edofs) = M(edofs,edofs)+ Ie(1)*Me0 + Ie(2)*Me1 + Ie(3)*Me2;
  end

  % Boundary conditions and sparsity
  K = sparse(K); 
  M = sparse(M);


  function [Ke0,Ke1,Ke2,Me0,Me1,Me2] = elementMatrices(l)
    
    % Elementary stiffness matrices
    Ke0 = 1/l*[1  0   0   -1  0   0;
              0   0   0   0   0   0;
              0   0   0   0   0   0;
              -1  0   0   1   0   0;
              0   0   0   0   0   0;
              0   0   0   0   0   0];

    Ke1 = 1/l*[0  0   -1  0   0   1;
              0   0   0   0   0   0;
              -1  0   0   1   0   0;
              0   0   1   0   0   -1;
              0   0   0   0   0   0;
              1   0   0   -1  0   0];

    Ke2 = 1/l^3*[0    0       0       0     0       0 ;
                 0    12     	6*l     0   	-12     6*l;
                 0    6*l     4*l^2   0     -6*l    2*l^2;
                 0    0       0       0     0       0 ;
                 0    -12     -6*l    0     12      -6*l;
                 0    6*l     2*l^2   0     -6*l    4*l^2];

    % Elementary mass matrices           
    Me0 = l/420*[140    0       0         70      0       0;
                 0      156     22*l      0       54      -13*l;
                 0      22*l    4*l^2     0       13*l    -3*l^2;
                 70     0       0         140     0       0;
                 0      54      13*l      0       156     -22*l;
                 0      -13*l   -3*l^2    0       -22*l   4*l^2];

    Me1 = 1/12*[0   6   -l    0     -6    l;
                6   0   0     6     0     0;
                -l  0   0     l     0     0;
                0   6   l     0     -6    -l;
                -6  0   0     -6    0     0;
                l   0   0     -l    0     0];

    Me2 = 1/(30*l)*[0   0     0       0   0     0;
                    0   36    3*l     0   -36   3*l;
                    0   3*l   4*l^2   0   -3*l  -l^2;
                    0   0     0       0   0     0;
                    0   -36   -3*l    0   36    -3*l;
                    0   3*l   -l^2    0   -3*l  4*l^2];

  end


end
