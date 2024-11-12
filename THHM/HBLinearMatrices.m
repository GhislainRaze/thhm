%
% Matrices associated with the linear dynamics for the harmonic balance
% method
%
% Input:
%   * M, C, K: linear mass, damping and stiffness matrices
%   * Nh: number of harmonics
% 
% Output:
%   * A0, A1, A2: matrices such that A(w) = A0 + w*A1 + w^2*A2 for the
%   harmonic balance method
%
function [A0,A1,A2] = HBLinearMatrices(M,C,K,Nh)
  
  % Harmonic operators
  Ih = speye(2*Nh+1);
  nabla = blkdiag(0,kron(sparse(1:1:Nh,1:1:Nh,1:1:Nh,Nh,Nh),[0,-1;1,0]));
  
  % Linear dynamics matrices
  A0 = kron(Ih,K);
  A1 = kron(nabla,C);
  A2 = kron(nabla^2,M);

end