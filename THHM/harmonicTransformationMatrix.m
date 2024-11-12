%
% Harmonic transformation matrix for m-subharmonic motion
%
% Input:
%   * m: period multiplier
%   * Ndof: number of degrees of freedom
%   * Nh1: number of harmonics in the original vector
%   * Nh2: number of harmonics in the transformed vector (default: Nh1)
%
% Output:
%   * T: harmonic transformation matrix
function T = harmonicTransformationMatrix(m,Ndof,Nh1,Nh2)

  if nargin < 4
    Nh2 = Nh1;
  end

  T = sparse([1,2*m*(1:1:Nh1),2*m*(1:1:Nh1)+1],[1,2*(1:1:Nh1),2*(1:1:Nh1)+1],...
    ones(2*Nh1+1,1),2*m*Nh1+1,2*Nh1+1);
  T = kron(T(1:(2*Nh2+1),:),speye(Ndof));
end