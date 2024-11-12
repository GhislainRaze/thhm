%
% Residual for the harmonic balance method
%
% Input:
%   * z: vector of harmonic coefficients
%   * w: frequency
%   * bext: external forcing vector
%   * A0, A1, A2: matrices such that A(w) = A0 + w*A1 + w^2*A2
%   * NL: nonlinearities
%   * Gamma: discrete Fourier transform matrix
%   * NFT: number of samples in the Fourier transform
%   * indCond: condensed harmonic coefficients (optional)
%
% Output:
%   * f: residual
%   * df: Jacobian
function [f,df] = HBResidual(z,w,bext,A0,A1,A2,NL,Gamma,NFT,indCond)

  % Linear dynamics matrix
  A = A0 + w*A1 + w^2*A2;
  dA = (2*w*A2+A1);

  if nargin > 9 && ~isempty(indCond)
    indNL = setdiff(1:1:length(A),indCond);
    Ai = decomposition(A(indCond,indCond));

    dbext = dA(indNL,indCond)*(Ai\bext(indCond)) - ...
      A(indNL,indCond)*(Ai\(dA(indCond,indCond)*(Ai\bext(indCond))));
    dA = dA(indNL,indNL) - dA(indNL,indCond)*(Ai\A(indCond,indNL)) - ...
      A(indNL,indCond)*(Ai\dA(indCond,indNL)) + ...
      A(indNL,indCond)*(Ai\(dA(indCond,indCond)*(Ai\A(indCond,indNL))));

    bext = bext(indNL) + A(indNL,indCond)*(Ai\bext(indCond));
    A = A(indNL,indNL) - A(indNL,indCond)*(Ai\A(indCond,indNL));
  else
    dbext = 0*bext;
  end
  
  % AFT
  x = Gamma*z;
  [fnl,dfnl] = NL(x);
  b = 2/NFT*Gamma.'*fnl;
  dbdz = 2/NFT*Gamma.'*dfnl*Gamma;

  % HB residual
  f = A*z+b-bext;
  df = [ A+dbdz , dA*z+dbext];

end