%
% Residual for the two-harmonic homotopy method
%
% Input:
%   * z: vector of harmonic coefficients
%   * w: frequency
%   * fm, fl: harmonic forcing amplitudes
%   * bm, bl: harmonic forcing vectors
%   * A0, A1, A2: matrices such that A(w) = A0 + w*A1 + w^2*A2
%   * NL: nonlinearities
%   * Gamma: discrete Fourier transform matrix
%   * NFT: number of samples in the Fourier transform
%   * cstrVec1
%   * cstrType2, cstrVec2:
%     * If cstrType2 = 1 (constant-force), cstrVec2 = [alpham , alphal]
%     * If cstrType2 = 2 (constant harmonic amplitude), cstrVec2 is such
%     that the mth harmonic coefficients of the degree of freedom of
%     interest are given by cstrVec2*z
%     * If cstrType2 = 3 (constant total amplitude), cstrVec2 is such that
%     the harmonic coefficients of the degree of freedom of interest are
%     given by cstrVec2*z
%     
% Output:
%   * f: residual
%   * df: Jacobian
%
function [f,df] = THHMResidual(z,w,fm,fl,bm,bl,A0,A1,A2,NL,Gamma,NFT,cstrVec1,cstrType2,cstrVec2,cstrVal2,indCond)

  % Harmonic balance residual
  if nargin < 17
    indCond = [];
    indNL = 1:1:length(z);
  else
    indNL = setdiff(1:1:length(A0),indCond);
  end

  bext = fm*bm + fl*bl;
  [fHB,dfHB] = HBResidual(z,w,bext,A0,A1,A2,NL,Gamma,NFT,indCond);

  % Frequency constraint
  c1 = cstrVec1.'*z;
  dc1 = [cstrVec1.',zeros(1,3)];

  % Forcing constraint
  switch cstrType2
    case 1
      c2 = cstrVec2*[fm;fl] - cstrVal2;
      dc2 = [zeros(1,length(z)+1) , cstrVec2];
    case 2
      u = cstrVec2*z;
      c2 = u(1)^2+u(2)^2-cstrVal2;
      dc2 = [[2*u(1),2*u(2)]*cstrVec2 , zeros(1,3)];
    otherwise
      [a,da] = FourierAmplitude(z,cstrVec2);
      % da = da*cstrVec2;
      c2 = a-cstrVal2;
      dc2 = [da , zeros(1,3)];
  end

  % Linear dynamics matrix & reduced
  if ~isempty(indCond)
    A = A0 + w*A1 + w^2*A2;
    Ai = decomposition(A(indCond,indCond));
    dfl = -bl(indNL)-A(indNL,indCond)*(Ai\bl(indCond));
    dfm = -bm(indNL)-A(indNL,indCond)*(Ai\bm(indCond));
  else
    dfl = -bl;
    dfm = -bm;
  end
  % Residual
  f = [ fHB ; c1 ; c2];
  df = [dfHB , dfl , dfm ; dc1 ; dc2];

end