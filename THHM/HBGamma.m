%
% Discrete Fourier transform matrix
%
% Input:
%   * Nh: number of harmonics
%   * NFT: number of sampling points
%   * Ndofs: number of degrees of freedom
% 
% Output:
%   * Gamma: discrete Fourier transform matrix
%
function Gamma = HBGamma(Nh,NFT,Ndofs)
  
  % Dimensionless time vector
  tau = linspace(0,2*pi,NFT+1).';
  tau = tau(1:end-1);

  % Harmonic functions
  Q = [1/sqrt(2)*ones(size(tau)) , kron(sin(tau*(1:1:Nh)),[1,0])+kron(cos(tau*(1:1:Nh)),[0,1])];

  % Discrete Fourier transform matrix
  Gamma = kron(Q,speye(Ndofs));
end