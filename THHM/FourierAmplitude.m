% 
% Amplitude of a Fourier series
% 
% Input:
%   * z: vector of harmonic coefficients, or matrix whose columns are
%   vectos of harmonic coefficients
%   * v: localization vector, such that v*z yields the harmonic
%   coefficients of a degree of freedom of interest
%
% Output:
%   * a: amplitude of the Fourier series (length(a) = size(z,2))
%   * da: derivative of the amplitude with respect to the coefficients
%
function [a,da] = FourierAmplitude(z,v)
  
  
  u = v*z;
  Nh = (size(u,1)-1)/2;
  a = zeros(size(u,2),1);
  if nargout > 1
    da = zeros(size(u,2),length(z));
  end
  d = 1i*(1:1:Nh).';
  for ii = 1 : size(u,2)
    poly = [-flipud(d).*flipud(u(3:2:end,ii)-1i*u(2:2:end,ii))/2;0;d.*(u(3:2:end,ii)+1i*u(2:2:end,ii))/2];
    tau = angle(roots(poly));
    xu = ones(size(tau))*u(1,ii)/sqrt(2) + sin(tau*(1:1:Nh))*u(2:2:end,ii) + cos(tau*(1:1:Nh))*u(3:2:end,ii); 
    [a(ii),ind] = max(abs(xu));
    if nargout > 1
      da(ii,:) = sign(xu(ind))*[1/sqrt(2) , kron(sin(tau(ind)*(1:1:Nh)),[1,0]) + kron(cos(tau(ind)*(1:1:Nh)),[0,1])]*v;
    end
  end
end