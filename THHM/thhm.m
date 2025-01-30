%
% Two-harmonic homotopy method
%
% Input:
%   * (X,Lambda): results of a continuation of the main NFR (Fourier
%   coefficients, frequency)
%   * m,l: incommensurate natural coefficients of the m:l resonance
%   * phir: estimated resonant phase lag
%   * cstrType2: second type of constraint
%     * 1: total forcing amplitude constraint
%     * 2: harmonic amplitude constraint
%     * other: total amplitude constraint
%   * f: forcing amplitude 
%   * fext: forcing vector
%   * Ndofs: number of degrees of freedom
%   * Nh: number of harmonics
%   * A0, A1, A2: HB linear structural matrices
%   * NL: nonlinearity
%   * Gamma: discrete Fourier transform matrix
%   * NFT: number of points in the Fourier transform
%   * SC: simpleContinuation object
%   * fBounds: bounds for the force for the THHM
%   * indCond: condensed harmonic coefficients (optional)
%   * uextr: localization vector for the dof of interest
% Output:
%   * (XTH,LambdaTH): THHM results:
%     * XTH = [zTH ; OmegaTH ; fl] (Fourier coefficients, frequency,
%     forcing on the lth harmonic)
%     * LambdaTH = fm (forcing on the mth harmonic)
%   * (bm,bl): forcing vectors of the mth and lth harmonics, respectively
%   * (xI, wI, fI): initial guesses for an isola
%     * xI: Fourier coefficients
%     * wI: frequency
%     * fI: forcing amplitude
%   
function [XTH,LambdaTH,bm,bl,xI,wI,fI] = thhm(X,Lambda,m,l,phir,cstrType2,f,fext,Ndofs,Nh,A0,A1,A2,NL,Gamma,NFT,SC,fBounds,indCond,uextr)
  
  if nargin < 18
    fBounds = [-0.1,2*f];
  end
  if nargin < 19
    indCond = [];
    uextr = fext;
  end

  % Find phase resonance
  phi = atan2(X(3,:),X(2,:));
  [~,ind] = min(abs(phi+pi/2));
  
  % Initial point
  T = harmonicTransformationMatrix(m,Ndofs,Nh,Nh);
  X0 = T*X(:,ind);
  W0 = Lambda(ind)/m;
  fm = f;
  fl = 0;
  
  % Phase lags/shifts
  if isempty(phir)
    if mod(m,2) == 0 || mod(l,2) == 0
      phir = 3*pi/(4*l);
    else
      phir = pi/2;
    end
  end
  phil = l*(phir-pi/2)/m;
  
  % Forcing vectors
  em = double(1:(2*Nh+1) == 2*m).';
  el1 = double(1:(2*Nh+1) == 2*l).';
  el2 = double(1:(2*Nh+1) == 2*l+1).';
  bm = kron(em,fext);
  bl = kron(el1*cos(phil)+el2*sin(phil),fext);
  
  % THHM constraints
  cstrVec1 = kron(em,uextr);
  switch cstrType2
    case 1
      cstrVec2 = [1,1];
      cstrVal2 = cstrVec2*[fm;fl];
    case 2
      em2 = double(1:(2*Nh+1) == 2*m+1).';
      cstrVec2 = [kron(em,uextr).' ; kron(em2,uextr).'];
      cstrVal2 = norm(cstrVec2*X0)^2;
    otherwise
      cstrVec2 = kron(eye(2*Nh+1),uextr.');
      cstrVal2 = FourierAmplitude(X0,cstrVec2);
  end
   
  
  % THHM
  [XTH,LambdaTH] = SC.continuation([X0;W0;fl],fm,...
    @ (z,lambda) THHMResidual(z(1:end-2),z(end-1),lambda,z(end),bm,bl,A0,A1,A2,NL,Gamma,NFT,cstrVec1,cstrType2,cstrVec2,cstrVal2,indCond),...
    fBounds,-1);


  if nargout > 4
    
    % Find fm=0 points
    ind = find(sign(LambdaTH(2:end))~=sign(LambdaTH(1:end-1)));
  
    xI = zeros(size(XTH,1)-2,length(ind));
    wI = zeros(1,length(ind));
    fI = zeros(1,length(ind));

    % Interpolate
    for ii = 1 : length(ind) 
      iInt = ind(ii)+[-1,0,1];
    
      xI(:,ii) = interp1(LambdaTH(iInt),XTH(1:end-2,iInt).',0,'pchip').';
      wI(ii) = interp1(LambdaTH(iInt),XTH(end-1,iInt).',0,'pchip');
      fI(ii) = interp1(LambdaTH(iInt),XTH(end,iInt),0,'pchip');
    end
  end
end