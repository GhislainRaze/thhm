% Two-degree-of-freedom example
% THHM with constant amplitude

%% Problem definition

% Nonlinear structure definition
M = eye(2);
C = [0.0244,-0.0144;-0.0144,0.0244];
K = [8,-7;-7,8];

Ndofs = length(M);
NL = @(x) cubicNonlinearity(x);

fext = [1;0];
f = 0.5;


% Harmonic balance quantities
Nh = 31;
NFT = 128;
[A0,A1,A2] = HBLinearMatrices(M,C,K,Nh);
Gamma = HBGamma(Nh,NFT,Ndofs);

%% Continuation object
SC = simpleContinuation();

% Numerical parameters
SC.tol = 1e-9;
SC.stepMax = 4e2;
SC.hMax = 0.1;


%% NFR
e2 = double(1:(2*Nh+1) == 2).';
bext = f*kron(e2,fext);


z0 = A0\bext;
lambda0 = 0;
[X,Lambda] = SC.continuation(z0,lambda0,...
  @ (z,lambda) HBResidual(z,lambda,bext,A0,A1,A2,NL,Gamma,NFT),[-0.1,2]);



a = FourierAmplitude(X,kron(eye(2*Nh+1),[1,0]));
figure
hold on
plot(Lambda,a,'-k')



%% THHM
m = 1;
l = 2;
% cstrType2 = 2;    % Constant harmonic amplitude
cstrType2 = 3;    % Constant total amplitude
[XTH,LambdaTH,bm,bl,xI,wI,fI] = thhm(X,Lambda,m,l,[],cstrType2,f,fext,Ndofs,Nh,A0,A1,A2,NL,Gamma,NFT,SC,[-1,3]*f);

%% Isola continuation
if ~isempty(fI)
  fl = fI(1);
  bextI = fl*bl;
  SC.detectIsola = true;
  [XI,LambdaI] = SC.continuation(xI(:,1),wI(1),...
    @ (z,lambda) HBResidual(z,lambda,bextI,A0,A1,A2,NL,Gamma,NFT),[-0.1,2]*wI(1));
  
  
  aI = FourierAmplitude(XI,kron(eye(2*Nh+1),[1,0]));
  plot(l*LambdaI,aI,'k')
end


box on
xlabel('Excitation frequency (-)')
ylabel('Amplitude (-)')


aTH = FourierAmplitude(XTH(1:end-2,:),kron(eye(2*Nh+1),[1,0]));
figure
plot(LambdaTH,aTH)
set(gca,'xDir','reverse')
xlabel('f_m (-)')
ylabel('Amplitude (-)')



function [fnl,dfnl] = cubicNonlinearity(x)
  fnl = 0*x;
  fnl(1:2:end) = 0.5*x(1:2:end).^3;
  dfnl = zeros(length(x),length(x));
  dfnl(1:2:end,1:2:end) = diag(1.5*x(1:2:end).^2);
end