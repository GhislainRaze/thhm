% Duffing oscillator example
% THHM with constant harmonic amplitude

%% Problem definition

% Nonlinear structure definition
M = 1;
C = 0.05;
K = 1;

Ndofs = length(M);
NL = @(x) cubicNonlinearity(x);

fext = 1;
f = 0.5;

% Harmonic balance quantities
Nh = 31;
NFT = 128;
[A0,A1,A2] = HBLinearMatrices(M,C,K,Nh);
Gamma = HBGamma(Nh,NFT,Ndofs);

%% Continuation object
SC = simpleContinuation();

% Numerical parameters
SC.tol = 1e-6;
SC.stepMax = 4e2;
SC.hMax = 0.1;


%% NFR
e2 = double(1:(2*Nh+1) == 2).';
bext = f*kron(e2,fext);


z0 = A0\bext;
lambda0 = 0;
[X,Lambda] = SC.continuation(z0,lambda0,...
  @ (z,lambda) HBResidual(z,lambda,bext,A0,A1,A2,NL,Gamma,NFT),[-0.1,10]);



a = FourierAmplitude(X,eye(size(X,1)));
figure
hold on
plot(Lambda,a,'-k')



%% THHM

m = 1;
l = 3;
cstrType2 = 2;    % Constant harmonic amplitude
% cstrType2 = 3;    % Constant total amplitude
[XTH,LambdaTH,bm,bl,xI,wI,fI] = thhm(X,Lambda,m,l,[],cstrType2,f,fext,Ndofs,Nh,A0,A1,A2,NL,Gamma,NFT,SC);

%% Isola continuation
if ~isempty(fI)
  fl = fI(1);
  bextI = fl*bl;
  SC.detectIsola = true;
  [XI,LambdaI] = SC.continuation(xI(:,1),wI(1),...
    @ (z,lambda) HBResidual(z,lambda,bextI,A0,A1,A2,NL,Gamma,NFT),[-0.1,2]*wI(1));
  
  
  aI = FourierAmplitude(XI,eye(2*Nh+1));
  plot(l*LambdaI,aI,'k')
end

aI = FourierAmplitude(XI,eye(size(X,1)));
plot(l*LambdaI,aI,'k')
box on
xlabel('Excitation frequency (-)')
ylabel('Amplitude (-)')


aTH = FourierAmplitude(XTH(1:end-2,:),eye(size(X,1)));
figure
plot(LambdaTH,aTH)
set(gca,'xDir','reverse')
xlabel('f_m (-)')
ylabel('Amplitude (-)')



function [fnl,dfnl] = cubicNonlinearity(x)
  fnl = x.^3;
  dfnl = diag(3*x.^2);
end