% Beam with clearance contact example

%% Problem definition

% Nonlinear structure definition
lx   = 500e-3; % Length, m
ly   = 20e-3;  % Width, m
lz   = 1e-3;   % Thickness, m
E    = 210e9;  % Young's modulus, Pa
rho  = 7800;   % Density, kg/m^3
Ne   = 10;     % Number of elements (-)
[M,K] = FEMBeam(lx,ly,lz,E,rho,Ne);

% Boundary conditions
M = M(4:end,4:end);
K = K(4:end,4:end);

% Resonance frequencies
[Phi,Om2] = eig(full(K),full(M));
xi = 5e-3;
c = [1 , Om2(1,1) ; 1 , Om2(2,2)]\[2*xi*sqrt(Om2(1,1));2*xi*sqrt(Om2(2,2))];
C = c(1)*M+c(2)*K;
Z2O = diag(diag(Phi.'*C*Phi));

Ndofs = length(M);
Ndofsr = 1;           % Reduced number of dofs

% Nonlinearity
I = ly*lz^3/12;       % Cross-section intertia (m^4)
kst = 3*E*I/lx^3;     % Static stiffnes at the beam tip (N/m)
knl = 99*kst;         % Nonlinear element stifness (N/m)
x0 = 1e-3;            % Gap (m)
l0 = 10e-6;           % Regularization legth (m)

NL = @(x) contactNL(x,knl,x0,l0);

% External forcing
fext = zeros(length(K),1);
fext(end-1) = 1;
f = 4.5e-4;

% Harmonic balance quantities
Nh = 51;
NFT = 512;
[A0,A1,A2] = HBLinearMatrices(M,C,K,Nh);
Gamma = HBGamma(Nh,NFT,1);
indNL = (length(M)-1):length(M):length(A0);
indCond = setdiff(1:1:length(A0),indNL);

%% Continuation object
SC = simpleContinuation();

% Numerical parameters
SC.tol = 1e-9;
SC.stepMax = 3e2;
SC.psi = 0.01;  % Scaling factor for the frequency, the continuation will take steps in lambda ~h/SC.psi
SC.hMin = 1e-9;
SC.hMax = 1e-2;

%% NFR
e2 = double(1:(2*Nh+1) == 2).';
bext = f*kron(e2,fext);

z0 = A0\bext;
z0 = z0(indNL);
lambda0 = 0;
[X,Lambda] = SC.continuation(z0,lambda0,...
  @ (z,lambda) HBResidual(z,lambda,bext,A0,A1,A2,NL,Gamma,NFT,indCond),[-0.1,30]);


a = FourierAmplitude(X,eye(size(X,1)));
figure
hold on
plot(Lambda,a,'k')


%% THHM
SC.psi = 1;
cstrType2 = 3;    % Constant total amplitude

m = 1;
l1 = 3;
[XTH,LambdaTH,bm,bl,xI,wI,fI] = thhm(X,Lambda,m,l1,[],cstrType2,f,fext,Ndofsr,Nh,A0,A1,A2,NL,Gamma,NFT,SC,[-0.1*f,3*f],indCond,1);

m = 1;
l2 = 5;
phir = -pi/2;
[XTH2,LambdaTH2,bm2,bl2,xI2,wI2,fI2] = thhm(X,Lambda,m,l2,phir,cstrType2,f,fext,Ndofsr,Nh,A0,A1,A2,NL,Gamma,NFT,SC,[-0.1*f,5*f],indCond,1);


%% Isola continuation
SC.psi = 0.01;
SC.detectIsola = true;
SC.dxIsola = 1e-3;
SC.hMax = 1e-4;

if ~isempty(fI)
  fl = fI(1);
  bextI = fl*bl;
  [XI,LambdaI] = SC.continuation(xI(:,1),wI(1),...
    @ (z,lambda) HBResidual(z,lambda,bextI,A0,A1,A2,NL,Gamma,NFT,indCond),[-0.1,2]*wI(1));

  aI = FourierAmplitude(XI,eye(size(XI,1)));
  plot(l1*LambdaI,aI,'k')
end

if ~isempty(fI2)
  fl2 = fI2(1);
  bextI2 = fl2*bl2;
  [XI2,LambdaI2] = SC.continuation(xI2(:,1),wI2(1),...
    @ (z,lambda) HBResidual(z,lambda,bextI2,A0,A1,A2,NL,Gamma,NFT,indCond),[-0.1,2]*wI2(1),-1);


  aI2 = FourierAmplitude(XI2,eye(size(XI2,1)));
  plot(l2*LambdaI2,aI2,'k')
end

box on
xlabel('Excitation frequency (rad/s)')
ylabel('Amplitude (m)')



function [fnl,dfnl] = contactNL(x,knl,x0,l0)
  xi = (x-x0)/l0;
  fnl = knl*l0/2*(log(cosh(xi)) + xi + log(2));
  dfnl = diag(knl/2*(1+tanh(xi)));
end