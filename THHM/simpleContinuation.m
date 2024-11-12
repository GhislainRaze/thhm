classdef simpleContinuation

  properties
    tol = 1e-6;       % Absolute tolerance
    iterMax = 15;     % Maximum number of corrector iterations
    hMax = 1;         % Maximum predictor step
    hMin = 1e-6;      % Minimum predictor step
    stepMax = 1e3;    % Maximum number of steps
    psi = 1;          % Scaling parameter for lambda (select it such that psi*lambda = norm(z))
    display = true;   % Display info about the continuation
    detectIsola = false;
    dxIsola = 5e-2;
  end

  methods


    function [X,Lambda] = continuation(obj,z,lambda,residual,LambdaBounds,dir)

      if nargin < 6
        dir = 1;
      end

      % Initializations
      X = zeros(length(z),obj.stepMax+2);
      Lambda = zeros(1,obj.stepMax+2);

      % Correct first guess at constant frequency
      x = [z;obj.psi*lambda];
      t = [zeros(size(z));dir];
      fun = @(dx) obj.pseudoArclengthFun(dx,x,residual,t);
      dx = zeros(size(x));
      [f,df] = fun(dx);
      iter = 0;
      while norm(f) > obj.tol && iter < obj.iterMax
        dx = dx-df\f;
        [f,df] = fun(dx);
        iter = iter+1;
      end
      x = x+dx;
      X(:,1) = x(1:end-1);
      Lambda(1) = lambda;

      % Tangent predictor
      t = [df(1:end-1,:);t.']\[zeros(length(z),1);1];
      t = t/norm(t);

      % Iteration counter
      step = 0;
      h = obj.hMax;

      %
      x0 = x;

      % Continuation procedure
      while step < obj.stepMax && x(end)/obj.psi > LambdaBounds(1) && x(end)/obj.psi < LambdaBounds(2)
        % Prediction
        x = x + t*h;

        % Pseudoarclength correction
        fun = @(dx) obj.pseudoArclengthFun(dx,x,residual,t);
        dx = zeros(size(x));
        [f,df] = fun(dx);
        iter = 0;
        while norm(f) > obj.tol && iter < obj.iterMax
          dx = dx-df\f;
          [f,df] = fun(dx);
          iter = iter+1;
        end

        if norm(f) < obj.tol
          x = x+dx;

          % Save solution
          X(:,step+2) = x(1:end-1);
          Lambda(step+2) = x(end)/obj.psi;
        else
          % Return to previous point and halve the step
          x = x - t*h;
          if h~=obj.hMin
            h = max([h/2,obj.hMin]);
          else
            disp('Unable to converge.')
            break
          end
          continue
        end

        % Adaptive step
        h = h*sqrt(3/max([iter,1]));
        h = max([min([h,obj.hMax]),obj.hMin]);

        % Tangent predictor
        t = [df(1:end-1,:);t.']\[zeros(length(z),1);1];
        t = t/norm(t);
        
        if obj.display
          disp(['Step ',num2str(step,'%.5d'),': lambda = ',num2str(x(end)/obj.psi,'%.5e')])
        end

        step = step+1;

        if obj.detectIsola == true && norm(x-x0) < obj.dxIsola*norm(x0) && step*h > obj.dxIsola*norm(x0)
          break
        end
      end

      % Free unused variables
      X(:,step+2:end) = [];
      Lambda(step+2:end) = [];
    end


    function [f,df] = pseudoArclengthFun(obj,dx,x,residual,t)
      xc = x+dx;
      [r,dr] = residual(xc(1:end-1),xc(end)/obj.psi);
      f = [r;t.'*dx];
      dr(:,end) = dr(:,end)/obj.psi;
      df = [dr ; t.'];
    end
    
  end

end
