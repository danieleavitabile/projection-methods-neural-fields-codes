function RunTest(rootPath,resultsPath)

  %% Adding directory
  close all; addpath(rootPath);

  disp('***********************************************************************');
  disp('Running tests for Finite Elements Galerkin scheme with Gauss quadrature');
  disp('***********************************************************************');

  %% List of problems to solve
  pList = {'P1','P2','P3','P4','P5','P6'};

  %% Number of gridpoints in each problem
  nVals = 2:10:300;

  fig = figure(1);
  %% For each problem
  for k = 1:length(pList)

    %% Load problem
    pFlag = pList{k};
    p = LoadProblem(pFlag);
    disp(['Problem ' pFlag '************************']);

    %% Error vectors and grid spacing values
    eVec = zeros(size(nVals));

    %% For every value of n
    for m = 1:length(nVals)

      disp(['n = ' num2str(nVals(m))]);

      % Spatial grid and integration weights
      nElem = nVals(m); x = linspace(-1,1,nElem+1)';  hx = 2/nElem;

      % Mass matrix
      d0 = [1/3; 2/3*ones(nElem-1,1); 1/3];
      d1 = 1/6*ones(nElem+1,1);
      M = hx*spdiags([d1 d0 d1],[-1 0 1],nElem+1,nElem+1);

      % Right-hand side function handle
      N = @(t,u) -u + M\NeuralField(t,u,p.wFun,p.f,p.xi,x);

      % Time step
      u0 = p.uAna(x,0); %u0 = M\InnerProduct(@(x) p.uAna(x,0),x);
      tspan = [0 3];
      [t,U] = ode45(N,tspan,u0);
   
      % Comptue error
      [X,T] = meshgrid(x,t(1:end));
      eVec(m) = max(max(abs(U(1:end,:)-p.uAna(X,T))));

    end


    figure(fig), hold on; plot(nVals,eVec,'.-','DisplayName',pFlag); hold off; drawnow;

  end
  set(gca,'XScale','log','YScale','log'); box on; grid on; legend;
  title('Galerkin Finite Elements (Gauss)');

  savefigure

end
