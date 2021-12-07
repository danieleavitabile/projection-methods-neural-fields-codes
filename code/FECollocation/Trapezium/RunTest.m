function RunTest(rootPath,resultsPath)

  %% Adding directory
  close all; addpath(rootPath);

  disp('***********************************************************************');
  disp('Running tests for FE Collocation scheme with trapezium quadrature rule');
  disp('***********************************************************************');

  %% List of problems to solve
  pList = {'P1','P2','P3','P4','P5','P6'};

  %% Number of gridpoints in each problem
  nVals = [2:10:100 200:100:1600];
  % nVals = 2:10:100;

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
      nx = nVals(m); x = linspace(-1,1,nx)'; hx = 2/(nx-1); rho = hx*[0.5; ones(nx-2,1); 0.5];

      % Integration weights and linear operators
      W = zeros(nx,nx);
      for i = 1:nx
	for j = 1:nx
	  W(i,j) = p.wFun(x(i),x(j))*rho(j);
	end
      end

      % Right-hand side function handle
      N = @(t,u) -u + W*p.f(u) + p.xi(x,t);

      % Time step
      u0 = p.uAna(x,0); 
      tspan = [0 3];
      [t,U] = ode45(N,tspan,u0);

      % Comptue error
      [X,T] = meshgrid(x,t);
      eVec(m) = max(max(abs(U-p.uAna(X,T))));

    end

    figure(fig), hold on; plot(nVals,eVec,'.-','DisplayName',pFlag); hold off; drawnow;

  end

  set(gca,'XScale','log','YScale','log'); box on; legend('NumColumns',2);
  title('Finite Elements Colloc. (Trapezium)');

  savefigure

end
