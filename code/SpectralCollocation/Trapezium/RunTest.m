
function RunTest(rootPath,resultsPath)

  %% Adding directory
  close all; addpath(rootPath);

  disp('***********************************************************************');
  disp('Running tests for Spectral Collocation scheme with trapezium quadrature rule');
  disp('***********************************************************************');

  %% List of problems to solve
  pList = {'P1','P2','P3','P4','P5','P6'};

  %% Number of gridpoints in each problem
  nVals = 2:10:300;
  %nVals = [2:10:100 200:100:1600];

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
      nx = nVals(m); x = cos(pi*(0:nx)'/nx); 

      % Integration grid, different from collocation grid
      ny = nx; y = linspace(-1,1,ny+1)'; hy = 2/ny; rho = hy*[0.5; ones(ny-1,1); 0.5];

      % Integration weights and linear operators
      W = zeros(nx+1,ny+1);
      for i = 1:nx+1
	for j = 1:ny+1
	  W(i,j) = p.wFun(x(i),y(j))*rho(j);
	end
      end

      % Right-hand side function handle
      N = @(t,u) -u + W*p.f(Lagrange(x',u',y'))' + p.xi(x,t);

      % Time step
      u0 = p.uAna(x,0); 
      tspan = [0 3];
      [t,U] = ode45(N,tspan,u0);

      % Comptue error
      [X,T] = meshgrid(x,t); 
      eVec(m) = max(max(abs(U-p.uAna(X,T))));

    end

    figure(fig), hold on; plot(nVals,eVec,'.-','DisplayName',pFlag); hold off;

  end
  % figure(fig), hold on; plot(nVals(1:4),1e2*nVals(1:4).^-2,'-','DisplayName','O(n^-2)');
  set(gca,'XScale','log','YScale','log'); box on; grid on;
  title('Spectral Colloc. (Trapezium)');

  savefigure;

end
