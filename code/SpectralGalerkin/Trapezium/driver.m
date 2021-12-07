function RunTest(rootPath,resultsPath)

  %% Adding directory
  close all; addpath(rootPath);

  disp('***********************************************************************');
  disp('Running tests for Spectral Galerkin scheme with trapezium quadrature rule');
  disp('***********************************************************************');

  %% List of problems to solve
  pList = {'P7p','P8p','P9p','P10p'}

  %% Number of gridpoints in each problem
  nVals = [2:4:100 102:100:1024];

  fig = figure(1);
  for k = 1:length(pList)

    pFlag = pList{k};
    p = LoadProblem(pFlag);
    disp(['Problem ' pFlag ' ************************']);

    %% Error vectors and grid spacing values
    eVec = zeros(size(nVals));

    %% For every value of n
    for m = 1:length(nVals)

      disp(['n = ' num2str(nVals(m))]);

      %% Spatial grid
      nx = nVals(m); Lx = pi; hx = 2*Lx/nx; x = -Lx +[0:nx-1]'*hx;

      % Integration weights and linear operators
      W = zeros(nx,nx);
      for i = 1:nx
	for j = 1:nx
	  W(i,j) = p.wFun(x(i),x(j));
	end
      end

      % Right-hand side function handle
      N = @(t,uHat) NeuralField(t,uHat,p.f,W,p.xi,x);

      % Time step
      u0Hat = fft(p.uAna(x,0)); 
      tspan = [0 3];
      [t,UHat] = ode45(N,tspan,u0Hat);

      % Comptue error
      [X,T] = meshgrid(x,t); 

      U = zeros(size(UHat));
      for i = 1:length(t)
	U(i,:) = real(ifft(UHat(i,:)));
      end
      eVec(m) = max(max(abs(U-p.uAna(X,T))));

    end
      figure(fig), hold on; plot(nVals,eVec,'.-','DisplayName',pFlag); hold off; drawnow;

  end

  set(gca,'XScale','log','YScale','log'); box on; grid on; legend;
  title('Spectral Galerkin (Trapezium)');

  savefigure

end
