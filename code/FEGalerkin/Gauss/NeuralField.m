function K = NeuralField(t,u,wFun,sFun,xiFun,x)

  % Number of elements
  nElem = length(u)-1;

  % Quadrature rule
  nq = 2;
  rhoq(1) = 1.0; zq(1) = -1/sqrt(3);
  rhoq(2) = 1.0; zq(2) =  1/sqrt(3);

  K = InnerProductPhiG();

  function F = InnerProductPhiG()

    F = zeros(size(u));

    % Loop over elements
    for elem =  1:nElem

      n1 = elem;
      n2 = elem+1;

      x1 = x(n1);
      x2 = x(n2);

      dx = x2 - x1;

      % Evaluate quadrature
      for q = 1:nq

	% Get location of quadrature point
	z = zq(q);

	% Calculate x location of quadrature point
	xq = x1 + 0.5*(1+z)*dx;

	% Calculate g
	g = xiFun(xq,t) + ComputeG(xq);

	% Calculate phi1 and phi2
	phi1 = 0.5*(1-z);
	phi2 = 0.5*(1+z);

	% Add term to n1 component
	F(n1) = F(n1) + rhoq(q)*0.5*phi1*g*dx;

	% Add term to n2 component
	F(n2) = F(n2) + rhoq(q)*0.5*phi2*g*dx;

      end

    end

  end

  function G = ComputeG(y)

    G = 0;

    for elem = 1:nElem

      n1 = elem;
      n2 = elem+1;

      x1 = x(n1);
      x2 = x(n2);

      dx = x2 - x1;

      for q = 1:nq

	% Get location of quadrature point
	z = zq(q);

	% Calculate x location of quadrature point
	xq = x1 + 0.5*(1+z)*dx;

	phi1 = 0.5*(1-z);
	phi2 = 0.5*(1+z);

	G = G + wFun(y,xq)*sFun(u(n1)*phi1 + u(n2)*phi2)*rhoq(q)*0.5*dx;

      end

    end

  end

end
