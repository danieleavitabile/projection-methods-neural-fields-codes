function F = NeuralField(t,u,f,W,xi,x)

  %% Rename variables
  n = size(x,1)-1;
  F = zeros(size(u));

  % Evaluate integral using Clenshaw-Curtis quadrature
  for i = 1:n+1
    fy = (W(i,:))'.*f(u)/(2*n);
    g = real(fft(fy([1:n+1 n:-1:2])));
    a = [g(1); g(2:n)+g(2*n:-1:n+2); g(n+1)];
    r = 0*a'; r(1:2:end) = 2./(1-(0:2:n).^2);
    F(i) =r*a;
  end

  % Add other contributions
  F = F -u + xi(x,t);

end
