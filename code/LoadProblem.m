function prob=LoadProblem(pFlag,CVal,LVal,gammaVal,thetaVal,kVal)

  %% Parameters
  if nargin == 1
    C     = 0.8;
    L     = 1;
    gamma = 0.5;
    theta = 0.3;
    k     = 5;
  elseif nargin > 1
    C     = CVal;
    L     = LVal;
    gamma = gammaVal;
    theta = thetaVal;
    k     = kVal;
  end

  %% Various function handles
  f     = @(u) 1./(1+exp(-k*(u-theta)));
  fInv  = @(u) theta - log(-1+u.^-1)/k;

  switch pFlag
    case 'P1'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2+y).*cos(y);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - ...
	       0.5*C*exp(-gamma*t-x.^2-1)*(-cos(1)+sin(1)+exp(2)*(cos(1)+sin(1)));
    case 'P2'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2).*y.^20;
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - 2/21*C*exp(-gamma*t-x.^2);
    case 'P3'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2)./(1+16*y.^2);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - 0.5*atan(4)*C*exp(-gamma*t-x.^2);
    case 'P4'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2)*exp(-y.^2);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - sqrt(pi)*erf(1)*C*exp(-gamma*t-x.^2);
    case 'P5'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2)*exp(-y); 
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - (exp(2)-1)*C*exp(-1-gamma*t-x.^2);
    case 'P6'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2)*abs(y).^3;
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - 0.5*C*exp(-gamma*t-x.^2);
    case 'P7p'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+cos(x).^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+cos(x).^2)./(C-exp(gamma*t+cos(x).^2));
      wFun  = @(x,y) exp(-cos(x).^2+cos(y).^2).*cos(y).^2;
      xi    = @(x,t) utAna(x,t) + uAna(x,t) -C*exp(-gamma*t-cos(x).^2)*pi;
    case 'P8p'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+cos(x).^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+cos(x).^2)./(C-exp(gamma*t+cos(x).^2));
      wFun  = @(x,y) exp(-cos(x).^2+cos(y).^2).*1./(1+16*cos(y).^2);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) -2*C*pi/sqrt(17)*exp(-gamma*t-cos(x).^2);
    case 'P9p'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+cos(x).^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+cos(x).^2)./(C-exp(gamma*t+cos(x).^2));
      wFun  = @(x,y) exp(-cos(x).^2+cos(y).^2).*abs(cos(y).^3);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) -8/3*C*exp(-gamma*t-cos(x).^2);
    case 'P10p'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+cos(x).^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+cos(x).^2)./(C-exp(gamma*t+cos(x).^2));
      wFun  = @(x,y) exp(-cos(x).^2+cos(y).^2).*(cos(y).^20);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) -46189/131072*pi*C*exp(-gamma*t-cos(x).^2);
    case 'P11'
      uAna  = @(x,t) (k*theta-log(-1+C^-1*exp(gamma*t+x.^2)))/k;
      utAna = @(x,t) k^-1 * gamma*exp(gamma*t+x.^2)./(C-exp(gamma*t+x.^2));
      wFun  = @(x,y) exp(-x.^2+y.^2).*exp(-1./y.^2);
      xi    = @(x,t) utAna(x,t) + uAna(x,t) - (2+2*exp(1)*sqrt(pi)*(-1+erf(1)))*C*exp(-1-gamma*t-x.^2);
    otherwise
      error('Unknown problem flag');
  end

 prob.C     = C     ;
 prob.L     = L     ;
 prob.gamma = gamma ;
 prob.theta = theta ;
 prob.k     = k     ;
 prob.f     = f     ;
 prob.fInv  = fInv  ;
 prob.uAna  = uAna  ;
 prob.utAna = utAna ;
 prob.wFun  = wFun  ;
 prob.xi    = xi    ;

end
