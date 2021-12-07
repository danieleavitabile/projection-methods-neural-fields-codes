function F = NeuralField(t,uHat,sFun,W,xiFun,x)

  hx = x(2)-x(1);
  u = real(ifft(uHat));
  K = hx*W*sFun(u);
  F = -uHat  + fft(K+xiFun(x,t));

end
