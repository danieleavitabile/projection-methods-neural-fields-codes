function p = lagrange(x,y,z)
%LAGRANGE    Lagrange interpolation.
%            LAGRANGE(X,Y,Z) Calculates the values of the Lagrange
%            interpolating polynomial, with interpolation points X and
%            values Y, at the points Z.
%            

% Get degree of polynomial from input.
    n = length(x) - 1;
% Check that x and y have the same length.
    if length(x) ~= length(y)
        error('x and y should be the same size')
    end
% Compute Lagrange basis.
    L = zeros(n+1,length(z));
    for k = 1:n+1
        L(k,:) = 1;
        for j = 1: length(x)
            if j ~= k
                L(k,:) = L(k,:) .* (z-x(j))/(x(k) - x(j));
            end
        end
    end
% Compute Lagrange interpolating polynomial.
    p = y*L;
end
