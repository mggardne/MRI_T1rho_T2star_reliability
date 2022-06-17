function [b,r2,df,sse,f,b1se,t] = regres3(x,y,p)
%REGRES3 Computes least squares linear regression for the model Y =
%        B_P*X^P+B_(P-1)*X^(P-1)+B_(P-2)*X^(P-2)+ ... +B_1*X+B_0.
%        B = REGRES3(X,Y,P) returns the polynomial coefficients in
%        vector B (1x(p+1)) given the independent X vector (nx1) and
%        dependent Y vector (nx1).
%
%        [B R2 DF SSE F] = REGRES3(X,Y,P) returns the coefficient of
%        determination (R2), the number of regression degrees of
%        freedom (DF), the sum of squared errors (SSE) and the F value
%        of the goodness of fit (F).
%
%        [B R2 DF SSE F B1SE T] = REGRES3(X,Y,1) returns the standard
%        error of the slope B(1,1) (B1SE) and t score (T) for the slope
%        B(1,1) of the linear fit.  NOTE that this is valid only for
%        a linear polynomial (P = 1).  NaNs are returned if P > 1.
%
%        NOTES:  1.  Polynomial power (P) is always a positive integer.
%
%                2.  Standard error (B1SE) and t score (T) of the slope
%                B(1,1) are only valid when polynomial power (P) equals
%                one (1) - a linear polynomial.
%
%                3.  The calculation of R2 assumes SS_reg+SS_err=SS_tot.
%                While correct for linear regression, this is not the
%                most general form for R2.   See p. 491 in "An
%                Introduction to Statistical Methods and Data Analysis"
%                by Lyman Ott, 3rd Ed., 1988 and:
%             http://en.wikipedia.org/wiki/Coefficient_of_determination.
%
%        25-Sep-09 * Mack Gardner-Morse
%

%#######################################################################
%
% Check for Inputs
%
if (nargin<3)
  error(' *** ERROR:  REGRES3 requires three input arguments.');
end
%
% Check that Both Data Inputs are Vectors
%
[n l] = size(x);
[n2 l2] = size(y);
%
m = min([n l]);
m2 = min([n2 l2]);
%
if ((m~=1)|(m2~=1))
  error(' *** ERROR:  REGRES3 only works with data vectors.')
end
%
x = x(:);
y = y(:);
%
[n l] = size(x);
[n2 l2] = size(y);
%
% Check that the Both Inputs have the Same Number of Rows
%
if (n~=n2)
  error(' *** ERROR:  REGRES3 requires X and Y have the same length.')
end
%
% Check Power
%
if (p<1)
  error(' *** ERROR:  REGRES3 polynomial power must be a postive integer.');
end
%
p = floor(p);
%
% Least Squares Linear Regression
%
b = polyfit(x,y,p);
%
% R^2 and t Score Statistics
%
if (nargout>1)
  ybar = mean(y);
  yhat = polyval(b,x);
  dy = y-ybar;
  tss = dy'*dy;
  dr = yhat-ybar;
  ssr = dr'*dr;
  if (tss>eps)
    r2 = ssr/tss;
  else
    r2 = 0;
  end
  df = n-p-1;
  se = norm(y-yhat)/sqrt(df);
  xbar = mean(x);
  sx = norm(x-xbar);
  b1se = se/sx;         % Standard error of slope
  if (se>eps)
    t = b(1,1)/b1se;
  else
    t = 0;
  end
  sse = se*se;
  if (sse>eps)
    f = (ssr/p)/sse;                   % Goodness of fit
  else
    f = 0;
  end
%
  if (p~=1)
    t = NaN;
    b1se = NaN;
  end
end
%
return