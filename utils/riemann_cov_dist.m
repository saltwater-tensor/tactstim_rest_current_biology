function [dist1] = riemann_cov_dist(C1,C2)
%% Distance between any two covariance matrices C1 and C2

          [E1] = eig(C1,C2);
          logE1 = log(E1);
          sqr = logE1 .* logE1;
          sumsqr = sum(sqr);
          sqrtsumsqr = sqrt(sumsqr);
          abssq = abs(sqrtsumsqr);
          dist1 = abssq;
  
end