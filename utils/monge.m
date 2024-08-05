function [A] = monge(covmat1,covmat2)
% monge map calculation for centered multidimensional Gaussians
    As_sqrt_inv = covmat1^(-0.5);
    As_sqrt = covmat1^(0.5);
    A_int = As_sqrt*covmat2*As_sqrt;
    A_int = A_int^(0.5);
    A = As_sqrt_inv*A_int*As_sqrt_inv;

end