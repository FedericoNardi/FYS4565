function [ xs, xps ] = generateBeam(epsilon, beta, alpha, N)
    gamma = (1+alpha^2)/beta; % Twiss gamma
    sigx = sqrt(beta*epsilon); % rms beam size
    sigxp = sqrt(gamma*epsilon); % rms divergence
    rho = - alpha/sqrt(1+alpha^2); % correlation
    us = randn(1,N); % Gaussian random variable #1
    vs = randn(1,N); % Gaussian random variable #2
    xs = sigx*us; % beam positions
    xps = sigxp*us*rho + sigxp*vs*sqrt(1 - rho^2); % beam angles
end