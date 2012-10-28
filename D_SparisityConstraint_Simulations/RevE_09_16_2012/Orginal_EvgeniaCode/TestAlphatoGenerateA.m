% Test different sizes of alpha to generate A
clear;
close all;
clc;

b = 1:200;
alpha = 2:0.1:12;
amp = 100;

for t=1:length(alpha)
    a = b.^2.*exp(-b.^2/(2*alpha(t)^2))/alpha(t)^3; % Generate Maxwell-Boltzman Curve
    a = (amp/max(a)).*a;                      % Normalize to defined amplitude
    A = tril(toeplitz(a));
    condA(t) = cond(A);
    FrobeniusNormA(t) = norm(A,'fro');
    subplot(3,1,1);plot(a);title(['\alpha = ',num2str(alpha(t))]);
    subplot(3,1,2);loglog(condA);grid on; title(['log(cond(A)) = ',num2str(log(condA(t)))]);
    subplot(3,1,3);loglog(FrobeniusNormA); grid on; title(['FrobeniusNorm(A) = ',num2str(FrobeniusNormA(t))]);
    pause(0.05);
end