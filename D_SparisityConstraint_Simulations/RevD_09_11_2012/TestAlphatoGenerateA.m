% Test different sizes of alpha to generate A
close all
clc

b = 1:200;
alpha = 1:0.1:14.9;
amp = 100;

for t=1:length(alpha)
    a = b.^2.*exp(-b.^2/(2*alpha(t)^2))/alpha(t)^3; % Generate Maxwell-Boltzman Curve
    a = (amp/max(a)).*a;                      % Normalize to defined amplitude
    plot(a);title(['alpha = ',num2str(alpha(t))]);
    pause;
end