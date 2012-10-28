% Plotting Scripts

% Log Plot of Lagrangian
figure(1);
semilogy(lagrangian);grid on;
title('Lagrangian Hunting (L*) per iteration');
xlabel('Lagrangian Adjustment iteration');
ylabel('log(L*)');

% Histogram of L0-norm error
figure(2);
hist(L0_norm_diff,max(L0_norm_diff));
h = findobj(gca,'Type','patch');
set(h,'FaceColor','b','EdgeColor','w');
title('\Delta L_0-norm for x_{est} and x');
xlabel('number of non-zero entries');ylabel('frequency of occurance');

% Plot L0-norm
figure(21);
plot(L0_norm_diff);
title('L_0-norm for x_{est} and x');
xlabel('Statistical Expirement');ylabel('Number of non-zero unmatched entries');

% Plot L0-norm error and Epsilon
figure(3);
scatter(L0_norm_diff,epsilon_track);

% Contour plot of lambda and y
figure(4);
subplot(1,2,1);
contour(lambda);xlabel('measurement (n)');ylabel('Spectra (m)');
title('\lambda for one Statistical Experiment (R)');
subplot(1,2,2);
contour(y);xlabel('measurement (n)');ylabel('Spectra (m)');
title('y for one Statistical Experiment (R)');

% Contour plot of x and X_sol
figure(5);
subplot(1,2,1);
contour(x);xlabel('measurement (n)');ylabel('Spectra (m)');
title('True x');
subplot(1,2,2);
contour(Xsol);xlabel('measurement (n)');ylabel('Spectra (m)');
title('x_{est} for one Statistical Experiment (R)');

% Stem plot true_x versus x_est
figure(6);
x_temp = (sum(Xsol.^2,2));
x_temp = x_temp.*(max(true_x)/max(x_temp)); % Normalize to true_x
x_temp(x_temp < 1e-3) = 0;
stem(true_x,'b');hold on;
stem(x_temp,'r');hold off
title('x_{true} and x_{estimate}');
xlabel('Spectra Index');

% Simple true_x versus x_est
figure(7);
wi = 4; % which index to look at
stem(true_x,'b');hold on;stem(x_est_saved(:,wi),'r');
title(['x_{true} and x_{estimate} with \Delta L_0-norm = ',num2str(L0_norm_diff(wi))]);
xlabel('Spectra Index');

