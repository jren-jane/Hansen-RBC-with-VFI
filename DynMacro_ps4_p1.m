
%%%%%%%%%%%%%%%%%%%%%%%% simulation and pruning %%%%%%%%%%%%%%%%%%%%%%%%%%

  clear
  close all
  clc
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  rho = 0.8;
  alpha = 0.5;
  sigma = 0.1;

  y0 = linspace(0 , 1);
  y1 = rho * y0 + alpha * y0.^2 + 2 * sigma;
  y2 = rho * y0 + alpha * y0.^2 - 2 * sigma;

  figure (1)
  plot (y0 , y1);
  hold on
  plot (y0 , y2);
  plot (y0 , y0);
  legend('+2sgima' , '-2sigma' , 'identity line')
  xlabel('yt-1')
  ylabel('yt')
  title('stochastic dynamics')
  hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% question 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  T = 500;
  sigma2 = sigma/10;
  sim_eps2 = sigma2 * randn(T , 1);
  y3 = zeros(T , 1);
  for t = 2 : T
      y3(t) = rho * y3(t - 1) + alpha * y3(t - 1) ^ 2 + sim_eps2(t);
  end
  
  figure (2)
  subplot(1 , 2 , 1)
  plot (y3)
  title('simulation: 2nd order, small sigma')

  y4 = zeros(T , 1);
  for t = 2 : T
      y4(t) = rho * y4(t - 1) + sim_eps2(t);
  end
  
  subplot(1 , 2 , 2)
  plot (y4)
  title('simulation: 1st order process, small sigma')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% question 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  sim_eps = sigma * randn(T , 1);
  y5 = zeros(T , 1);
  for t = 2 : T
      y5(t) = rho * y5(t - 1) + alpha * y5(t - 1) ^ 2 + sim_eps(t);
  end
  
  figure (3)
  subplot(1 , 2 , 1)
  plot (y5)
  title('simulation: 2nd order, large sigma')

  y6 = zeros(T , 1);
  for t = 2 : T
      y6(t) = rho * y6(t - 1) + sim_eps(t);
  end
  
  subplot(1 , 2 , 2)
  plot (y6)
  title('simulation: 1st order, large sigma')
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% question 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  ytilde = zeros(T , 1);
  for t = 2 : T
      ytilde(t) = rho * ytilde(t - 1) + sim_eps(t);
  end
  
  ypruned = zeros(T , 1);
  for t = 2 : T
      ypruned(t) = rho * ypruned(t - 1) + alpha * ytilde(t - 1) ^ 2 + sim_eps(t);
  end
  
  figure (4)
  plot(ytilde)
  hold on
  plot(ypruned)
  title('linear component and pruned')
  legend('linear component','pruned')
  hold off