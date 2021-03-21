%%%%%%%%%%%%%%%%%%%% Hansen's RBC model, VFI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  clear
  close all
  clc

%%%%%%%%%%%%%%%%%%%%%%%%% parameter values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  beta = 0.99;
  gamma = 0.0045;
  eta = 1.0039;
  theta = 0.2342;
  A = 6.0952;
  delta = 0.025;
  rho = 0.95;
  sigmasqr = 0.00025;
  sigma = sqrt(sigmasqr);
  
%%%%%%%%%%%%%%%%%%%%%%% steady state values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  k_ss = 2.2629e+04;
  c_ss = 3.1183e+03;
  y_ss = 3.7722e+03;
  inv_ss = 653.9894;
  h_ss = 205.8690;
  
%%%%%%%%%%%%%%%%%%%%%%%%%%% state spaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  nk = 100;
  k = ones(nk , 1); % state space for capital
  k(1) = k_ss * 0.99; % lower bound of the grid, ss value was computed to be 22629
  k(nk) = k_ss * 1.01; % upper bound of the grid
  for i = 2 : nk
      k_dist = (k(nk) - k(1))/(nk - 1);
      k(i) = k(1) + (i - 1) * k_dist; % evenly spaced grid for capital
  end
  
  
  N = 7;
  m = 3;
  [Trans , s] = markovapprox(rho , sigma , m , N); % Trans is the transition matrix and s is a discretized vector

%%%%%%%%%%%%%%%%%%%%%%%%%% compute for labor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  h = zeros(nk , N , nk); % The first dimension is for today's capital
                          % The second dimension is for the shocks
                          % The third dimension is for tomorrow's capital
                          % The same holds for all the other 3-D arrays
                   
  
  for i = 1 : nk
    for j = 1 : nk
        for x = 1 : N
            FOC = @(h0) ((1 - theta) * (A * exp(s(x))) * k(i)^theta * h0^(1 - theta))...
                /(gamma * ((A * exp(s(x))) * k(i)^theta * h0^(1 - theta) + (1 - delta) * k(i) - eta * k(j))) - h0;
            h(i , x , j) = fzero(FOC , h_ss);     
        end
    end
  end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%% return matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  R = ones(nk , N , nk); % 3 dimensional return matrix
  C = ones(nk , N , nk); % 3 dimensional consumption matrix

for i = 1 : nk
    for j = 1 : nk
        for x = 1 : N
            
            c = A * exp(s(x)) * (k(i)^theta) * h(i , x , j)^(1 - theta)...
            + (1 - delta) * k(i) - eta * k(j);
            C(i , x , j) = c;
            R(i , x , j) = log(c.*(c>0)).*(c>=0)-999*(c<0) - gamma * h(i , x , j);
            % setting a large negative value to the util function when c<0
           
        end
    end
end
      
%%%%%%%%%%%%%%%%%%%%%%% value function iteration %%%%%%%%%%%%%%%%%%%%%%%%%

  eps = 0.001;              % convergence criteria
  MxIt = 1000;              % maximum numbers of iteration
  
  V = ones(nk , N);         % state space for value function
  Vnew = ones(nk , N);      % value function to be updated
  
  Z = ones(nk , N , nk);    % Z temporarily stores the discounted expected future value
  
  for x = 1 : MxIt
      
      M = beta * Trans * V'; % the discounted expected future value
      
      for f = 1 : nk
          Z(f , : , :) = M; % replicate the matrix along the first dimension, as it is the dimension for today's capital
      end
      
      W = R + Z;            % the RHS of the Bellman Equation
      
      [W , index] = max(W , [] , 3); % maximize over the 3rd dimension, which is tomorrow's capital
      
      Vnew = W;
      Pol = index';
      
  if max(abs(Vnew - V)) < eps
      fprintf(1 , 'The value function has CONVERGED in %2.0f iterations.\n' , x);
      break;
  end
  V = Vnew;          % Update the value function if criteria not met
  end
  
%%%%%%%%%%%%% plot the value function and the policy function %%%%%%%%%%%%

  V_low = V(: , 1);
  V_high = V(: , 7);
  
  k_Pol_low = Pol(1 , :);   % policy function for capital in low state
  k_Pol_high = Pol(7 , :);  % policy function for capital in high state
  
  h_low = squeeze(h(: , 1 , :)); 
  h_high = squeeze(h(: , 7 , :));
  
  for i = 1 : nk
      h_Pol_low(i) = h_low(i , Pol(1 , i)); % policy function for labor in low state
  end
  
  for i = 1 : nk
      h_Pol_high(i) = h_high(i , Pol(7 , i)); % policy function for labor in high state
  end
  
  c_low = squeeze(C(: , 1 , :));
  c_high = squeeze(C(: , 7 , :));
  
  for i = 1 : nk
      c_Pol_low(i) = c_low(i , Pol(1 , i)); % policy function for consumption in low state
  end
  
  for i = 1 : nk
      c_Pol_high(i) = c_high(i , Pol(7 , i)); % policy function for consumption in high state
  end

  subplot(2 , 2 , 1)
  plot(k , V_low)
  hold on
  plot(k , V_high)
  legend('low state','high state')
  title('value function');
  xlabel('current capital');
  
  subplot(2 , 2 , 2)
  plot(k , k(k_Pol_low))
  hold on
  plot(k , k(k_Pol_high))
  legend('low state','high state')
  title('policy function for capital')
  xlabel('current capital')
  
  subplot(2 , 2 , 3)
  plot(k , h_Pol_low)
  hold on
  plot(k , h_Pol_high)
  legend('low state','high state')
  title('policy function for labor')
  xlabel('current capital')
  
  subplot(2 , 2 , 4)
  plot(k , c_Pol_low)
  hold on
  plot(k , c_Pol_high)
  legend('low state','high state')
  title('policy function for consumption')
  xlabel('current capital')
  
 