% Program by Sreevalsan S Menon(sm2hm@mst.edu)

clc                                                                         % Clear command window
clear                                                                       % clear workspace
rng('default')                                                              % For reproducibility
%% N=50 independent random numbers
u = rand(50,1);                                                             % Generating 50 random numbers from uniform distribution between 0 and 1
fprintf('The mean for 50 random numbers is: %d \n', mean(u));               % Outputs the mean
fprintf('The standard deviation for 50 random numbers is: %d \n\n', std(u));% Output standard deviation
subplot(1,3,1)                                                              % Add subplot
histogram(u,10)                                                             % Plot histogram with 10 bins
title('N=50')                                                               % Adding title

%% N=500 independent random numbers
u = rand(500,1);                                                            % Generating 50 random numbers from uniform distribution between 0 and 1
fprintf('The mean for 500 random numbers is: %d \n', mean(u));              % Outputs the mean
fprintf('The standard deviation for 500 random numbers is: %d \n\n', std(u));% Output standard deviation
subplot(1,3,2)                                                              % Add subplot
histogram(u,10)                                                             % Plot histogram with 10 bins
title('N=500')                                                              % Adding title

%% N=5000 independent random numbers
u = rand(5000,1);                                                           % Generating 50 random numbers from uniform distribution between 0 and 1
fprintf('The mean for 5000 random numbers is: %d \n', mean(u));             % Outputs the mean
fprintf('The standard deviation for 5000 random numbers is: %d \n\n', std(u));% Output standard deviation
subplot(1,3,3)                                                              % Add subplot
histogram(u,10)                                                             % Plot histogram with 10 bins
title('N=5000')                                                             % Adding title
