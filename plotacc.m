%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotacc.m plots data generated in main_error.m 
% data for second order accuracy in time is stored in errortime.txt
% errortime.txt is organized as follows:
% first column error, second column normalization, third column energy
% data for spectral accuracy in phase-space is stored in errorspace.txt
% errorspace.txt is organized as follows:
% first column error, second column normalization, third column energy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot 2.order accuracy in time
t = (0.1:-0.01:0.01)';  
data = dlmread('errortime.txt'); 
errtime = data(:,1); 
norm = data(:,2); 
energy = data(:,3); 
figure(1)
fit = polyfit(log(t),log(errtime),1); 
ausglgerad = exp(fit(1)*log(t) + fit(2));  
loglog(t,errtime,'-s', t,ausglgerad); grid on; 
ylabel('$\| W^n - W^n_{ex}\|$','Interpreter','LaTex','FontSize',15); 
xlabel('$\delta t$', 'Interpreter','LaTex','FontSize',15); 

%% Plot Spectral accuracy in space
x = (6:2:20)'; 
data_space = dlmread('errorspace.txt'); 
errspace = data_space(:,1); 
normspace = data_space(:,2); 
energyspace = data_space(:,3); 
figure(2)
loglog(x,errspace,'-s'); grid on; 
ylabel('$\| W^n - W^n_{ex}\|$','Interpreter','LaTex','FontSize',15); 
xlabel('$N,M$', 'Interpreter','LaTex','FontSize',15);
