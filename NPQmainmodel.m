function [qNo, alpha, tau] = NPQmainmodel(qN, t)
% t:       array of times.  One value for each pulse 
% qN:      array of measured npq values (1 - F'v/Fv)
% qNmodel: array of modeled npq
% These parameters are adjusted to minimize the sum of square errors: 
% param = [npq0, alpha, tau]  
% The human-friendly parameter vector param is 1-1 mapped to an
% optimization-friendly vector u to make the optimization problem
% well-posed. 
% u0 = initial guess for u
% 

% initial param guess 

qNo_0 = qN(1);
alpha_0 = [0.15, 0.75, 0.10];
tau_0 = [0.5, 8, 120];  % units: minutes

nExp = length(tau_0); 
assert( length(alpha_0) == nExp );

% initial guess vector
u0 = param2u( qNo_0, alpha_0, tau_0 ); 

% lower and upper bounds for u 

ub = [1.25*qNo_0, 2,   60, 360, ones(1,nExp)];
lb = [0.65*qNo_0, 0.01, 1, 35,zeros(1,nExp)];

% Nonlinear least squares minimization of the fit error
 

% Make sure the model works without crashing; 
% Get the initial sum of square errors 
SSE0 = sum( (qN - NPQmodel(u0, t)).^2 ); 


if exist('lsqnonlin.m','file')>0
    options = optimset('lsqnonlin');
    options.Display = 'none';
    options.TolX = 1e-3;
    options.TolFun = 1e-4;
    tic;
    u = lsqnonlin( @(u)(qN - NPQmodel(u, t)), u0, lb, ub, options );
    runTime = toc;
elseif exist('lsqcurvefit.m','file')>0
    options = optimset('lsqcurvefit');
    options.Display = 'i';
    options.TolX = 1e-3;
    options.TolFun = 1e-4;
    tic;
    u = lsqcurvefit(@(u)NPQmodel(u, t),u0,t,qN,lb,ub,options);
    runTime = toc;
else
    error('lsqnonlin or lsqcurvefit is required')
end

[qNo, alpha, tau] = u2param( u ); 
modelFit = NPQmodel(u, t); 
SSE = sum( (qN-modelFit).^2 ); 

% print the results
fprintf('Run time = %.2e s\n', runTime); 
fprintf('SSE = %.2e -> %.2e\n', SSE0, SSE);
fprintf('qNo = %.2e -> %.2e\n', qNo_0, qNo);
fprintf('alpha  tau\n');
for i=1:nExp
    fprintf('%.4f, %.2e\n', alpha(i), tau(i) );
end

plot(t, qN,'ko'); hold on
plot(t, modelFit,'k-','linewidth',2)
plot(t, qNo*alpha(1)*exp(-t/tau(1)), 'r-', 'linewidth',1)
plot(t, qNo*alpha(2)*exp(-t/tau(2)), 'b-', 'linewidth',1)
plot(t, qNo*alpha(3)*exp(-t/tau(3)), '-', 'linewidth',1,'color',[0, 0.5, 0])
legend('obs.', 'model','qE', 'qT', 'qI')
hold off

xlabel('Time (min)')
ylabel('qN')

end
%==========================================================================

function qN = NPQmodel( u, t) 

[qNo, alpha, tau] = u2param( u );
qN = qNo*alpha(1)*exp(-t/tau(1)) + qNo*alpha(2)*exp(-t/tau(2)) + qNo*alpha(3)*exp(-t/tau(3));

end
%==========================================================================
