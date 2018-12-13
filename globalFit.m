function [p_1,p_2] = globalFit(X,Y,tint) 
% function [p_1,p_2] = globalFit(X,Y,tint) 
%
% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Fit culmulative residence time distribution (CRTD) of binding events
% measured at multiple intervals. The fitting functions are:
%   1. Mono-exponential distribution: 
%       y = Aexp(-(kb*tint/ttl + koff1)t)
%   2. Bi-exponential distribution:
%       y = A(B*exp(-(kb*tint/ttl + koff1)t) + 
%                   (1-B)*exp(-(kb*tint/ttl + koff2)t))
%   3. Tri-exponential distribution:
%       y = A(B1*exp(-(kb*tint/ttl + koff1)t) + 
%                   B2*exp(-(kb*tint/ttl + koff2)t) +
%                       (1-B1-B2)*exp(-(kb*tint/ttl + koff3)t))
%   where A represents the number of molecules, kb and koff represent the 
%   photobleaching rate and off rate respectively (s-1), t is time (s), 
%   tint and ttl are integration time (s) and interval time (s) respectively.
% % %  Global parameters: B, kb, koff1 and koff2
% % %  A is not global
%
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X - matrix containing real time
%       + row: real time in each interval (ttl)
%       + column: multiples of ttl
% Y - matrix containing CRTD
%       + row: the number of counts in CRTDs for each interval
%       + column: counts as a function of time
% tint - integration time
% 
% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% p1 - fitting outcomes by globally fitting CRTDs using mono-exponential
% function
%       + rows represent fitting outcomes from each bootstrapped samples
%       + column 1: kb
%       + column 2: koff
%       + column 3: 1 (indicating single population)
%       + columns 4-7: zeros (to facilitate comparison with results from
%       fitting using multiple-exponential functions)
%       + columns 8-last: A at the corresponding interval
% p2 - fitting outcomes by gloablly fitting CRTDs using bi-exponential
% function
%       + rows represent fitting outcomes from each bootstrapped samples
%       + column 1: kb
%       + column 2: koff1
%       + column 3: B1
%       + column 4: koff2
%       + column 5: B2 (1 - B1)
%       + columns 6-7: zeros (to facilitate comparison with results from
%       fitting using tri-exponential functions)
%       + columns 8-last: A at the corresponding interval
%
% Authors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Han N. Ho & Daniel Zalami - December 2018
%
%--------------------------------------------------------------------------

%% Known Parameters
ttl = X(:,1); % vector containing all interval times
lower_B = 1/min(Y(:,1)); % the lower bound for the amplitudes
upper_koff = 1/tint;     % the upper bound for off rates

%% Initialization
% vector containing fitting results by globally fitting 
% using mono-exponential function
p_1 = zeros(1,length(ttl)+7); 
% vector containing fitting results by globally fitting 
% using bi-exponential function
p_2 = zeros(1,length(ttl)+7); 
a_para = 100*ones(size(ttl)); % Initialize the vector for A
weights = ones(size(X));      % fitting weights 
%n_para = 2: mono-exponential distribution
%n_para = 4: bi-exponential distribution
for n_para = 2:2:4 
    if n_para == 2   % fitting using mono-exponential function                             
            para = [1, 1, a_para']; % para  = [kb, koff, A];
            lb   = [0, 0, zeros(size(ttl))']; % lower bounds
            ub   = [Inf, upper_koff, Inf*ones(size(ttl))']; % upper bounds
            f1 = @(p)(   (model(n_para,p,X,tint,ttl)-Y).*weights ); % fitting function
            opts = optimset('Display','off');
            % Fit
            [p] = lsqnonlin(f1,para,lb,ub,opts);  
            p_1(1,:) = [p(1:2),1,zeros(1,4),p(3:end)];
    elseif n_para == 4 % fitting using bi-exponential function        
            para = [1, 1, 0.5, 2, a_para']; % para = [kb, koff1, B1, koff2, A];
            lb   = [0, 0, lower_B, 0, zeros(size(ttl))']; % lower bounds
            ub   = [Inf, upper_koff, 1-lower_B, upper_koff, Inf*ones(size(ttl))']; % upper bounds
            f1 = @(p)(   (model(n_para,p,X,tint,ttl)-Y).*weights ); % fitting function
            opts = optimset('Display','off');
            % Fit
            [p] = lsqnonlin(f1,para,lb,ub,opts);
            % assign the smaller off rate to be the slow off rate (koff1)
            if p(2) > p(4)
                assign = p(2);
                p(2) = p(4);
                p(4) = assign;
                p(3) = 1 - p(3);
            end
            p_2(1,:) = [p(1:4),1-p(3),zeros(1,2),p(5:end)];
    end
end

end

function f = model(n_para,para,X,tint,ttl)
% function f = model(n_para,para,X,tint,ttl)
%
% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Define fitting functions:
%   1. Mono-exponential distribution: 
%       y = Aexp(-(kb*tint/ttl + koff1)t)
%   2. Bi-exponential distribution:
%       y = A(B*exp(-(kb*tint/ttl + koff1)t) + 
%                   (1-B)*exp(-(kb*tint/ttl + koff2)t))
%   3. Tri-exponential distribution:
%       y = A(B1*exp(-(kb*tint/ttl + koff1)t) + 
%                   B2*exp(-(kb*tint/ttl + koff2)t) +
%                       (1-B1-B2)*exp(-(kb*tint/ttl + koff3)t))
% where A represents the number of molecules, kb and koff represent the 
% photobleaching rate and off rate respectively (s-1), t is time (s), 
% tint and ttl are integration time (s) and interval time (s) respectively.
% 
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% n_para - 2:fitting using mono-exponential function
%          4:fitting using bi-exponential function 
% para - vector containing fitting outcomes%           
% X - matrix containing real time
%       + row: real time in each interval (ttl)
%       + column: multiples of ttl
% tint - integration time
% ttl - interval time
%
% Authors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Han N. Ho & Daniel Zalami - December 2018
%
%--------------------------------------------------------------------------
p = tint./ttl; p = p(:);
if n_para == 2
    kb = para(1);
    koff1 = para(2);
    ampl = para(3:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (exp(-((kb.*p + koff1)*ones(1,size(X,2))).*X));   
elseif n_para == 4
    kb = para(1);    
    koff1 = para(2);
    B1 = para(3);
    koff2 = para(4);
    ampl = para(5:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (B1.*exp(-((kb.*p + koff1)*ones(1,size(X,2))).*X)+...
        (1-B1) .* exp( -(kb.*p + koff2)*ones(1,size(X,2)).*X ));
elseif n_para ==6

    kb = para(1);
    koff1 = para(2);    B1 = para(3);
    koff2 = para(4);    B2 = para(5);
    koff3 = para(6);
    ampl = para(7:end);
    f = (ampl'*ones(1,size(X,2))).*...
        (B1.*exp(-((kb.*p + koff1) * ones(1,size(X,2))).*X)...
        + B2.* exp( -(kb.*p + koff2)*ones(1,size(X,2)).*X )+...
        (1-B1-B2).* exp( -(kb.*p + koff3)*ones(1,size(X,2)).*X ));
end   
end
