function residenceTimeDistSimulation
% function residenceTimeDistSimulation
% This function is to be used in combination with the globalFit.m file.
%
% Background %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Dissociation kinetics of single DNA-bound protein molecules are measured
% by single-molecule live-cell imaging coupled to interval imaging
% (Gebhardt et al. 2013). In this method, the dissociation of proteins 
% can be modelled to mono- or multiple-exponential distribution as below:
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
% photobleaching rate and off rate respectively (s-1), t is time (s), tint
% and ttl are integration time (s) and interval time (s) respectively.
%
% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simulate culmulative residence time distribution (CRTD) of single-molecule
% fluorescence imaging data acquired using the interval imaging approach.
% Mono-exponential or multiple-exponential distributions with specified kb,
% koff, B1 and B2, and the number of counts are simulated using 
% the custom-written function 'simulate_res_time'. The simulated
% population is bootstrapped by randomly sampling 80% of the simulated
% population. From bootstrapped samples, CRTDs are constructed. Then, two 
% fitting procedures are performed:
%   
%   1. Individual CRTD for each interval is fitted to obtained the effective
%   rate (keff), from which keff*tau_tl plots are constructed. The shape of
%   the plot informs on the fitting function to be used in global fitting.
%   
%   2. CRTDs across all intervals are subject to global fitting, using the
%   custom-written function 'globalFit', in which non-linear least-squares
%   solver (trust-region-reflective algorithm) is used.
% 
% For each simulation, mean and standard deviation of fitting outcomes from
% bootstrapped samples are reported. Simulation using the same inputs are
% repeated for n_repeat_all times. 
%
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% kb - the photobleaching rate (s-1);
% B_all - matrix in which each row contains B1 and B2
%           B1 - the amplitude of the slowly dissociating population;
%               For mono-exponential distribution, B1 = 1;
%           B2 - the amplitude of the second dissociating population;
%               For mono-exponential distribution, B2 = 0;
%               For bi-exponetial distribution, B1 + B2 = 1;   
%               For tri-exponential distribution, the B2 is the amplitude
%               of the population with the intermediate rate. 
%               The amplitude of the fast dissociating population 
%               is (1 - B1 - B2);
% koff1 - the off rate in mono-exponential distributions;
%       - the slow off rate in bi-exponential distribution;
%       - the slow off rate in tri-exponential distribution;
% koff2 - the slow off rate in bi-exponential distributions;
%       - the intermediate off rate in tri-exponential distributions;
% koff3 - the slow off rate in tri-exponential distributions;
% tint  - integration time (or camera exposure time);
% ttl   - vector containing all interval times;
% n_count_all - array containing numbers of observations;
% n_repeat_all - the number of simulations performed using the same inputs;
%
% Outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% d(i).data - all randomly generated values at the interval ttl(i)
% frequency_write - CRTDs generated in each simulation across all
% intervals;
% M - keff*tau_tl table 
%       + column 1: tau_tl
%       + column 2: mean keff*tau_tl
%       + column 3: standard deviations of keff*tau_tl
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
% tau - matrix containing time constants obtained from globally fitting 
%       using mono-exponential distribution
%       + row: values of tau from different simulations using the same 
%           inputs 
%       + column: values of tau from different simulations with different
%           number of counts
% koff1_all_count - matrix containing off rates obtained from globally
%       fitting using mono-exponential distribution
%       + row: values of off rates from different simulations using the 
%           same 
%       + column: values of off rates from different simulations with 
%           different number of counts

% Examples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% tint = 0.1;
% ttl = [0.1; 0.2; 0.3; 0.4; 0.6; 1; 2; 3; 5; 8; 10];  
% kb = 7;
% koff1 = 0.01;
% koff2 = 0.1;
% koff3 = 1;
% n_count_all = [1000; 3000; 10000; 30000];
% n_repeat_all = 100;
%
% Authors %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Han N. Ho - Apr 2019
% 
%--------------------------------------------------------------------------
tic
clear
% specify directory where simulated CRTDs are stored
dirpath = 'C:\Users';
[filepath] = uigetdir(dirpath);
%% Specify inputs of the simulation
% kb is the photobleaching rate (unit: per second)
kb = 7;
% tint is integration time (unit: second)
tint = 0.1;
%ttl include all intervals to be used in the simulation (unit: second)
ttl = [0.1; 0.3; 0.7; 1; 3; 7; 10; 30; 70; 100];  
% B_all is a vector whose rows represent the amplitudes of koff1 and koff2
% (B1 and B2). The amplitude of koff3 is (1 - B1 - B2).
% koff1, koff2 and koff3 represent the off rates of sub-populations (unit:
% per second). 
%For simulations of single-exponential distribution: B_all = [1 0];
%B_all = [0.9 0.1; 0.75 0.25; 0.5 0.5; 0.25 0.75; 0.1 0.9];
B_all = [1/13 9/13];
koff1 = 0.01;
koff2 = 0.1;
koff3 = 1;
% n_count_all is an array whose values determine the number of counts (n) 
% in each simulation. n corresponds to the value of the first bin in the 
% cumulative residence time distribution and serve as the stopping 
% condition of the simulation. 
%n_count_all = [1000; 3000; 10000; 30000; 100000];
n_count_all = [1000; 3000; 10000; 30000; 100000; 1000000];
% n_repeat_all is the number of simulations using the identical inputs.
n_repeat_all = 100;
fit_model = [1; 2; 3];
n_bootstrap = 10; % the number of bootstrap repeats
A_global = 1; % A_global = 1: fit with amplitude as global parameter

% save simulation inputs as csv file
summary = [kb koff1 koff2 koff3 tint ttl' n_count_all' B_all(:)'];
dlmwrite(strcat(filepath, filesep,'Summary.csv'), summary);
%% Simulation
% Simulate single or multiple exponential distribution with koff1, koff2 
% and koff3 and their corresponding amplitudes B1, B2 and (1 - B1 - B2). 
for i_B = 1:size(B_all,1)
    B = B_all(i_B,:); B = B'; %B = [B1; B2];    
    % time constants obtained from global fitting using mono-exponential
    % function
    for u = 1:length(n_count_all)
    n_count_total = n_count_all(u); % the number of observation
    % summary of fitting outcomes by global fitting 
    koff1_summary = zeros(n_repeat_all,13); %mono-exponential
    koff2_summary = koff1_summary; %di-exponential
    koff3_summary = koff1_summary; %tri-exponential
    for n_repeat = 1:n_repeat_all
        disp(strcat('Repeat ',num2str(n_repeat),'_counts_',...
            num2str(n_count_total)));
        % CRTD table. 
        % The odd columns represent real time.
        % The even columns represent counts.
        frequency_write = zeros(10, length(ttl)*2);        
        bin = zeros(10,length(ttl));
        %% simulate CRTD for each interval
        for i = 1:length(ttl)              
            % calculate real time
            time = ttl(i)*(0:10)';            
            frequency_write(:,(i-1)*2+1) = time(1:end-1);            
            % define exponential distribution for each sub-population
            keff1 = (kb*tint/ttl(i) + koff1);
            mu1 = 1/keff1;
            keff2 = (kb*tint/ttl(i) + koff2);
            mu2 = 1/keff2;
            keff3 = (kb*tint/ttl(i) + koff3);
            mu3 = 1/keff3;                                 
            % determine the number of counts for each sub-population based
            % on the amplitudes B1 and B2
            n_count1 = round(B(1)*n_count_total);
            n_count2 = round(B(2)*n_count_total);             
            n_count3 = n_count_total - n_count1 - n_count2;
            % bin1, bin2 and bin3 represent sub-population CRTDs of koff1,
            % koff2 and koff3 sub-populations.
            % population1, population2 and population3 contain all values
            % in the corresponding sub-populations. These values are 
            % generated randomly by the simulate_res_time function.
            bin2 = zeros(10,1); population2 = [];
            bin3 = zeros(10,1); population3 = [];
            [bin1,population1] = simulate_res_time(mu1,time,n_count1);
            if n_count2 > 1
                [bin2,population2] = simulate_res_time(mu2,time,n_count2);            
            end
            if n_count3 > 1
                [bin3, population3] = simulate_res_time(mu3,time,n_count3);
            end
            % bin represents the sum CRTD from sub-population CRTDs
            bin(:,i) = bin1 + bin2 + bin3;
            % d(i).data contains the entire simulated population that
            % contributes to the sum CRTD
            d(i).data = [population1; population2; population3];            
            % write sum CRTD to the CRTD table
            frequency_write(:,(i-1)*2+2) = bin(:,i);
        end
        % save CRTD table across all intervals to csv file
        filename1 = strcat(filepath, filesep, num2str(n_count_total),...
            '_CRTD_',num2str(n_repeat),'_B1_', num2str(B(1)),'_koff1_',...
            num2str(koff1), '_B2_', num2str(B(2)), '_koff2_',...
            num2str(koff2),'.csv');
        dlmwrite(filename1, frequency_write);         
        %% Bootstrap the simulated population        
        keff_tautl  = zeros(n_bootstrap,length(ttl)); % keff_tautl from all bootstrapped samples
        time_all    = zeros(length(ttl),9); % real time for global fitting
        p1          = zeros(n_bootstrap, length(ttl) + 7);
        p2          = p1;
        p3          = p1;
        p1_A_global = zeros(n_bootstrap, 8);
        p2_A_global = p1_A_global;
        p3_A_global = p1_A_global; 
        for i_bt = 1:n_bootstrap
            counts = zeros(10,length(ttl));
            for i = 1:length(ttl)
                time = ttl(i)*(0:10)';
                time_all(i,:) = time(2:end-1)';
                clear a
                a = d(i).data;          
                % randomly select 80% of the simulated population
                if numel(a)>2                             
                    bt_length = ceil(0.8*length(a)); % length of bootstrapped samples         
                    bootstrap_sample = datasample(a,bt_length);                    
                    [counts(:,i),~] = histcounts(bootstrap_sample,time); % CRTD of the bootstrapped sample
                end
            %Fit 1 single exponential to obtain keff
            f = fit(time(2:end-1),counts(2:end,i),'exp1');
            keff_tautl(i_bt,i) = ttl(i)*f.b*(-1);
            end
%%          Global Fitting     
            for kk = 1:length(fit_model)
                if fit_model(kk) == 1      %fit mono-exponential model 
                    [p1(i_bt,:)] = globalFit(1, time_all, counts(2:end,:)', tint);            
                elseif fit_model(kk) == 2  %fit di-exponential model
                    [p2(i_bt,:)] = globalFit(2, time_all, counts(2:end,:)', tint);            
                elseif fit_model(kk) == 3  %fit tri-exponential model
                    [p3(i_bt,:)] = globalFit(3, time_all, counts(2:end,:)', tint);            
                end
            end
            if A_global == 1
                for kk = 1:length(fit_model)
                    if fit_model(kk) == 1      %fit mono-exponential model 
                        [p1_A_global(i_bt,:)] = globalFit2(1, time_all, counts(2:end,:)', tint);            
                    elseif fit_model(kk) == 2  %fit di-exponential model
                        [p2_A_global(i_bt,:)] = globalFit2(2, time_all, counts(2:end,:)', tint);            
                    elseif fit_model(kk) == 3  %fit tri-exponential model
                        [p3_A_global(i_bt,:)] = globalFit2(3, time_all, counts(2:end,:)', tint);            
                    end
                end
            end
         end   
%%      keff*tautl plot
        % calculate keff*tautl from all bootrstrapped samples
        M = zeros(length(ttl),3);        
        for i = 1:length(ttl)    
            keff_tautl_temp = keff_tautl(:,i);
            keff_tautl_temp = keff_tautl_temp(keff_tautl_temp>0);
            M(i,1) = ttl(i);                % tau_tl
            M(i,2) = mean(keff_tautl_temp); % mean keff*tautl from all bootstrapped samples
            M(i,3) = std(keff_tautl_temp);  % standard deviation of keff*tautl from all bootstrapped samples
        end
        % save tau_tl and keff*tautl as csv file
        filename2 = strcat(num2str(n_count_total),'_keff_tautl_',...
            num2str(n_repeat),'_B1_',num2str(B(1)),'_koff1_',...
            num2str(koff1), '_B2_',num2str(B(2)), '_koff2_',num2str(koff2));
        dlmwrite(strcat(filepath, filesep, filename2,'.csv'), M); 
        % plot keff*tautl plot
        figure(1);        
        shadedErrorBar(M(:,1),M(:,2),M(:,3),'-g',1);
        saveas(gcf, fullfile(filepath,filesep,strcat(filename2,'.fig')));
    %% Summary of Global Fitting
        try
        % save fitting outcomes from individual simulation to csv files
        % mono-exponential    
            filename3 = strcat(filepath, filesep, num2str(n_count_total),...
            '_1koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            
            dlmwrite(filename3, p1);
            % fitting outcomes using mono-exponential function
            kb_1        = mean(p1(:,1)); % mean photobleaching rate
            kb_sd_1     = std(p1(:,1)); % std of photobleaching rate
            koff_1      = mean(p1(:,2)); % mean off rate
            koff_sd_1   = std(p1(:,2)); % std of off rate
            koff1_summary(n_repeat,:) = [kb_1 kb_sd_1 koff_1 koff_sd_1 1 ...
                zeros(1,8)];
            
            % bi-exponential
            filename4 = strcat(filepath, filesep, num2str(n_count_total),...
            '_2koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            
            dlmwrite(filename4, p2);
            % fitting outcomes using bi-exponential function
            kb_2        = mean(p2(:,1)); % mean photobleaching rate
            kb_sd_2     = std(p2(:,1)); % std of photobleaching rate
            koff1_2     = mean(p2(:,2)); % mean of the slow off rate
            koff1_sd_2  = std(p2(:,2)); % std of the slow off rate
            B1_2        = mean(p2(:,3)); % mean amplitude of the slow off rate
            B1_sd_2     = std(p2(:,3)); % std of the amplitude of the slow off rate
            koff2_2     = mean(p2(:,4)); % mean of the fast off rate
            koff2_sd_2  = std(p2(:,4)); % std of the fast off rate
            B2_2        = 1 - B1_2; % the amplitude of the fast off rate
            koff2_summary(n_repeat,:) = [kb_2 kb_sd_2 koff1_2 koff1_sd_2 ...
                B1_2 B1_sd_2 koff2_2 koff2_sd_2 B2_2 zeros(1,4)]; 
            
            % tri-exponential
            filename5 = strcat(filepath, filesep, num2str(n_count_total),...
            '_3koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            
            dlmwrite(filename5, p3);
            
            kb_3        = mean(p3(:,1)); % mean photobleaching rate
            kb_sd_3     = std(p3(:,1)); % std of photobleaching rate
            koff1_3     = mean(p3(:,2)); % mean of the slow off rate
            koff1_sd_3  = std(p3(:,2)); % std of the slow off rate
            B1_3        = mean(p3(:,3)); % mean amplitude of the slow off rate
            B1_sd_3     = std(p3(:,3)); % std of the amplitude of the slow off rate
            koff2_3     = mean(p3(:,4)); % mean of the intermediate off rate
            koff2_sd_3  = std(p3(:,4)); % std of the fast intermediate off rate
            B2_3        = mean(p3(:,5)); % mean amplitude of the intermediate off rate
            B2_sd_3     = std(p3(:,5)); % mean amplitude of the intermediate off rate
            koff3_3     = mean(p3(:,6)); % mean of the fast off rate
            koff3_sd_3  = std(p3(:,6)); % std of the fast off rate
            B3_3        = 1 - B1_3 - B2_3; % amplitude of the fast off rate
            koff3_summary(n_repeat,:) = [kb_3 kb_sd_3 koff1_3 koff1_sd_3 ...
                B1_3 B1_sd_3 koff2_3 koff2_sd_3 B2_3 B2_sd_3 koff3_3 ...
                koff3_sd_3 B3_3];
        catch
        end
        
    %% Summary of Global Fitting with amplitude as global parameter
    if A_global == 1    
        try
        % save fitting outcomes from individual simulation to csv files
        % mono-exponential    
            filename3 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_1koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            
            dlmwrite(filename3, p1_A_global); 
            % fitting outcomes using mono-exponential function
            kb_1        = mean(p1_A_global(:,1)); % mean photobleaching rate
            kb_sd_1     = std(p1_A_global(:,1)); % std of photobleaching rate
            koff_1      = mean(p1_A_global(:,2)); % mean off rate
            koff_sd_1   = std(p1_A_global(:,2)); % std of off rate
            A_koff1_summary(n_repeat,:) = [kb_1 kb_sd_1 koff_1 koff_sd_1 1 ...
                zeros(1,8)];
            
            filename4 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_2koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            % bi-exponential
            dlmwrite(filename4, p2_A_global);
            % fitting outcomes using bi-exponential function
            kb_2        = mean(p2_A_global(:,1)); % mean photobleaching rate
            kb_sd_2     = std(p2_A_global(:,1)); % std of photobleaching rate
            koff1_2     = mean(p2_A_global(:,2)); % mean of the slow off rate
            koff1_sd_2  = std(p2_A_global(:,2)); % std of the slow off rate
            B1_2        = mean(p2_A_global(:,3)); % mean amplitude of the slow off rate
            B1_sd_2     = std(p2_A_global(:,3)); % std of the amplitude of the slow off rate
            koff2_2     = mean(p2_A_global(:,4)); % mean of the fast off rate
            koff2_sd_2  = std(p2_A_global(:,4)); % std of the fast off rate
            B2_2        = 1 - B1_2; % the amplitude of the fast off rate
            A_koff2_summary(n_repeat,:) = [kb_2 kb_sd_2 koff1_2 koff1_sd_2 ...
                B1_2 B1_sd_2 koff2_2 koff2_sd_2 B2_2 zeros(1,4)]; 
            
            filename5 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_3koff_',num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
            % tri-exponential
            dlmwrite(filename5, p3_A_global);
            
            % fitting outcomes using tri-exponential function 
            kb_3        = mean(p3_A_global(:,1)); % mean photobleaching rate
            kb_sd_3     = std(p3_A_global(:,1)); % std of photobleaching rate
            koff1_3     = mean(p3_A_global(:,2)); % mean of the slow off rate
            koff1_sd_3  = std(p3_A_global(:,2)); % std of the slow off rate
            B1_3        = mean(p3_A_global(:,3)); % mean amplitude of the slow off rate
            B1_sd_3     = std(p3_A_global(:,3)); % std of the amplitude of the slow off rate
            koff2_3     = mean(p3_A_global(:,4)); % mean of the intermediate off rate
            koff2_sd_3  = std(p3_A_global(:,4)); % std of the fast intermediate off rate
            B2_3        = mean(p3_A_global(:,5)); % mean amplitude of the intermediate off rate
            B2_sd_3     = std(p3_A_global(:,5)); % mean amplitude of the intermediate off rate
            koff3_3     = mean(p3_A_global(:,6)); % mean of the fast off rate
            koff3_sd_3  = std(p3_A_global(:,6)); % std of the fast off rate
            B3_3        = 1 - B1_3 - B2_3; % amplitude of the fast off rate
            A_koff3_summary(n_repeat,:) = [kb_3 kb_sd_3 koff1_3 koff1_sd_3 ...
                B1_3 B1_sd_3 koff2_3 koff2_sd_3 B2_3 B2_sd_3 koff3_3 ...
                koff3_sd_3 B3_3];
        catch
        end
    end    
    end
    % save fitting outcomes from all simulations to csv file
    try
        % summary of tau and amplitude across counts
        % from mono-exponential fit
        filename5 = strcat(filepath, filesep, num2str(n_count_total),...
            '_summary_1koff_', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename5, koff1_summary);
        tau_all(:,u) = 1./koff1_summary(:,3);
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_',num2str(B(2)),...
        '_Summary_1koff_tau.csv'), tau_all);
    
        filename6 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_summary_1koff_', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename6, A_koff1_summary);
        A_tau_all(:,u) = 1./A_koff1_summary(:,3);
        dlmwrite(strcat(filepath, filesep,'AG', num2str(B(1)),'_',num2str(B(2)),...
        '_Summary_1koff_tau.csv'), A_tau_all);
         
    catch
    end
    % from bi-exponential fit
    try
        filename6 = strcat(filepath, filesep, num2str(n_count_total),...
            '_summary_2koff', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename6, koff2_summary);
        
        tau1_2_all(:,u) = 1./koff2_summary(:,3);
        B1_2_all(:,u)   =    koff2_summary(:,5);
        tau2_2_all(:,u) = 1./koff2_summary(:,7);
        
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_tau1.csv'), tau1_2_all);
    % save time constants (from globaly fitting using mono-exponential function) as csv file
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_B1.csv'), B1_2_all);
    % save time constants (from globaly fitting using mono-exponential function) as csv file
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_tau2.csv'), tau2_2_all);
        
        filename7 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_summary_2koff', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename7, A_koff2_summary);
        
        A_tau1_2_all(:,u) = 1./A_koff2_summary(:,3);
        A_B1_2_all(:,u)   =    A_koff2_summary(:,5);
        A_tau2_2_all(:,u) = 1./A_koff2_summary(:,7);
        dlmwrite(strcat(filepath, filesep,'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_tau1.csv'), A_tau1_2_all);
    % save time constants (from globaly fitting using mono-exponential function) as csv file
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_B1.csv'), A_B1_2_all);
    % save time constants (from globaly fitting using mono-exponential function) as csv file
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_2koff_tau2.csv'), A_tau2_2_all);   
    catch
    end
    
    try
        filename7 = strcat(filepath, filesep, num2str(n_count_total),...
            '_summary_3koff', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename7, koff3_summary);
        tau1_3_all(:,u) = 1./koff3_summary(:,3);
        B1_3_all(:,u)   =    koff3_summary(:,5);
        tau2_3_all(:,u) = 1./koff3_summary(:,7);
        B2_3_all(:,u)   =    koff3_summary(:,9);
        tau3_3_all(:,u) = 1./koff3_summary(:,11);
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_B1.csv'), B1_3_all);
    
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_B2.csv'), B2_3_all);
    
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau1.csv'), tau1_3_all);
    
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau2.csv'), tau2_3_all);
    
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau3.csv'), tau3_3_all);
        
        filename8 = strcat(filepath, filesep, 'A', num2str(n_count_total),...
            '_summary_3koff', num2str(n_repeat),'_B1_',...
            num2str(B(1)),'_koff1_',num2str(koff1), '_B2_',...
            num2str(B(2)), '_koff2_',num2str(koff2),'.csv');
        dlmwrite(filename8, A_koff3_summary);
        A_tau1_3_all(:,u) = 1./A_koff3_summary(:,3);
        A_B1_3_all(:,u)   =    A_koff3_summary(:,5);
        A_tau2_3_all(:,u) = 1./A_koff3_summary(:,7);
        A_B2_3_all(:,u)   =    A_koff3_summary(:,9);
        A_tau3_3_all(:,u) = 1./A_koff3_summary(:,11);
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_B1.csv'), A_B1_3_all);
    
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_B2.csv'), A_B2_3_all);
    
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau1.csv'), A_tau1_3_all);
    
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau2.csv'), A_tau2_3_all);
    
        dlmwrite(strcat(filepath, filesep, 'AG', num2str(B(1)),'_', num2str(B(2)),...
        '_Summary_3koff_tau3.csv'), A_tau3_3_all);
    
    catch
    end 
    end
%% Bee swarm plot for distribution of tau
    try    
        Colors = zeros(length(n_count_all),3);
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_1koff_Tau_distribution.fig')));
    
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau1_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_2koff_Tau1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau2_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_2koff_Tau2_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(B1_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_2koff_B1_distribution.fig')));
 
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau1_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_3koff_Tau1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau2_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_3koff_Tau2_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(tau3_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_3koff_Tau3_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(B1_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_3koff_B1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(B2_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL',num2str(B(1)), '_',...
            num2str(B(2)),'_3koff_B2_distribution.fig')));
    catch
    end
    close all;
           %% standard deviations of tau values and amplitudes
    try
 
        tau_temp    = tau_all - 1/koff1;
        tau_temp    = tau_temp.*tau_temp;
        sd_tau      = sqrt(sum(tau_temp,1)./n_repeat_all);
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_',...
            num2str(B(2)),'_1koff_Summary_tau_sd.csv'), sd_tau);
   
        tau1_2_temp = tau1_2_all - 1/koff1; 
        tau1_2_temp = tau1_2_temp.*tau1_2_temp;
        tau2_2_temp = tau2_2_all - 1/koff2; 
        tau2_2_temp = tau2_2_temp.*tau2_2_temp;
        sd_tau1_2   = sqrt(sum(tau1_2_temp,1)./n_repeat_all);
        sd_tau2_2   = sqrt(sum(tau2_2_temp,1)./n_repeat_all);
        B1_2_temp   = B1_2_all - B(1); 
        B1_2_temp   = B1_2_temp.*B1_2_temp;
        sd_B1_2     = sqrt(sum(B1_2_temp,1)./n_repeat_all);
        sd_2_table  = [n_count_all sd_B1_2' sd_tau1_2' sd_tau2_2'];
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_',...
            num2str(B(2)),'_2koff_Summary_tau1_tau2_amp_sd.csv'), sd_2_table);
    
        tau1_3_temp = tau1_3_all - 1/koff1; 
        tau1_3_temp = tau1_3_temp.*tau1_3_temp;
        tau2_3_temp = tau2_3_all - 1/koff2; 
        tau2_3_temp = tau2_3_temp.*tau2_3_temp;
        tau3_3_temp = tau3_3_all - 1/koff3; 
        tau3_3_temp = tau3_3_temp.*tau3_3_temp;
        sd_tau1_3   = sqrt(sum(tau1_3_temp,1)./n_repeat_all);
        sd_tau2_3   = sqrt(sum(tau2_3_temp,1)./n_repeat_all);
        sd_tau3_3   = sqrt(sum(tau3_3_temp,1)./n_repeat_all);
        B1_3_temp   = B1_3_all - B(1); 
        B1_3_temp   = B1_3_temp.*B1_3_temp;
        B2_3_temp   = B2_3_all - B(2); 
        B2_3_temp   = B2_3_temp.*B2_3_temp;
        sd_B1_3     = sqrt(sum(B1_3_temp,1)./n_repeat_all);
        sd_B2_3     = sqrt(sum(B2_3_temp,1)./n_repeat_all);
        sd_3_table = [n_count_all  sd_B1_3' sd_B2_3' sd_tau1_3'...
            sd_tau2_3' sd_tau3_3'];
        dlmwrite(strcat(filepath, filesep,'AL',num2str(B(1)),'_',...
            num2str(B(2)),'_3koff_Summary_tau_amp_sd.csv'), sd_3_table);
    catch
        
    end
    %% Bee swarm plot - Amplitude as a global parameter
    try
        Colors = zeros(length(n_count_all),3);
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep, strcat('AG', num2str(B(1)),...
            '_', num2str(B(2)),'_1koff_Tau_distribution.fig')));
    
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau1_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep, strcat('AG', num2str(B(1)),...
        '_', num2str(B(2)),'_2koff_Tau1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau2_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep, strcat('AG', num2str(B(1)),...
            '_', num2str(B(2)),'_2koff_Tau2_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_B1_2_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG', num2str(B(1)),...
            '_',num2str(B(2)),'_2koff_B1_distribution.fig')));
   
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau1_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG',num2str(B(1)),...
            '_',num2str(B(2)),'_3koff_Tau1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau2_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG',num2str(B(1)),...
            '_',num2str(B(2)),'_3koff_Tau2_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_tau3_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG',num2str(B(1)), ...
            '_', num2str(B(2)),'_3koff_Tau3_distribution.fig')));
        
        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_B1_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG',num2str(B(1)),...
            '_',num2str(B(2)),'_3koff_B1_distribution.fig')));   

        figure;
        label = sprintfc('%.f',n_count_all);
        UnivarScatter(A_B2_3_all,'Label',label,'MarkerFaceColor',Colors);
        saveas(gcf, fullfile(filepath,filesep,strcat('AG',num2str(B(1)),...
            '_',num2str(B(2)),'_3koff_B2_distribution.fig')));
    catch
    end
        close all;
%% standard deviations of tau values and amplitudes
    try     
        A_tau_temp    = A_tau_all - 1/koff1;
        A_tau_temp    = A_tau_temp.*A_tau_temp;
        A_sd_tau      = sqrt(sum(A_tau_temp,1)./n_repeat_all);
        dlmwrite(strcat(filepath, filesep,'AG',num2str(B(1)),'_',...
            num2str(B(2)),'_1koff_Summary_tau_sd.csv'), A_sd_tau);
        sd_ratio_1 = sd_tau./A_sd_tau;
        figure;
        bar(sd_ratio_1);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL_over_AG',num2str(B(1)),...
            '_',num2str(B(2)),'_1koff_tau_distribution.fig')));
        
        A_tau1_2_temp = A_tau1_2_all - 1/koff1; 
        A_tau1_2_temp = A_tau1_2_temp.*A_tau1_2_temp;
        A_tau2_2_temp = A_tau2_2_all - 1/koff2; 
        A_tau2_2_temp = A_tau2_2_temp.*A_tau2_2_temp;
        A_sd_tau1_2   = sqrt(sum(A_tau1_2_temp,1)./n_repeat_all);
        A_sd_tau2_2   = sqrt(sum(A_tau2_2_temp,1)./n_repeat_all);
        A_B1_2_temp   = A_B1_2_all - B(1); 
        A_B1_2_temp   = A_B1_2_temp.*A_B1_2_temp;
        A_sd_B1_2     = sqrt(sum(A_B1_2_temp,1)./n_repeat_all);
        A_sd_2_table  = [n_count_all A_sd_B1_2' A_sd_tau1_2' A_sd_tau2_2'];
        dlmwrite(strcat(filepath, filesep,'AG',num2str(B(1)),'_',...
            num2str(B(2)),'_2koff_Summary_tau1_tau2_amp_sd.csv'), A_sd_2_table);
        sd_ratio_2 = sd_2_table(:,2:end)./A_sd_2_table(:,2:end);
        figure;
        bar(sd_ratio_2);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL_over_AG',num2str(B(1)),...
            '_',num2str(B(2)),'_2koff.fig')));
        
        
        A_tau1_3_temp = A_tau1_3_all - 1/koff1; 
        A_tau1_3_temp = A_tau1_3_temp.*A_tau1_3_temp;
        A_tau2_3_temp = A_tau2_3_all - 1/koff2; 
        A_tau2_3_temp = A_tau2_3_temp.*A_tau2_3_temp;
        A_tau3_3_temp = A_tau3_3_all - 1/koff3; 
        A_tau3_3_temp = A_tau3_3_temp.*A_tau3_3_temp;
        A_sd_tau1_3   = sqrt(sum(A_tau1_3_temp,1)./n_repeat_all);
        A_sd_tau2_3   = sqrt(sum(A_tau2_3_temp,1)./n_repeat_all);
        A_sd_tau3_3   = sqrt(sum(A_tau3_3_temp,1)./n_repeat_all);
        A_B1_3_temp   = A_B1_3_all - B(1); 
        A_B1_3_temp   = A_B1_3_temp.*A_B1_3_temp;
        A_B2_3_temp   = A_B2_3_all - B(2); 
        A_B2_3_temp   = A_B2_3_temp.*A_B2_3_temp;
        A_sd_B1_3     = sqrt(sum(A_B1_3_temp,1)./n_repeat_all);
        A_sd_B2_3     = sqrt(sum(A_B2_3_temp,1)./n_repeat_all);
        A_sd_3_table = [n_count_all A_sd_B1_3' A_sd_B2_3' A_sd_tau1_3'...
            A_sd_tau2_3' A_sd_tau3_3'];
        dlmwrite(strcat(filepath, filesep,'AG',num2str(B(1)),'_',...
            num2str(B(2)),'_3koff_Summary_tau_amp_sd.csv'), A_sd_3_table);   
    catch 
    end
    
    try
        sd_ratio_3 = sd_3_table(:,2:end)./A_sd_3_table(:,2:end);
        figure;
        bar(sd_ratio_3);
        saveas(gcf, fullfile(filepath,filesep,strcat('AL_over_AG',num2str(B(1)),...
            '_',num2str(B(2)),'_3koff.fig')));
        close all; 
    catch
    end
end

end
 
function [counts, each_molecule] = simulate_res_time(mu,edges,n_count)
% function [counts, each_molecule] = simulate_res_time(mu,edges,n_count)
%
% Purpose %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate random numbers from the exponential distribution with mean
% parameter mu, using the built-in exprnd function. The simulation was
% carried out in several rounds, and the histogram was constructed 
% following each round with edges being multiple of tau_tl. The simulation
% stops when the number of numbers in the first bin (0 to 1*tau_tl) is
% larger than the specified n_count.
%
% Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mu - 1/keff or 1/(kb*tint/ttl + koff)
% edges - vector of 11 values, (0:10).*ttl;
% n_count - the number of counts in the first bin
%
% Han N. Ho & Joris M. H. Goudsmits - December 2018
%
%--------------------------------------------------------------------------
each_molecule = [];
counts = zeros(10,1);
nx = 0;
while counts(1) < n_count
    sim = exprnd(mu,round(n_count/2.71),1);
    each_molecule = [each_molecule; sim];
    [N,~] = histcounts(sim,edges);
    counts = counts + N';
    nx = nx + 1;
end
nx
end


function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
% function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
%
% Purpose
% Makes a 2-d line plot with a pretty shaded error bar made
% using patch. Error bar color is chosen automatically.
%
% Inputs
% x - vector of x values [optional, can be left empty]
% y - vector of y values or a matrix of n observations by m cases
%     where m has length(x);
% errBar - if a vector we draw symmetric errorbars. If it has a size
%          of [2,length(x)] then we draw asymmetric error bars with
%          row 1 being the upper bar and row 2 being the lower bar
%          (with respect to y). ** alternatively ** errBar can be a
%          cellArray of two function handles. The first defines which
%          statistic the line should be and the second defines the
%          error bar.
% lineProps - [optional,'-k' by default] defines the properties of
%             the data line. e.g.:
%             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
% transparent - [optional, 0 by default] if ==1 the shaded error
%               bar is made transparent, which forces the renderer
%               to be openGl. However, if this is saved as .eps the
%               resulting file will contain a raster not a vector
%               image.
%
% Outputs
% H - a structure of handles to the generated plot objects.
%
%
% Examples
% y=randn(30,80); x=1:size(y,2);
% shadedErrorBar(x,mean(y,1),std(y),'g');
% shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});
% shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});
%
% Overlay two transparent lines
% y=randn(30,80)*10; x=(1:size(y,2))-40;
% shadedErrorBar(x,y,{@mean,@std},'-r',1);
% hold on
% y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
% shadedErrorBar(x,y,{@mean,@std},'-b',1);
% hold off
%
%
% Rob Campbell - November 2009



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error checking
error(nargchk(3,5,nargin))


%Process y using function handles if needed to make the error bar
%dynamically
if iscell(errBar)
    fun1=errBar{1};
    fun2=errBar{2};
    errBar=fun2(y);
    y=fun1(y);
else
    y=y(:)';
end

if isempty(x)
    x=1:length(y);
else
    x=x(:)';
end


%Make upper and lower error bars if only one was specified
if length(errBar)==length(errBar(:))
    errBar=repmat(errBar(:)',2,1);
else
    s=size(errBar);
    f=find(s==2);
    if isempty(f), error('errBar has the wrong size'), end
    if f==2, errBar=errBar'; end
end

if length(x) ~= length(errBar)
    error('length(x) must equal length(errBar)')
end

%Set default options
defaultProps={'-k'};
if nargin<4, lineProps=defaultProps; end
if isempty(lineProps), lineProps=defaultProps; end
if ~iscell(lineProps), lineProps={lineProps}; end

if nargin<5, transparent=0; end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot to get the parameters of the line
H.mainLine=plot(x,y,lineProps{:});


% Work out the color of the shaded region and associated lines
% Using alpha requires the render to be openGL and so you can't
% save a vector image. On the other hand, you need alpha if you're
% overlaying lines. There we have the option of choosing alpha or a
% de-saturated solid colour for the patch surface .

col=get(H.mainLine,'color');
edgeColor=col+(1-col)*0.55;
patchSaturation=0.15; %How de-saturated or transparent to make patch
if transparent
    faceAlpha=patchSaturation;
    patchColor=col;
    set(gcf,'renderer','openGL')
else
    faceAlpha=1;
    patchColor=col+(1-col)*(1-patchSaturation);
    set(gcf,'renderer','painters')
end


%Calculate the error bars
uE=y+errBar(1,:);
lE=y-errBar(2,:);


%Add the patch error bar
holdStatus=ishold;
if ~holdStatus, hold on,  end


%Make the patch
yP=[lE,fliplr(uE)];
xP=[x,fliplr(x)];

%remove nans otherwise patch won't work
xP(isnan(yP))=[];
yP(isnan(yP))=[];


H.patch=patch(xP,yP,1,'facecolor',patchColor,...
    'edgecolor','none',...
    'facealpha',faceAlpha);


%Make pretty edges around the patch.
H.edge(1)=plot(x,lE,'-','color',edgeColor);
H.edge(2)=plot(x,uE,'-','color',edgeColor);

%Now replace the line (this avoids having to bugger about with z coordinates)
delete(H.mainLine)
H.mainLine=plot(x,y,lineProps{:});


if ~holdStatus, hold off, end


if nargout==1
    varargout{1}=H;
end
end

function [xPositions, yPositions, Label, RangeCut] = UnivarScatter(data, varargin)
%--------------------------------------------------------------------------
% UnivariateScatter - Draw an Univariate Scatter plot out of a nx2 table with a
%                     categorical/string and a numerical variable, or out of a
%                     numerical array, with groups/categories corresponding
%                     to the columns of the numerical array.
%                     Many custom options are available as Name,Value
%                     pairs. For optimal visualization of your data I
%                     recommend to play with 'RangeCut', and also with pbaspect of the plot, it
%                     really changes the appearance.
%                     You need the following functions for this function to
%                     work:
%                                              
%                       +CatTable2Array: included in the file
%
%                     Also, a very simple function to assign colors is
%                     provided, ColorCoder, you can see an example in the
%                     script attached
%
%  Input: data       - Can take two types of input:
%
%                       + a nx2 table with a categorical/string and a numerical
%                       variable, the groups will be made according to the
%                       categories
%
%                       + a numerical array in which the different columns
%                       correspond to the groups/categories of the Univariate 
%                       Scatter plot in x axis. When plotting data with
%                       different number of points for the different
%                       groups, the empty elements of the array should be
%                       filled with nan. In order to do this much faster I
%                       recommend the function padcat(available in
%                       mathworks, Copyright (c) 2009, Jos van der Geest), which already
%                       concatenates row/column vectors of different
%                       lengths and fills with nans the missing values.
%
%         Name,Value  - Name-Value pairs arguments to customize the plot:
%
%                       +'Label' - Only when data is a numerical array
%                       A cell array of strings containing the
%                       labels that you want to give to categories/groups
%                       in the plot, by default they are 1,2,3...n The
%                       number of columns of this array must be the same as
%                       the number of columns in data, example:
%                       UniVarScatter(data,'Label',{'a','b','c'})
%
%                       +'Whiskers' - A string which can take 3 values
%                           *'lines': std, 95% SEM, mean are represented as
%                           lines
%                           *'boxes':std, 95% SEM, mean are represented as
%                           boxes. Default value
%                           *'none':none of these are shown
%
%                       +'PointStyle' - A String containing the different
%                       values that can be given as markertype for scatter,
%                       such as 'o','*''.', etc. Default is 'o'
%
%                       +'Width' - A number that determines the spread of
%                       each subset of points in x axis. Default is 0.4
%
%                       +'Compression' - A number that determines how
%                       separated from each other the points are in x
%                       dimension. Default is 5
%
%                       +'RangeCut' - A number or column/row vector with
%                       length= size(data,2) that determines how broad
%                       groups of points are made in y dimension, groups
%                       width in y dimension corresponds to 1/RangeCut of the
%                       95% confidence interval of a normal distribution
%                       with std and mean of each subset of points.
%                       Default is 2+2*sqrt(number of points) for each subset
%                       of points, this has shown to be a good value, but
%                       it does not perform optimally for every dataset, so
%                       I strongly recommend not to discard the graph and
%                       try different combinations of values if you don't
%                       obtained a nice result.
%
%                       +'WhiskersWidthRatio' - Lines/boxes that represent
%                       the mean,SEM,std will be Width/WhiskersWidthRatio
%                       long in x dimension. Default is 4
%
%                       +'PointSize' - A number indicating the point size.
%                       Default is 36
%
%                       +'LineWidth' - A number indicating the width of the
%                       edge of the points. Default is 0.5
%
%                       +'MarkerEdgeColor' color(s) of the edges of the points; def. 'k'
%                       +'MarkerFaceColor' color(s) of the points; def. 'flat'
%                       +'MeanColor' color(s) of the mean line; def. black
%                       +'SEMColor'  color(s) of the SEM box; def. dark gray
%                       +'StdColor'  color(s) of the Std box; def. light gray
%
%                       All this last Color Properties behave similarly,
%                       they can take several types of values
%
%                           *RGB vector 1x3: all the representations
%                           are made in the same colour represented by that vector
%                     
%                           *A string: it may contain the name of the color
%                           for instance 'k' or 'black', also it can take
%                           the values 'none', or 'flat'
%
%                           *RGB array nx3 where n>=size(data,2): each row
%                           of this array corresponds to the color of the n
%                           group/category in x. You can use the function
%                           ColorCoder included in the zip to Easily create
%                           a custom RGB array by chosing the colors from
%                           the matlab palette
%
%  Output: the figure
%
%          xPositions - An Array with the x  positions of the points, with
%                       the same size as yPositions
%                           
%          
%          yPositions - a numerical array in which the different columns
%                       correspond to the groups/categories of the Univariate 
%                       Scatter plot in x axis, its the same as data if data was
%                       already a numerical array, and if it was a table,
%                       it generates a numerical array, and the columns
%                       correspond to the Label groups.
%       
%          Label      - A string cell array with the names of the groups
%                       specially useful when using a table as input, so
%                       you know the order in which the scatter plots are
%                       showed.
%
%          RangeCut   - The value of RangeCut(explained at the Name,Values
%                       section), specially useful to optimize the
%                       representation of your data if the RangeCut value
%                       generated by the function does not perform well for
%                       your data.
%
% version: <1.00> from 30/11/2015
% 
%       Manuel Lera Ramrez: manulera14@gmail.com
%
%       Developed during Master Thesis Internship in Marcos Gonzlez-Gaitn's Lab
%       Departments of Biochemistry and Molecular Biology, Sciences II, 
%       30 Quai Ernest-Ansermet, CH-1211 Geneva 4, Switzerland
%--------------------------------------------------------------------------


%% Variable handling and default values generation

if nargin==0
    error(sprintf('There\''s no input'))
end
if mod(length(varargin),2)~=0
    error(sprintf('The arguments should be in pairs,in which the first one is a string, as in (data, \''Label\'',{\''A\'',\''B\''})'))
end

stringVars=varargin(1:2:end);
valueVars=varargin(2:2:end);
if any(~cellfun(@isstr,stringVars))
    error(sprintf('The arguments should be in pairs,in which the first one is a string, as in (data, \''Label\'',{\''A\'',\''B\''})'))
end

%data_istable is a logical value used to suppress the default Label value
%assignation

if istable(data)
    data_istable=true;
    [data,Label]=CatTable2Array(data);
    warning('The label will be taken from the table, any label input will be ignored');
else
    data_istable=false;
end

%% String Variables
if ~data_istable
LabelInd=strmatch('Label',stringVars);
if ~isempty(LabelInd) && data_istable==false
    %Check that the Label is a string cell array and has a propper length
    if iscellstr(valueVars{LabelInd}) && length(valueVars{LabelInd})==size(data,2) 
        Label=valueVars{LabelInd};
    else
        error('wrong Label, Label should be a linear cell array of strings, its length should be the same as the number of columns of data')
    end
else
    Label=[];
end
end

Ind=strmatch('Whiskers',stringVars);
if ~isempty(Ind)
    %Check that Whiskers is a string with a valid value
    if ischar(valueVars{Ind}) && any([strcmp('box',valueVars{Ind}),strcmp('none',valueVars{Ind}),strcmp('lines',valueVars{Ind})]) 
        Whiskers=valueVars{Ind};
    else
        error(sprintf('wrong Whiskers, it should be a string saying \''none\'', \''box\'' or \''lines\'''))
    end
else
    Whiskers='box';
end


Ind=strmatch('PointStyle',stringVars);
if ~isempty(Ind)
    %Check that PointStyle is a string
    if ischar(valueVars{Ind}) 
        PointStyle=valueVars{Ind};
    else
        error(sprintf('wrong PointStyle, it should be a string defining the shape of the point \''o\'', \''.\'', \''*\'', etc.'))
    end
else
    PointStyle='o';
end

%% Numeric Variables
Ind=strmatch('Width',stringVars);
if ~isempty(Ind)
    %Check that Width is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        Width=valueVars{Ind};
    else
        error(sprintf('wrong Width, it should be a number'))
    end
else
    Width=0.4;
end

Ind=strmatch('Compression',stringVars);
if ~isempty(Ind)
    %Check that Compression is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        Compression=valueVars{Ind};
    else
        error(sprintf('wrong Compression, it should be a number'))
    end
else
    Compression=5;
end

Ind=strmatch('RangeCut',stringVars);
if ~isempty(Ind)
    %Check that RangeCut is a number
    if isnumeric(valueVars{Ind})
        if size(valueVars{Ind})==1
            RangeCut=valueVars{Ind}*ones(size(data,2),1);
        elseif all(size(valueVars{Ind})== [size(data,2),1]) || all(size(valueVars{Ind}) == [1,size(data,2)])
            RangeCut=valueVars{Ind};
        else
            error(sprintf('wrong RangeCut, it should be a number, or a column/row vector with length=size(data,2)'))
        end
    else
        error(sprintf('wrong RangeCut, it should be a number, or a column/row vector with length=size(data,2)'))
    end
else
    numPoints=sum(~isnan(data));
    RangeCut=(log(numPoints)/log(150)+3*(numPoints).^(1/3))';
end

Ind=strmatch('WhiskersWidthRatio',stringVars);
if ~isempty(Ind)
    %Check that WhiskersWidthRatio is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        WhiskersWidthRatio=valueVars{Ind};
    else
        error(sprintf('wrong WhiskersWidthRatio, it should be a number'))
    end
else
   WhiskersWidthRatio=4;
end

Ind=strmatch('PointSize',stringVars);
if ~isempty(Ind)
    %Check that PointSize is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        PointSize=valueVars{Ind};
    else
        error(sprintf('wrong PointSize, it should be a positive number'))
    end
else
   PointSize=36;
end

Ind=strmatch('LineWidth',stringVars);
if ~isempty(Ind)
    %Check that LineWidth is a number
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind})== [1,1])
        LineWidth=valueVars{Ind};
    else
        error(sprintf('wrong LineWidth, it should be a positive number'))
    end
else
   LineWidth=0.5;
end


%% Mixed Variables
%Some of the following variables are the classic scatter personalizing tools,
%please revisit the scatter documentation for deeper understanding.

Ind=strmatch('MarkerEdgeColor',stringVars);
EdgeIsArray=false;
if ~isempty(Ind)
    %Check that MarkerEdgeColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MarkerEdgeColor=valueVars{Ind};
        EdgeIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MarkerEdgeColor=valueVars{Ind};
    else
        error(sprintf('wrong MarkerEdgeColor input'))
    end
else
    MarkerEdgeColor='k';
end

Ind=strmatch('MarkerFaceColor',stringVars);
%This will be used later when plotting
FaceIsArray=false;

if ~isempty(Ind)
    %Check that MarkerFaceColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MarkerFaceColor=valueVars{Ind};
        FaceIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MarkerFaceColor=valueVars{Ind};
    else
        error(sprintf('wrong MarkerFaceColor input')) 
    end
else
    MarkerFaceColor='flat';
end

Ind=strmatch('MeanColor',stringVars);
MeanColorIsArray=false;
if ~isempty(Ind)
    %Check that MeanColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        MeanColor=valueVars{Ind};
        MeanColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        MeanColor=valueVars{Ind};
    else
        error(sprintf('wrong MeanColor input'))
    end
else
    MeanColor='k';
end

Ind=strmatch('SEMColor',stringVars);
SEMColorIsArray=false;
if ~isempty(Ind)
    %Check that SEMColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        SEMColor=valueVars{Ind};
        SEMColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        SEMColor=valueVars{Ind};
    else
        error(sprintf('wrong SEMColor input'))
    end
elseif strcmp(Whiskers,'box')
    SEMColor=[1 1 1]*0.95;
else
    SEMColor='k';
end

Ind=strmatch('StdColor',stringVars);
StdColorIsArray=false;
if ~isempty(Ind)
    %Check that StdColor takes valid arguments
    if isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) >= [size(data,2),3])
        StdColor=valueVars{Ind};
        StdColorIsArray=true;
    elseif (isnumeric(valueVars{Ind}) && all(size(valueVars{Ind}) == [1,3])) || ischar(valueVars{Ind})
        StdColor=valueVars{Ind};
    else
        error(sprintf('wrong StdColor input'))
    end
elseif strcmp(Whiskers,'box')
    StdColor=[1 1 1]*0.85;
else
    StdColor='k';
end



%% Plotting data
HoldWasOn=ishold;
%This is just a trick to not change anything if there was a plot with hold on already,
%to create a new figure if there was none, and to clear the figure if it
%had hold off
plot(nan,nan)
hold on

xPositions=nan(size(data));
yPositions=data;
for i=1:size(data,2)

%% Sort the data and determine the limits of the different groups in which we divide the points 

% yValues is a column vector containing one of the columns of the matrix
% data, it contains the y values of one of the subset of points we want to
% represent.
% We sort it in order to make groups according to its value, however, we
% want to keep the information of the order of the yvalues, just in case
% its needed for something, we used SortingIndex for that

% eliminate the nans and keep the information to later put it in
% xPositions, this info is kept in yValues_nanIndex
yValues_nanIndex=~isnan(data(:,i));
yValues=data(yValues_nanIndex,i);

[yValues,y_SortingIndex]=sort(yValues);
[~,x_SortingIndex]=sort(y_SortingIndex);

%If one of the columns is empty, we don't plot anything
if isempty(yValues)
    continue
end
% xValues is a column vector as big as yValues, but with all values equal to i. If we
% represented plot(xValues, yValues), we would get all the points with the
% same x,a dot blot, and therefore there would be a massive overlapping
% depending on how many points we had. Thats why the while loops that we
% find afterwards change the x of each point, so they spread and do not
% overlap.
xValues=ones(size(yValues))*i;

%% Plotting the error bars and the standard deviation bars

%Color type managing

if MeanColorIsArray
    MeanIndex=i;
else
    MeanIndex=1;
end

if SEMColorIsArray
    SEMIndex=i;
else
    SEMIndex=1;
end

if StdColorIsArray
    StdIndex=i;
else
    StdIndex=1;
end

%If we want boxes for the mean and std
if strcmp(Whiskers,'box')
    
    yMean=mean(yValues);
    yStd=std(yValues);
    ySem=yStd/sqrt(size(yValues,1));
    yCI=ySem*1.96;
    %plot the standard deviation box
    rectangle('Position',[i-Width/WhiskersWidthRatio,yMean-yStd,2*Width/WhiskersWidthRatio,2*yStd ],'FaceColor',StdColor(StdIndex,:),'EdgeColor', StdColor(StdIndex,:),'LineWidth',0.1);
    %plot the mean+-SEM box
    rectangle('Position',[i-Width/WhiskersWidthRatio,yMean-yCI,2*Width/WhiskersWidthRatio,2*yCI ],'FaceColor',SEMColor(SEMIndex,:),'EdgeColor', SEMColor(SEMIndex,:),'LineWidth',0.1);
    %plot the mean line 
    plot([i-Width/WhiskersWidthRatio i+Width/WhiskersWidthRatio],[yMean yMean],'Color',MeanColor(MeanIndex,:),'LineWidth',2)

end
%% Changing the xValue of each point
% This is done the following way: The points are sorted in equidistant
% groups according to yValues, this distance is stablished as 1/RangeCut of the
% width of the 95% interval of the data. Of course this makes sense for
% more or less normally distributed data, but it can be changed easily, and also it can
% be changed easily changing RangeCut value

range_val=norminv([0.025 0.975],mean(yValues),std(yValues));

%range_val=quantile(yValues,[0.025 0.975]);
cuts=abs(range_val(1)-range_val(2))/RangeCut(i);

% In case one of the variables has equal values for all its points, the
% norminv will return a nan, since std(yValues) will be 0, in that case we
% arbitrarily give cuts the value of cuts=mean(yValues)/2, this value does not
% really matter because there will only be one group anyway
if isnan(cuts)
    cuts=mean(yValues)/2;
end
% The "seed" group contains the values which are in the range of the
% mean+-1/2cuts. cutUp is the higher border of the interval. So we
% take this value, and from it we move up and down in this two loops(steps is a variable 
% that allows us to use the same while loop to move in both directions,
% selecting the subset of yValues that are in each group range, and modifying
% their xValues, so that they are spread and do not overlap.

for steps=[-1,1]
    cutUp=mean(yValues)+cuts/2;
    subsetind= yValues<=cutUp & yValues>cutUp-cuts; 
    keep_going=true;
while keep_going
    % subsetind is a logical variable that represents how many points are
    % in that group
    if all(sum(subsetind)~=[1 0]) %If there's only one point, we represent it in the middle, and xValues remains equal to i
        
        % distmaker is a variable that is equal to Width when sum(subsetind)=1 , and gets smaller when the number of points
        % gets bigger in a group, exponentially, and it tends to zero in
        % the infinite. I don't have strong arguments for the choice of
        % this particular expression, but for me it works very well, and
        % you can tune with Width and Compression very well the appereance of
        % the plot
        distmaker=Width./(exp((sum(subsetind).^2-1)./(sum(subsetind)*Compression)));

        % xSubset is a column vector with all values equal to i, and as
        % long as the subset of the number of values in yValues that belong
        % to this particular group range
        xSubset=ones(sum(subsetind),1)*i;
        
        %oddie is a variable that indicates whether the number of points in
        %the group is odd or even, it's very important for the indexing
        %afterwards
        oddie=mod(size(xSubset,1),2);
        
        % xb is a symmetrical vector with mean=0, and the values in it have
        % the displacements in x that we want to apply to xValues, but if
        % we just applied it increasingly, the smaller yValue would get the most negative 
        % xValues displacement,etc. and what we would have is a series of lines of points, and
        % in publications, what you usually get is that the points with
        % either the highest or
        % lowest values get the positions that are more outside, and then the lowest or 
        % the highest get the central position. In this function, the lowest values of yValues
        % get the external positions, and the highest values get the central positions, making
        % some sort of eyebrow or sad mouth shape)
        % The range of xb gets bigger as the number of points inside the
        % group increases due to the action of distmaker, which gets
        % smaller as number of points increases.
        xb=linspace(1-Width+distmaker,1+Width-distmaker,size(yValues(subsetind),1))-1;
        
        % Since xb is symmetrical it's easy to add the more extreme values
        % of xb to the xValues with the lower yValues. oddie is needed to
        % index properly depending on wether the number of elements is even
        % or odd. If we sorted the values of xb by their absolute value,
        % their index would be (1,end,2,end-1,3,end-2,etc.) If you think of
        % an example or you debug the function and observe the values of xb
        % it's very easy to understand the indexing done here. Also, it's
        % here where oddie is useful
        
        xSubset(1:2:end-oddie)=xSubset(1:2:end-oddie)+xb(1:round(end/2)-oddie)';
        xSubset(2:2:end)=xSubset(2:2:end)-xb(1:round(end/2)-oddie)';
        xValues(subsetind)= xSubset;    
        
        %This following code line is very useful to see how the groups are made,
        %for watching it, comment the scatter lines in the end, and
        %uncomment this line, however groups with one point won't be shown
        %scatter(xValues(subsetind),yValues(subsetind),PointStyle)
    end
    % advance in the while loop as long as the cutUp value is out of the
    % range of the points
    keep_going=~(cutUp>max(yValues) || cutUp<min(yValues));
    cutUp=cutUp+steps*cuts;
    subsetind= yValues<cutUp & yValues>cutUp-cuts;
end
end

%% Drawing each subset of points
% This is to overcome the fact that MarkerEdgeColor and MarkerFaceColor can
% be either strings,1x3 vectors or nx3 arrays

if EdgeIsArray
    EdgeIndex=i;
else
    EdgeIndex=1;
end

if FaceIsArray
    FaceIndex=i;
else
    FaceIndex=1;
end

scatter(xValues,yValues,PointSize,PointStyle,'MarkerEdgeColor', MarkerEdgeColor(EdgeIndex,:),'MarkerFaceColor',MarkerFaceColor(FaceIndex,:),'LineWidth',LineWidth)

%If we want lines to represent the SEM and the std
if strcmp(Whiskers,'lines')
    %plot the mean line
    yMean=mean(yValues);
    plot([i-Width/WhiskersWidthRatio i+Width/WhiskersWidthRatio],[yMean yMean],'Color',MeanColor(MeanIndex),'LineWidth',1.5)
    
    %plot the sd 
    yStd=std(yValues);
    plot([i-Width/WhiskersWidthRatio*0.5 i+Width/WhiskersWidthRatio*0.5],[yMean+yStd yMean+yStd],'Color',StdColor(StdIndex,:),'LineWidth',1.5)
    plot([i-Width/WhiskersWidthRatio*0.5 i+Width/WhiskersWidthRatio*0.5],[yMean-yStd yMean-yStd],'Color',StdColor(StdIndex,:),'LineWidth',1.5)
    plot([i i],[yMean-yStd yMean+yStd],'Color',StdColor(StdIndex,:))
    %plot the conf. interval of the mean
    ySem=yStd/sqrt(size(yValues,1));
    yCI=ySem*1.96;
    plot([i-Width/WhiskersWidthRatio*0.7 i+Width/WhiskersWidthRatio*0.7],[yMean+yCI yMean+yCI],'Color',SEMColor(SEMIndex,:),'LineWidth',1.5)
    plot([i-Width/WhiskersWidthRatio*0.7 i+Width/WhiskersWidthRatio*0.7],[yMean-yCI yMean-yCI],'Color',SEMColor(SEMIndex,:),'LineWidth',1.5)
    plot([i i],[yMean-yCI yMean+yCI],'Color',SEMColor(SEMIndex,:))

end

xPositions(yValues_nanIndex,i)=xValues(x_SortingIndex);
end
set(gca,'xtick',1:size(data,2))

if ~isempty(Label)
set(gca,'xticklabel',Label)
pbaspect([size(data,2)+1,4,1])
end
if ~HoldWasOn
hold off
end
end
