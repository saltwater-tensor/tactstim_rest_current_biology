function [group_spectra] = trial_level_spectral_analysis(data,T,Gamma,options_mt)
    
%% This section has been commented out due to MEMORY ISSUES 
% A single trial produces a 30 ~ 40 MB spectral output and for 4000 plus
% trials the output will be approximately 200GB
    % Calculating on a trial by trial basis
    N = size(data,1);
    trial_spectra = cell(N,1);
    d = length(options_mt.embeddedlags) - 1;
    acc = 0 ;
    N = size(data,1);
    
    for n = 1:N
        h = waitbar(n/N);
        X = data{n};
        gamma = Gamma(acc + (1:(sum(T{n})-length(T{n})*d)),:);
        acc = acc + size(gamma,1);
        trial_spectra{n} = hmmspectramt(X,T{n},gamma,options_mt);
        trial_spectra{n}.state = rmfield(trial_spectra{n}.state,'ipsd');
        trial_spectra{n}.state = rmfield(trial_spectra{n}.state,'pcoh');
        trial_spectra{n}.state = rmfield(trial_spectra{n}.state,'phase');
        disp(['Trial ' num2str(n)])
        save
    end
    
%% Group level
    
    % Calculating on the group level
%     T = cell2mat(T');
%     % T = T';
%     data = cell2mat(data);
%     group_spectra = hmmspectramt(data,T,Gamma,options_mt);
    
end