function [group_spectra] = single_trial_spectral_analysis(data,T,Gamma,options_mt,trial_name)
    
%% This section has been commented out due to MEMORY ISSUES 
% A single trial produces a 30 ~ 40 MB spectral output and for 4000 plus
% trials the output will be approximately 200GB
    % Calculating on a trial by trial basis

    d = length(options_mt.embeddedlags) - 1;
    acc = 0 ;
    
    X = data;
    clear data
    gamma = Gamma(acc + (1:(sum(T)-length(T)*d)),:);
    acc = acc + size(gamma,1);
    trial_spectra = hmmspectramt(X,T,gamma,options_mt);
    trial_spectra.state = rmfield(trial_spectra.state,'ipsd');
    trial_spectra.state = rmfield(trial_spectra.state,'pcoh');
    trial_spectra.state = rmfield(trial_spectra.state,'phase');
    
    save(trial_name,"trial_spectra")
    
    
%% Group level
    
    % Calculating on the group level
%     T = cell2mat(T');
%     % T = T';
%     data = cell2mat(data);
%     group_spectra = hmmspectramt(data,T,Gamma,options_mt);
    
end