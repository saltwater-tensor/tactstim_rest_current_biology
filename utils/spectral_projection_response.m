function [state] = spectral_projection_response(response_spectra, supply_profiles)

% Of course the code can be far more efficient 
% But this version is easier to understand for people who would use it for
% the first time

coh = response_spectra.state(1).coh;
ndim = size(coh,2);
ndim2 = ndim*(ndim-1)/2;
Nf = size(coh,1);
Ncomps = size(supply_profiles,2);
K = size(response_spectra.state,2);
ind_offdiag = triu(true(ndim),1)==1;
coh_comps = zeros(1,K,Nf,ndim,ndim);

for k = 1: size(response_spectra.state,2)

    coh_comps(1,k,:,:,:) = response_spectra.state(k).coh;

    for rgn = 1:42
        psd_comps(1,k,:,rgn) = response_spectra.state(k).psd(:,rgn,rgn);
    end
%     coh = response_spectra.state(k).coh;


%     for rgn1 = 1:42
%         for rgn2 = 1:42
%             c = coh(:,rgn1,rgn2);
%             c_proj = c'*supply_profiles;
%             coh_proj(:,rgn1,rgn2) = c_proj;
%         end
%     end
% 
%     Nan_mat = NaN(size(coh,2),1);
%     Nan_mat = diag(Nan_mat);
% 
%     for cpr = 1:4
%         cprj = squeeze(coh_proj(cpr,:,:));
%         cprj = cprj + Nan_mat;
%         cprj(isnan(cprj)) = 1.00;
%         coh_proj(cpr,:,:) = cprj;
%     end
    
end
    
    %PSD
    n = 1;
    Xpsd = zeros(Nf,ndim*K);
    for k=1:K
        ind = (1:ndim) + (k-1)*ndim;
        Xpsd(:,ind)= squeeze(abs(psd_comps(n,k,:,:)));
    end

    opt = statset('maxiter',1);
    [~,b] = nnmf(Xpsd,Ncomps,'algorithm','als',...
        'w0',supply_profiles,'Options',opt);
    psd_b = b'; 
    clear b
    
    % Coherence
    n = 1;
    Xcoh = zeros(Nf,K*ndim2);
    for k = 1:K
        ind = (1:ndim2) + (k-1)*ndim2;
        ck = squeeze(abs(coh_comps(n,k,:,:,:)));
        Xcoh(:,ind) = ck(:,ind_offdiag);
    end

    opt = statset('maxiter',1);
    [~,b] = nnmf(Xcoh,Ncomps,'algorithm','als',...
    'w0',supply_profiles,'Options',opt);
    coh_b = b';

    
    % Reshape stuff
    for k = 1:K
        state(k).psd= zeros(Ncomps,ndim,ndim);
        state(k).coh = ones(Ncomps,ndim,ndim);
        ind = (1:ndim) + (k-1)*ndim;
        for i = 1:Ncomps
            state(k).psd(i,:,:) = diag(psd_b(ind,i));
        end
        ind = (1:ndim2) + (k-1)*ndim2;
        for i = 1:Ncomps
            graphmat = zeros(ndim);
            graphmat(ind_offdiag) = coh_b(ind,i);
            graphmat=(graphmat+graphmat') + eye(ndim);
            state(k).coh(i,:,:) = graphmat;
        end
    end


end
