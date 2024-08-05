function [options] = format_options(options)

% data options
if ~isfield(options,'Fs'), options.Fs = 1; end
if ~isfield(options,'onpower'), options.onpower = 0; end
if ~isfield(options,'leida'), options.leida = 0; end
if ~isfield(options,'embeddedlags'), options.embeddedlags = 0; end
if ~isfield(options,'pca'), options.pca = 0; end
if ~isfield(options,'pca_spatial'), options.pca_spatial = 0; end
if ~isfield(options,'rank'), options.rank = 0; end
if ~isfield(options,'firsteigv'), options.firsteigv = 0; end
if ~isfield(options,'varimax'), options.varimax = 0; end
if ~isfield(options,'pcamar'), options.pcamar = 0; end
if ~isfield(options,'pcapred'), options.pcapred = 0; end
if ~isfield(options,'vcomp') && options.pcapred>0, options.vcomp = 1; end
if ~isfield(options,'filter'), options.filter = []; end
if ~isfield(options,'detrend'), options.detrend = 0; end
if ~isfield(options,'downsample'), options.downsample = 0; end
if ~isfield(options,'leakagecorr'), options.leakagecorr = 0; end
if ~isfield(options,'standardise'), options.standardise = 1; end
if ~isfield(options,'As'), options.As = []; end
if ~isfield(options,'standardise_pc') 
    options.standardise_pc = length(options.embeddedlags)>1; 
end



end