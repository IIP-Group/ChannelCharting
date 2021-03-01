% =========================================================================
% -- Function to Generate Line-of-Sight Channel Vectors
% -------------------------------------------------------------------------
% -- (c) 2016-2016 Christoph Studer (studer@ethz.ch)
% =========================================================================

function H = channel_los(par)

% evaluate all pairwise distances with spherical wave model
dist_swm = sqrt((par.bs.x*ones(1,par.U)-ones(par.B,1)*par.ue.x').^2 + ...
    (par.bs.y*ones(1,par.U)-ones(par.B,1)*par.ue.y').^2 + ...
    (par.bs.z*ones(1,par.U)-ones(par.B,1)*par.ue.z').^2);

% simple distance scaling
H_swm = 1./dist_swm .* exp(-1i*2*pi*dist_swm/par.lambda);

% output results
H = H_swm;