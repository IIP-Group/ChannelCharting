% =========================================================================
% -- Function to compute Sammon's Mapping using FASTA
% -------------------------------------------------------------------------
% -- (c) 2016-2021 Christoph Studer (studer@ethz.ch)
% =========================================================================

function mappedX = SM_fasta(par,d_cov,no_dims)

% random initialization
Y = randn(no_dims,par.U);
x0 = Y(:);

% -- set up fasta parameters
opts = [];
opts.tol = 1e-3;  % use strict tolerance
opts.maxIters = 1000;
opts.recordObjective = true; % record the objective function so we can plot it
opts.verbose=true;
opts.stringHeader='    '; % append a tab to all text output from FISTA.  This option makes formatting look a bit nicer.
opts.accelerate = false;
opts.adaptive = true;

At = @(x) x;
A = @(x) x;

% define weight function
weight_fun = @(x) 1./(x+0.00001); % Sammon's mapping deweighting factor

% define objective funtion f and gradient
f = @(z) f_SM_fun(z,par,d_cov,no_dims,weight_fun);
grad = @(z) grad_SM_fun(z,par,d_cov,no_dims,weight_fun);

% define normalization projection and proximal operator
g = @(x) 0;
prox = @(x,t) proxg_SM_fun(x,par,d_cov,no_dims); %proxg_MDS_fun(x,par,d_cov,no_dims,location); %shrink(x,t*mu);

% call FASTA
[solution, outs] = fasta(A,At,f,grad,g,prox,x0,opts);
mappedX = reshape(solution,no_dims,par.U)';

end

%% auxiliary functions for SM

%% compute objective f function
function objfun = f_SM_fun(x,par,d_cov,no_dims,weight_fun)
Y = reshape(x,no_dims,par.U);

% scaled MDS cost
G = real(Y'*Y);
v = diag(G);
Xpdist = sqrt(max(bsxfun(@plus,v',bsxfun(@plus,(-2)*G,v)),0));
Q = weight_fun(d_cov).*((d_cov-Xpdist).^2);
Q(1:par.U+1:end) = 0;
objfun = sum(Q(:));

if isnan(objfun)||isinf(objfun)
    error('this should not happen');
end

end

%% compute gradient
function grad = grad_SM_fun(x,par,d_cov,no_dims,weight_fun)
Y = reshape(x,no_dims,par.U);

% scaled SM part
G = real(Y'*Y);
v = diag(G);
Xpdist = sqrt(bsxfun(@plus,v',bsxfun(@plus,(-2)*G,v)));
Q = 2*weight_fun(d_cov).*(d_cov./Xpdist-1);
Q(1:par.U+1:end) = 0;

y_grads = Y*(Q-diag(sum(Q, 2)));
grad = y_grads(:);

if any(isnan(grad(:)))||any(isinf(grad(:)))
    error('this should not happen');
end

end

%% compute objective g function
function objfun = g_SM_fun(x,par,d_cov,no_dims)
% do nothing... EPIC
end

%% compute proximal function to centering function
function proxg = proxg_SM_fun(x,par,d_cov,no_dims)
Y = reshape(x,no_dims,par.U);
Y = bsxfun(@minus, Y, mean(Y, 2) ); % center the datapoints
proxg = Y(:);
end