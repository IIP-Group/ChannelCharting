% =========================================================================
% -- Simple Channel Charting Simulator
% -------------------------------------------------------------------------
% -- (c) 2016-2021 Christoph Studer, Emre Gonultas, and Said Medjkouh
% -- e-mail: studer@ethz.ch, eg566@cornell.edu, sm2685@cornell.edu
% -------------------------------------------------------------------------
%
% -- If you use this simulator or parts of it, then you must cite our
% -- original Channel Charting paper:
% 
% -- Christoph Studer, Said Medjkouh, Emre Gonultas, Tom Goldstein, Olav
% -- Tirkkonen, "Channel charting: Locating users within the radio 
% -- environment using channel state information," IEEE Access, Vol. 6,
% -- pp. 47862-47598, Aug. 2018
%
% -- Erratum: The above paper claims to be using the plane-wave model 
% -- for the vanilla line-of-sight channels; this is incorrect and
% -- a spherical wave model was used in our simulations. 
%
% -- In case you decide to use the Sammon's mapping solver, then you must 
% -- cite the original FASTA paper:
%
% -- Tom Goldstein, Christoph Studer, and Richard G. Baraniuk, "A field 
% -- guide to forward-backward splitting with a FASTA implementation," 
% -- Technical Report, arXiv preprint:1411.3406, Nov. 2014; available
% -- at https://arxiv.org/pdf/1411.3406.pdf
%
% -- For more details, please watch the following YouTube video:
% -- https://www.youtube.com/watch?v=dQw4w9WgXcQ
%
% =========================================================================

function channel_charting

clc
close all
addpath('include');
addpath('channels');

% initialize random seed (make results reproducible)
rng(1); % results depend on choice of random seed

%% set up simulation parameters

% physics
par.c = 3e8;    % speed of light [m/s]
par.fc = 2e9;   % carrier frequency [Hz]
par.lambda =  par.c/par.fc; % carrier wavelength [m]

% system specification
par.avg = 10; % number of channel averages
par.SNRdB = 0; % signal to noise ratio in dB
par.B = 32; % number of basestation (BS) antennas
par.U = 2048; % number of user-equipment locations
par.plot = true; % plot scenario?
par.channel.model = 'LoS'; % 'LoS' = vanilla LoS; 'QLoS' = quadriga LoS; 'QNLoS' = quadriga nLoS

% output settings
par.plotMapping = true;  % plot mapping in 3D or 2D
par.boundary = true; % highlight boundary

% feature scaling parameter
par.sigma = 16; % decay factor in attenuaten term: beta = 1+1/alpha

% channel charting method 'PCA' = PCA; 'SM' = Sammon's mapping
DRmethodList = {'PCA','SM'};

%% determine antenna and UE locations

% generate uniform linear antenna (ULA) array at origin
par.bs.x = ones(par.B,1)*0; % no change in x axis
par.bs.y = par.lambda/2*((1:par.B)'-(1+par.B)/2); % linear array in y
par.bs.z = ones(par.B,1)*10; % antenna array is 10m above ground

% generate random UE positions in square
par.ue.x = 100+rand(par.U,1)*500; % between 100m and 600m of array
par.ue.y = (2*rand(par.U,1)-1)*500; % between 500m and -500m or array
par.ue.z = ones(par.U,1)*1.5; % 1.5m above ground

%load pre-defined trajectory
load('./include/vip_drawing.mat');
par.Circ = size(coords,1);
par.ue.x(1:par.Circ) = -coords(:,2)*2.4+950;
par.ue.y(1:par.Circ) = coords(:,1)*2.4-560;
par.Circ = 234;
par.shape = 0;
par.coordShape = par.Circ+1:par.Circ+par.shape;

%collect all coordinates in a matrix
location = [ par.ue.x' ; par.ue.y' ; par.ue.z' ] ;

%% channel vector generation

% one could generate new quadriga channels here 
% but this requires quadriga to be installed
% call this here if everything is set up properly:
% H = quadriga_gen(par,location)

% generate channel vectors
switch par.channel.model
    case 'LoS' % vanilla line-of-sight channel
        H = channel_los(par);
    case 'QLoS' % pregenerated Quadriga line-of-sight channel
        % Important: 
        %  - works only for U=2048
        %  - overwrites location 
        load('./channels/BERLIN_UMA_LOS'); 
    case 'QNLoS' % pregenerated Quadriga non-line-of-sight channel
        % Important: 
        %  - works only for U=2048
        %  - overwrites location         
        load('./channels/BERLIN_UMA_NLOS'); % works only for U=2048
    otherwise
        error('par.channel.model undefined');
end

%% define properties for plotting/visualization

% define colormap 
color1 = (par.ue.x-min(par.ue.x))/max(par.ue.x) ;
color2 = (par.ue.y-min(par.ue.y))/max(par.ue.y)  ;
color3 = zeros(length(color2),1) ;
par.colorMap = [color2 color1 color3] ;

% find boundary of UE locations
if (par.boundary)
    par.boundIdx = boundary(par.ue.x,par.ue.y);
end

% show simulation scenario
if par.plot
    plotMapping([],par,'channel')
    drawnow; % force plot scenario
end

%% Feature extraction

% normalize channel gains to average norm 1 per entry to match SNR definition
H = (sqrt(par.B*par.U)/norm(H,'fro'))*H;

% compute covariance matrices & beamspace transform
D = zeros(par.U,par.B*par.B);
for uu=1:par.U
    
    % averaging over par.avg samples
    N = sqrt(10^(-par.SNRdB/10)*0.5)*(randn(par.B,par.avg)+1i*randn(par.B,par.avg));
    HN = bsxfun(@plus,H(:,uu),N);    
    Hmean = mean(HN,2);
        
    % feature scaling
    v = (par.B^(1/(2*par.sigma))/norm(Hmean,2)^(1+1/(2*par.sigma)))*Hmean;
    
    % beamspace transform & take the absolute value
    K = abs(fft2(v*v')/sqrt(par.U^2));
    
    % vectorize CSI features
    D(uu,:) = K(:);
    
end

%% perform channel charting via different methods

for methIdx = 1:length(DRmethodList)
    
    par.DRmethod = char(DRmethodList(methIdx)); % extract current method        
    no_dims = 2; % all methods use 2-dimensional embeddings    
    par.nameParam = [par.DRmethod '_' par.channel.model]; % name to be used when naming figures
    
    %% choose a method and its parameters
    switch par.DRmethod
        case 'PCA' % principal component analysis
            disp('Compute PCA...');
            mappedX = pca(D.', 'NumComponents', no_dims,'Algorithm','eig'); % use MATLAB's PCA function
                        
        case 'SM' % Sammon's mapping using FASTA
            disp('Calculating pairwise l2-distances between features...');            
            sum_D = sum(D .^ 2, 2);
            d_cov = bsxfun(@plus, sum_D, bsxfun(@plus, sum_D', -2 * (D * D')));
            d_cov = sqrt(d_cov); % extract distances
            d_cov = d_cov/sqrt(var(d_cov(:))); % normalize
            
            disp('Compute Sammon''s mapping with FASTA...');
            mappedX = SM_fasta(par,d_cov,no_dims);
            
        otherwise
            error('Invalid par.DRmethod selected')
    end
    
    % plot and save resulting channel chart 
    if (par.plotMapping)
        plotMapping(mappedX,par,'embedding2D')
    end
    
    % calculate performance metrics and plot TW and CT 
    TW_CT(location,mappedX,par)
    
end