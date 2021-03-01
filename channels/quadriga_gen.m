% =========================================================================
% -- Script to Generate Quadriga Channel Vectors
% -------------------------------------------------------------------------
% -- (c) 2016-2021 Christoph Studer (studer@ethz.ch)
% =========================================================================

function [H] = quadriga_gen(par,location)

% *** requires the quadriga toolbox to be in MATLAB's path ***
addpath(genpath('quadriga'));

%quadriga parameters
par.channel.BW = 10e6; % Bandwidth in MHz
par.channel.Rxdistance = 200; % max user distance to array in meters
par.channel.fc = 2000e6; % carrier frequency in Hz
par.channel.BW = 10e6; % Bandwidth in MHz
par.channel.Rxdistance = 2000; % max user distance to array in meters
%par.channel.scenario = 'BERLIN_UMa_LOS'; % LoS scenario
par.channel.scenario = 'BERLIN_UMa_NLOS'; % non-LoS scenario
par.channel.N = 1; %subcarriers

s = qd_simulation_parameters;
s.center_frequency = par.channel.fc; % channel frequency
s.sample_density = 2; % number of (spatial) samples per lambda/2. **** WHAT DOES THAT MEAN?
s.use_absolute_delays = 1; 

l = qd_layout(s);
l.no_rx = par.U; %number of users

%% place users
% generate random users in square
l.rx_position(1,:) = par.ue.x; % *** assuming first coordinate is x
l.rx_position(2,:) = par.ue.y; % *** assuming second coordinate is y
l.rx_position(3,:) = par.ue.z; % *** assuming third coordinate is z

% one track is sufficient
for i=1:l.no_rx
    l.track(i).no_snapshots = 1;
end

% set scenario
l.set_scenario(par.channel.scenario);

disp('Antenna placement');

l.tx_array.generate( 'omni' );  
l.tx_array.no_elements = par.B;

ant_pos_mat = zeros(3,l.tx_array.no_elements);
lambda = s.wavelength;
% positioning of array elements 
for i=1:l.tx_array.no_elements
    % only on y-axis
    ant_pos_mat(2,i) = lambda/2*(i) - lambda/2*(l.tx_array.no_elements)/2;
end

l.tx_array.element_position = ant_pos_mat;
l.tx_position(3)= 10; % antenna array is 10m above ground
l.rx_array.generate('omni');  

% shows the position of users and antennas (sanity check)
disp('visualization');
figure(2)
l.visualize([],[],0);
view(-33, 60);
drawnow;

disp('channel generation');

p =l.init_builder;
p.plpar =[];

p.scenpar.SC_lambda = 0;
p.scenpar.NumClusters = 1; % number of paths
p.gen_ssf_parameters; %small scaled fading

% get channel coefficients
c = l.get_channels;

coeff = squeeze( cat( 1, c.coeff ) );
delay = permute( cat(3,c.delay) , [3,1,2] );

% % Visualizing the clusters
% cb.visualize_clusters(1);

%% Calculating the frequency response for this channel
disp('started frequency response');

bandwidth = par.channel.BW/par.channel.N;
subcarriers = par.channel.N;
i_snapshot = 1;

% Using built-in method to evaluate the frequency response from the channel
% coefficients
%freq_response = c.fr(bandwidth,subcarriers,i_snapshot);
for ll = 1:l.no_rx
    freq_response(ll,:,:)=c(ll).fr(bandwidth,subcarriers,i_snapshot); % DFT of each channel vector
end
% Writing the frequency response in the desired way: a matrix corresponding
% to the user grid, each point with a vector with the number of antennas
H = zeros(par.B,par.U,par.channel.N);

for uu=1:par.U
    H(:,uu,:) = freq_response(uu,:,:);
end

% save channel plus locations
save(['./channels/',par.channel.scenario,'.mat'], 'H', 'location');

end