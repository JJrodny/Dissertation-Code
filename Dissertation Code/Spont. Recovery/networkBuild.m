function network = networkBuild(params)
% This function will take in parameters and output the data structures
% required to run a network.
% network parameters.
% %Number of excitatory neurons
% params.Ne = 800;
% %Number of inhibitory neurons
% params.Ni = 200;
% params.N = params.Ne + params.Ni;
% params.exA = 0.02;
% params.exB = 0.2;
% params.exC = -65;
% params.exD = 8;
% params.inA = 0.1;
% params.inB = 0.2;
% params.inC = -65;
% params.inD = 2;
% % Number of outgoing connections per neuron.
% params.numConnectionsPerNeuron = 100;
% % Maximum delay amount. Izhi&Szat have 20ms max, uniform.
% params.delayRange = [1 20];
% params.initExWeight = 5;
% params.initInWeight = -4;
% % weight bounds, applied every millisecond.
% params.weightUpperbound = 8;
% params.weightLowerbound = -8;
% % Izhikevich iters two .5 ms to add to 1ms. This is not really a parameter.
% params.timeStep = .5; % 0.5 ms.

% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012
% Further modified by Jeffrey Rodny jrodny@ucmerced.edu 2013-2017
% Here, Jeff Rodny modified only those labeled: 'JJ CODE'


params.a=[params.exA*ones(params.Ne,1);    params.inA*ones(params.Ni,1)];
params.b=[params.exB*ones(params.Ne,1);    params.inB*ones(params.Ni,1)];
params.c=[params.exC*ones(params.Ne,1);    params.inC*ones(params.Ni,1)];
params.d=[params.exD*ones(params.Ne,1);    params.inD*ones(params.Ni,1)];


% connectivityMatrix is a matrix of zeros and ones, with ones implying the
% existence of a connection. It is used to impose the network connectivity
% on any matrix of equal size.
connectivityMatrix = zeros(params.N);
for i = 1:params.N
    connectionsSoFar = [];
    for j = 1:params.numConnectionsPerNeuron
        if i <= params.Ne || params.Ne == 0
            randomConnection = ceil(rand*params.N);
        else
            randomConnection = ceil(rand*params.Ne);
        end
        while nnz(connectionsSoFar == randomConnection) > 0
            if i <= params.Ne || params.Ne == 0
                randomConnection = ceil(rand*params.N);
            else
                randomConnection = ceil(rand*params.Ne);
            end
        end
        connectionsSoFar = [connectionsSoFar randomConnection];
        connectivityMatrix(i,randomConnection) = 1;
    end
end
clear randomConnection
clear connectionsSoFar

%% Jeff Code:
 % Implemented 'lateral inhibition'
if ...
        params.numConnectionsPerInhibNeuron == params.Ne && ...
        params.numConnectionsPerExciteNeuron == params.Ni
    display([params.Name, ' implemented with lateral inhibition.']);
    for ii = 1:params.N
        if ii <= params.Ne % if it's excitatory
            for jj = params.Ne+1:params.N
                connectivityMatrix(jj,ii) = 1;% * round(rand);
            end
        else % otherwise it's inhibitory
            for jj = 1:params.Ne
                connectivityMatrix(jj,ii) = 1;
            end
%             for jj = params.Ne+1:params.N
%                 connectivityMatrix(jj,ii) = 1;% * round(rand);
%             end
        end
    end
    
end

%%


% with uniform weights.
network.S=[params.initExWeight*ones(params.Ne,params.N);params.initInWeight*ones(params.Ni,params.N)];
network.S = network.S .* connectivityMatrix; %.* (rand(size(connectivityMatrix))*.1-.05+1);

%Izhi-Szat paper has short-term STDP values which scale input for each
%neuron, unique to each synapse.
%    gatingValues = zeros(size(S));

% The range of the millisecond delays between spikes, uniformly chosen.
delayRange = params.delayRange(1):params.delayRange(2);
conductanceDelays = delayRange(ceil(rand(params.N)*length(delayRange)));
% set inhibitory connections to minimum delay.
if params.N > params.Ne
	conductanceDelays((params.Ne+1):params.N,:) = 1;
end
network.conductanceDelays = conductanceDelays .* connectivityMatrix;

% Initially, there are NO conducting potentials (assumption).
network.conductingPotentials = [];

%% JJ Code
network.eligibilityTrace = zeros(size(connectivityMatrix));

%% END JJ

% Initial neuron model values.
network.membranePotential=params.exC*ones(params.Ne+params.Ni,1);
network.equalibriumForce=params.b .* network.membranePotential;
network.firings=[]; % spike timings

% Initialize all neurons to have their last fire over 50 timesteps ago.
network.lastFire = NaN(params.N);
network.connectivityMatrix = connectivityMatrix;
network.params = params;

end