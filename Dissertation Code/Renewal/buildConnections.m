function wholeNetwork = buildConnections(wholeNetwork, ffConnect, ...
    delayRange, numConnectionsPerNeuron, weightLowerbound, ...
    weightUpperbound, weightMultiplicand, isRandomWeights)
% function wholeNetwork = buildConnections(wholeNetwork, ffConnect, ...
%     baseDelay, numConnectionsPerNeuron, weightLowerbound, weightUpperbound)
%
% build feed forward connections from nuclei ffConnect(1) to ffConnect(2)
% with the baseDelay minimum conductance delay and numConnectionsPerNeuron
% connections per neuron in nuclei ffConnect(1).
%
% weightLowerbound is the lower bound to be enforced for the connection weights,
% 		   corresponds to the inhibitory minimum.
% weightUpperbound is the upper bound to be enforced for the connection weights,
% 		   corresponds to the excitatory maximum.

% connectivityMatrix is a matrix of zeros and ones, with ones implying the
% existence of a connection. It is used to impose the network connectivity
% on any matrix of equal size.
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012
% Further modified by Jeffrey Rodny jrodny@ucmerced.edu 2013-2017
% Here, Jeff Rodny modified only those labeled: 'JJ CODE'

%% JJ CODE:
%% If its: hippocampus bins to bins
 % loop through for each subseciton of INTh all-to-all connect each
 %    subsection of HIPP to the correpsonding subsestionc in INTh
if wholeNetwork.nucIndex.HIPP == ffConnect(1) ...
        && wholeNetwork.nucIndex.INTh == ffConnect(2)
    params1 = wholeNetwork.nuclei{ffConnect(1)}{ffConnect(1)}.params;
    params2 = wholeNetwork.nuclei{ffConnect(2)}{ffConnect(2)}.params;
    
    
    
    connectivityMatrix = zeros(params1.N, params2.N);
    pSize = size(connectivityMatrix);
    segmentHSize = floor(pSize(1)/15);
    segmentISize = floor(pSize(2)/15);
    for kk = 0:14
        for ii = (kk*segmentHSize)+1:(kk*segmentHSize)+segmentHSize
            for jj = (kk*segmentISize)+1:(kk*segmentISize)+segmentISize
                connectivityMatrix(ii,jj) = 1;
            end
        end
    end
   
    S=[params1.initExWeight*ones(params1.Ne,params2.N); ...
       params1.initInWeight*ones(params1.Ni,params2.N)];

    S = S .* connectivityMatrix .* weightMultiplicand + (isRandomWeights .* (rand(size(S))-.5));

    for kk = [1,4,7,8,12];
        for ii = (kk*segmentHSize)+1:(kk*segmentHSize)+segmentHSize
            for jj = (kk*segmentISize)+1:(kk*segmentISize)+segmentISize
                S(ii,jj) = 10;
            end
        end
    end
    
    %% or if it's SEN to INT then connect to correct INTh subsections
elseif wholeNetwork.nucIndex.INTh == ffConnect(2) ...
        && (wholeNetwork.nucIndex.SENa == ffConnect(1) || ...
        wholeNetwork.nucIndex.SENb == ffConnect(1) || ...
        wholeNetwork.nucIndex.SENc == ffConnect(1))
    
    % Hippocampus subsections for reference:
    %     HippStim.BL.index = 1;
    %     HippStim.BF.index = 2;
    %     HippStim.BRa.index = 3;
    %     HippStim.BRb.index = 4;
    %     HippStim.LF.index = 5;
    %     HippStim.LRa.index = 6;
    %     HippStim.LRb.index = 7;
    %     HippStim.FRa.index = 8;
    %     HippStim.FRb.index = 9;
    %     HippStim.RaRb.index = 10;
    %     HippStim.TB.index = 11;
    %     HippStim.TL.index = 12;
    %     HippStim.TF.index = 13;
    %     HippStim.TRa.index = 14;
    %     HippStim.TRb.index = 15;
    % Bell:  1,2,3,4,11
    % Light: 1,5,6,7,12
    % Food:  2,5,8,9,13

    if wholeNetwork.nucIndex.SENa == ffConnect(1)
        tmpChannels = [1,2,3,4,11] - 1;
    elseif wholeNetwork.nucIndex.SENc == ffConnect(1)
        tmpChannels = [1,5,6,7,12] - 1;
    elseif wholeNetwork.nucIndex.SENb == ffConnect(1)
        tmpChannels = [2,5,8,9,13] - 1;
    end
    
    params1 = wholeNetwork.nuclei{ffConnect(1)}{ffConnect(1)}.params;
    params2 = wholeNetwork.nuclei{ffConnect(2)}{ffConnect(2)}.params;
    
    
    
    connectivityMatrix = zeros(params1.N, params2.N);
    pSize = size(connectivityMatrix);
    segmentHSize = floor(pSize(1)/15);
    segmentISize = floor(pSize(2)/15);
    for kk = tmpChannels
        for ii = 1:params1.N
            for jj = (kk*segmentISize)+1:(kk*segmentISize)+segmentISize
                connectivityMatrix(ii,jj) = 1;
            end
        end
    end
    %% loop through for each subseciton of INTh all-to-all connect each subsection of HIPP to the correpsonding subsestionc in INTh
    
    S=[params1.initExWeight*ones(params1.Ne,params2.N); ...
       params1.initInWeight*ones(params1.Ni,params2.N)];

    S = S .* connectivityMatrix .* weightMultiplicand + (isRandomWeights .* (rand(size(S))-.5));

    for kk = [1,4,7,8,12];
        for ii = (kk*segmentHSize)+1:(kk*segmentHSize)+segmentHSize
            for jj = (kk*segmentISize)+1:(kk*segmentISize)+segmentISize
                S(ii,jj) = 10;
            end
        end
    end
    
   %% END JJ CODE 
%% Original code: for all other connections:
else

    params1 = wholeNetwork.nuclei{ffConnect(1)}{ffConnect(1)}.params;
    params2 = wholeNetwork.nuclei{ffConnect(2)}{ffConnect(2)}.params;


    connectivityMatrix = zeros(params1.N, params2.N);
    for i = 1:params1.N
        connectionsSoFar = [];
        for j = 1:numConnectionsPerNeuron
            if i <= params1.Ne
                randomConnection = ceil(rand*params2.N);
            else
                randomConnection = ceil(rand*params2.Ne);
            end
            while nnz(connectionsSoFar == randomConnection) > 0
                if i <= params1.Ne
                    randomConnection = ceil(rand*params2.N);
                else
                    randomConnection = ceil(rand*params2.Ne);
                end
            end
            %% JJ CODE EDIT (Discovered and mentioned to ben: nuclei with 0
            %%   excitatory connections)
            if randomConnection == 0
                continue;
            end
            %% end JJ EDIT
            connectionsSoFar = [connectionsSoFar randomConnection];
            connectivityMatrix(i,randomConnection) = 1;
        end
    end
clear randomConnection
clear connectionsSoFar

% with uniform weights.
S=[params1.initExWeight*ones(params1.Ne,params2.N); ...
   params1.initInWeight*ones(params1.Ni,params2.N)];
%% JJ CODE Added "weightMultiplicand"
% S = S .* connectivityMatrix;
S = S .* connectivityMatrix .* weightMultiplicand + (isRandomWeights .* (rand(size(S))-.5));
%% End JJ code
end

%Izhi-Szat paper has short-term STDP values which scale input for each
%neuron, unique to each synapse.
%    gatingValues = zeros(size(S));

% The millisecond delays between spikes.
% The range of the millisecond delays between spikes, uniformly chosen.
delayRange = delayRange(1):delayRange(2);
conductanceDelays = delayRange(ceil(rand(params1.N,params2.N)*length(delayRange)));

% this line normally sets all inhibitory connections to delay 1, as in the
% izhikevich paper. However, for nuclei to nuclei we expect atleast a
% minimum delay of baseDelay+1+1 (the original 1, but for each nuclei).
if params1.N > params1.Ne
	conductanceDelays((params1.Ne+1):params1.N,:) = delayRange(1);
end
conductanceDelays = conductanceDelays .* connectivityMatrix;

% Initially, there are NO conducting potentials (assumption).
conductingPotentials = [];

wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.conductingPotentials = conductingPotentials;
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.conductanceDelays = conductanceDelays;
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.S = S;
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.connectivityMatrix = connectivityMatrix;
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.lastFire = ...
    NaN(wholeNetwork.nuclei{ffConnect(1)}{ffConnect(1)}.params.N, ...
    wholeNetwork.nuclei{ffConnect(2)}{ffConnect(2)}.params.N);
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.params.weightLowerbound = weightLowerbound;
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.params.weightUpperbound = weightUpperbound;

%% JJ code
% reset the eligibility traces
wholeNetwork.nuclei{ffConnect(1)}{ffConnect(2)}.eligibilityTrace = zeros(size(connectivityMatrix));
%%


end
