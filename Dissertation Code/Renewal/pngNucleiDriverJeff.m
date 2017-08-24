% This is a driver file to build and run a many nuclei neural network.
% It uses the functions...
%    networkBuild(params)
%    buildConnections(wholeNetwork, ffConnect, baseDelay, numConnectionsPerNeuron, weightLowerbound, weightUpperbound)
% and the scripts...
%     pngSimStimuliInit
%   networkIterate  ... which uses function
%       STDP(lastFire1, fired, t, Ne2, lastFire2)
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012
% Further modified by Jeffrey Rodny jrodny@ucmerced.edu 2013-2017

%% to Run the code in parallel
% if matlabpool('size') == 0 % checking to see if my pool is already open
%     matlabpool(6)
% end

%% Actual code
clear all;
close all;

currentTime = now;
wholeNetwork.expId = datestr(currentTime,'mm-dd HH.MM.SS.FFF');
clear currentTime;
disp(wholeNetwork.expId);
wholeNetwork.t = 0;
wholeNetwork.firings = [];

%% %% Parameter exploration: tmpSENINTMult, tmpMult, tmpUpp, tmpPFCSTR, tmpInitHipp

%% Strength of SEN to INT
% [.5, 1.5, 2.5, 3.5]
wholeNetwork.tmpSENINTMult = .5 + 1 * floor(rand() * 4);

%% Strength of Hipp to SNC connection (a non-learning connection)
% % [.04, .08, .12, .16]
wholeNetwork.tmpMult = 0.04 + 0.04 * floor(rand() * 4);

%% Hipp to INTHipp connection upper bound on connection strength
   % (a learning connection/weight)
% [12, 15, 18, 21]
wholeNetwork.tmpUpp = 12 + 3 * floor(rand() * 4);

%% PFC to STR connection upper bound on connection strength
   % (a learning connection/weight)
% [2, 5, 8, 11]
wholeNetwork.tmpPFCSTR = 2 + 3 * floor(rand() * 4);

%% Initial Strength of unlearned HIPP to INTHIPP connections
% [ 0, 0.05, 0.10, 0.15]
wholeNetwork.tmpInitHipp = 0.0 + 0.05 * floor(rand() * 4);


%% Set the number of learning and extinction trials
% 100 learning trials
wholeNetwork.Num_Learning_Trials = 100;

% 100 BL learning trials
wholeNetwork.Num_Extinction_Trials = 100;
wholeNetwork.Extinction_Trials_in_A = 1;%min(1,0 + 1 * floor(rand() * 2));

% 100 back and forth trials of just B and just L renewal trial of light only (in A or B)
wholeNetwork.Num_Renewal_Trials = 100;
wholeNetwork.Renewal_Trials_in_A = 0;
% 

% 
% % 100 post extinction trials (in A or B)
% wholeNetwork.Num_Post_Ext_Trials = 50;
% wholeNetwork.Post_Ext_Trials_in_A = min(1,0 + 1 * floor(rand() * 2));
% % 
% % 70% of the trials should all be in same context
% if rand() > 0.3
%     wholeNetwork.Post_Ext_Trials_in_A = 1;
%     wholeNetwork.Renewal_Trials_in_A = 1;
%     wholeNetwork.Ext_Trials_in_A = 1;
% end

% ' H-iHmax=',num2str(wholeNetwork.tmpUpp),...
% '_H-iHinit=',num2str(wholeNetwork.tmpInitHipp),...
% '_iH-DA=',num2str(wholeNetwork.tmpMult),...
% '_Sen-INT=',num2str(wholeNetwork.tmpSENINTMult),...
% '_PFC-STR=',num2str(wholeNetwork.tmpPFCSTR),...
% '_ExTrA=',num2str(wholeNetwork.Ext_Trials_in_A),...
% '_ReTrA=',num2str(wholeNetwork.Renewal_Trials_in_A),...
% '_PExTrA=',num2str(wholeNetwork.Post_Ext_Trials_in_A)


wholeNetwork.BorLbool = 0;%floor(rand() * 2);
% wholeNetwork.Post_Ext_Trials_in_A = 1;
% wholeNetwork.Renewal_Trials_in_A = 1;
% wholeNetwork.Extinction_Trials_in_A = 1;
wholeNetwork.tmpSENINTMult = 1.5;
wholeNetwork.tmpInitHipp = .15;
wholeNetwork.tmpPFCSTR = 11;
wholeNetwork.tmpUpp = 21;
wholeNetwork.tmpMult = .08;
wholeNetwork.tmpSENINTMult = 1.5



% reset DA1, DA2 average values (for purposes of average max DA values
% graph)
wholeNetwork.DA1avg = [];
wholeNetwork.DA2avg = [];

% Each nucleus is managed within wholeNetwork.nuclei. Each nucleus has a nucleus
% index. globalNuclei is a cell array that is arranged like a connectivity
% matrix, but the elements are connectivity matrices between the nuclei.
% Since each connectivity matrix goes in two directions, wholeNetwork.nuclei is an
% upper triangular array. This means that for any two nuclei indices, the
% lower index must be the first dimension, and the higher index must be the
% second dimension. Along the diagonal, there are self-connection
% connectivity matrices, as well as a parameter list for that nucleus.

nucIndex.PFC = 1;
nucIndex.STRstr = 2;
nucIndex.SNc = 3;
nucIndex.SENa = 4;
nucIndex.SENb = 5;
nucIndex.SENc = 6;
nucIndex.HIPP = 7;
nucIndex.INTh = 8;

% unused areas (indices kept for ease of modularity)
nucIndex.HIPP2 = 9;
nucIndex.INTa = 10;
nucIndex.INTb = 11;
nucIndex.INTc = 12;
nucIndex.STRdi = 13;
nucIndex.STRddi = 14;
nucIndex.GPi = 15;
nucIndex.GPe = 16;
nucIndex.STN = 17;
nucIndex.SNr = 18;
nucIndex.THAL = 19;
nucIndex.SUB = 20;
nucIndex.EC = 21;
nucIndex.DG = 22;
nucIndex.CA3 = 23;
nucIndex.CA1 = 24;

wholeNetwork.numNuclei = 8;
wholeNetwork.nucIndex = nucIndex;

%% reset values for SimNoiseStimulation file

HippStim.BL.b = 0;
HippStim.BL.index = 1;
HippStim.BF.b = 0;
HippStim.BF.index = 2;
HippStim.BRa.b = 0;
HippStim.BRa.index = 3;
HippStim.BRb.b = 0;
HippStim.BRb.index = 4;

HippStim.LF.b = 0;
HippStim.LF.index = 5;
HippStim.LRa.b = 0;
HippStim.LRa.index = 6;
HippStim.LRb.b = 0;
HippStim.LRb.index = 7;

HippStim.FRa.b = 0;
HippStim.FRa.index = 8;
HippStim.FRb.b = 0;
HippStim.FRb.index = 9;

HippStim.RaRb.b = 0;
HippStim.RaRb.index = 10;

HippStim.TB.b = 0;
HippStim.TB.index = 11;
HippStim.TL.b = 0;
HippStim.TL.index = 12;
HippStim.TF.b = 0;
HippStim.TF.index = 13;
HippStim.TRa.b = 0;
HippStim.TRa.index = 14;
HippStim.TRb.b = 0;
HippStim.TRb.index = 15;

HippStim.vars.T = 0;
    
%% All possible combinations of context A, context B, light, bell, and food

myStimType.aBLF = 0;
myStimType.aBF = 1;
myStimType.aBL = 2;
myStimType.aLF = 3;
myStimType.aB = 4;
myStimType.aL = 5;
myStimType.aF = 6;
myStimType.bBLF = 7;
myStimType.bBF = 8;
myStimType.bBL = 9;
myStimType.bLF = 10;
myStimType.bB = 11;
myStimType.bL = 12;
myStimType.bF = 13;
myStimType.aNothing = 14;
myStimType.bNothing = 15;

stimType = myStimType.aNothing;


%%
wholeNetwork.plasticConnections = zeros(wholeNetwork.numNuclei);

%% EDIT: reset the last fire values (after conversation with Ben)
for ii = 1:wholeNetwork.numNuclei
    for jj = 1:wholeNetwork.numNuclei
        wholeNetwork.nuclei{ii}{jj}.lastFire = [];
    end
end
%% EDIT end

clear tmp;

%% %% BUILD NUCLEI

%% Template
% network parameters.
params.Name = 'template';
%Number of excitatory neurons
params.Ne = 1000;
%Number of inhibitory neurons
params.Ni = 0;
params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0.1;
params.percentConnectivity = 0;
% neuron model parameters for the excitatory neurons
params.exA = 0.02;
params.exB = 0.2;
params.exC = -65;
params.exD = 8;
% neuron model parameters for the inhibitory neurons.
params.inA = 0.1;
params.inB = 0.2;
params.inC = -65;
params.inD = 2;
% Number of outgoing connections per neuron.
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.numConnectionsPerNeuron = 100;
% a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
% Izhi&Szat have 20ms max, uniform.
params.delayRange = [1 20];
params.initExWeight = 5;
params.initInWeight = -4;
% weight bounds, applied every millisecond.
params.weightUpperbound = 8;
params.weightLowerbound = -8;
% Izhikevich iters two .5 ms to add to 1ms. This is not really a parameter.
params.timeStep = .5; % 0.5 ms.
params.noiseAmplitude = 6.5;

% WTA local inhibition
params.numConnectionsPerInhibNeuron = 0;
params.numConnectionsPerExciteNeuron = 0;


%% PFC
params.Name = 'PFC';
params.Ne = 1500;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.PFC;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% STR striosomes
params.Name = 'STR';
params.Ne = 150;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.STRstr;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% SNc
params.Name = 'DA';
params.Ne = 100;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.SNc;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% SENa
params.Name = 'SEN-Bell';
params.Ne = 50;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.SENa;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% SENb
params.Name = 'SEN-Food';
params.Ne = 50;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.SENb;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% SENc
params.Name = 'SEN-Light';
params.Ne = 50;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.SENc;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% HIPP
params.Name = 'HIPP';
params.Ne = 20 * (numel(fieldnames(HippStim))-1); %% 20 units per section, 15 diff sections (6 choose 2 is 15)
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 5.5;
tmp = nucIndex.HIPP;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% HIPPINT
params.Name = 'INTh';
params.Ne = 300; %% 20 units per section, 15 diff sections (6 choose 2 is 15)
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.INTh;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% Reset Model DA level
wholeNetwork.DA = 0;

%% %% Create connections

tmpPFC = nucIndex.PFC;
tmpSTRstr = nucIndex.STRstr;

% JJ: PFC-STRstr
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRstr}{tmpSTRstr}.S)*0.1);
ffConnect = [tmpPFC tmpSTRstr];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = wholeNetwork.tmpPFCSTR;%%% FROM 10 JJ EDIT 2-21-2017
weightLowerbound = 0;
weightMultiplicand = 0.12;
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

clear tmpPFC;
clear tmpSTRstr;


tmpSTRstr = nucIndex.STRstr;
tmpSNc = nucIndex.SNc;

% JJ: STRstr-SNc
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNc}{tmpSNc}.S)*1);
ffConnect = [tmpSTRstr tmpSNc];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 0;
weightLowerbound = -10;
weightMultiplicand = -1;
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

clear tmpSNc;
clear tmpSTRstr;

%% Sensory Pathway

tmpSENa = nucIndex.SENa;
tmpSENb = nucIndex.SENb;
tmpSENc = nucIndex.SENc;
tmpINTh = nucIndex.INTh;
tmpSNc = nucIndex.SNc;

% JJ: SENA-INTh
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTh}{tmpINTh}.S)*1);
ffConnect = [tmpSENa tmpINTh];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = wholeNetwork.tmpSENINTMult; %%%% CHANGED ALL THESE FROM 10 TO 5
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = 0;
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

% JJ: SENB-INTh
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTh}{tmpINTh}.S)*1);
ffConnect = [tmpSENb tmpINTh];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = wholeNetwork.tmpSENINTMult;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = 2; % 5 * 2 = 10
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);


% JJ: SENC-INTh
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTh}{tmpINTh}.S)*1);
ffConnect = [tmpSENc tmpINTh];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = wholeNetwork.tmpSENINTMult;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = 0; % 5 * 2 = 10
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);


clear tmpSENa;
clear tmpSENb;
clear tmpINTa;
clear tmpINTb;
clear tmpSNc;


%% Hippocampus

tmpHIPP = nucIndex.HIPP;
tmpINTh = nucIndex.INTh;
tmpSNc = nucIndex.SNc;

% JJ: HIPP-SNc
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNc}{tmpSNc}.S)*1);
ffConnect = [tmpINTh tmpSNc];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = wholeNetwork.tmpMult;%0.06; %%%%%%% set to zero
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

% JJ: HIPP-INTh
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTh}{tmpINTh}.S)*1);
ffConnect = [tmpHIPP tmpINTh];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = wholeNetwork.tmpUpp;%5;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = wholeNetwork.tmpInitHipp; % 5 * 2 = 10
% disp(weightMultiplicand);
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);


clear ffConnect delayRange numConnectionsPerNeuron
clear weightUpperbound weightLowerbound weightMultiplicand

%%

%Init stimuli. It will yield a cell stim with stimNumber of elements, each
%referring to a 1000ms period. Each spike stimuli is a row element of
%stim{stimNumber}, consisting of [neuron timeOfSpike stimmedNuclei]
onset = 400;


% JJ:
tstart = tic;
treset = tic;

% JJ
stimType = 0;
% bParams = [];

% JJ
clear tmpDAWeightChanges tmpDA;
spkPerStep = [];
spkPerSec = [];
timestepPerSecond = [];
s1=RandStream.create('mt19937ar','seed',1);
s2=RandStream.create('mt19937ar','seed',2);
s3=RandStream.create('mt19937ar','seed',3);

previousTime = wholeNetwork.t;
%% Experimenter Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLearningTrials = 100;
numExtinctionTrials = 100;
numTestTrialsAfterExtinction = 100;

trialLength = 3000;

numTrials = wholeNetwork.Num_Learning_Trials ...
    + wholeNetwork.Num_Extinction_Trials ...
    + wholeNetwork.Num_Renewal_Trials% ...
    %+ wholeNetwork.Num_Post_Ext_Trials;

timeBLFStart = wholeNetwork.Num_Learning_Trials*trialLength;
timeBorLStart = timeBLFStart + (wholeNetwork.Num_Extinction_Trials * trialLength);
% timePostExtinctionStart = timeRenewalStart + (wholeNetwork.Num_Renewal_Trials * trialLength);

runTime = numTrials*trialLength;

% % % % 
% % % % wholeNetwork.Num_Extinction_Trials = 100;
% % % % wholeNetwork.Extinction_Trials_in_A = min(1,0 + 1 * floor(rand() * 2));
% % % % 
% % % % % 100 back and forth trials of just B and just L renewal trial of light only (in A or B)
% % % % wholeNetwork.Num_Renewal_Trials = 100;
% % % % wholeNetwork.Renewal_Trials_in_A = min(1,0 + 1 * floor(rand() * 2));

% side note: triallength much be larger than 1000+1750= 2750
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

tmpDA = zeros(runTime,1);

% (fyi: not actually time of first onset) - after that time, then randomly
% choose a time for the next onset as if a new trial is occuring
onset = trialLength+1;

% offset = the difference between when SENA gets input and when PFC gets
% input
offset = 100;

%% Write a file in the folder with info on this run
% fclose(fopen([...
%     'data/',myIdentifier,'/'...
%     'A=',num2str(A),...
%     ', B=',num2str(B),...
%     ', C=',num2str(C),...
%     ', D=',num2str(D),...
%     ', offset=',num2str(offset),'.txt'], 'w'));
myFileIdentifier = [' H-iHmax=',num2str(wholeNetwork.tmpUpp),...
    '_H-iHinit=',num2str(wholeNetwork.tmpInitHipp),...
    '_iH-DA=',num2str(wholeNetwork.tmpMult),...
    '_Sen-INT=',num2str(wholeNetwork.tmpSENINTMult),...
    '_PFC-STR=',num2str(wholeNetwork.tmpPFCSTR),...
    '_BorLbool=',num2str(wholeNetwork.BorLbool),...
    ];
%     '_ExTrA=',num2str(wholeNetwork.Extinction_Trials_in_A),...
%     '_ReTrA=',num2str(wholeNetwork.Renewal_Trials_in_A),...
%     '_PExTrA=',num2str(wholeNetwork.Post_Ext_Trials_in_A)];

myIdentifier = [wholeNetwork.expId, myFileIdentifier];

% % other possible variables:
%     '_NLearnT=',num2str(wholeNetwork.Num_Learning_Trials),...
%     '_NExtT=',num2str(wholeNetwork.Num_Ext_Trials),...
%     '_NRenewalT=',num2str(wholeNetwork.Num_Renewal_Trials),...
%     '_Str-IH-DA=',num2str(wholeNetwork.tmpMult),...
%     '_Init-H-IH=',num2str(wholeNetwork.tmpInitHipp),...
    
mkdir('data', myIdentifier);
disp(myFileIdentifier);

HippStim.vars.T = 1;
currentContext = 1;

% Run the network for the entirety of runTime.
for t = (previousTime+1):(previousTime+runTime)
    wholeNetwork.t = t;%% >> fixed below line
   
    % if there are more than two points, and it's at the end of a trial, or
    % its the end of the run:
    if length(spkPerStep) > 1 && (mod(t-(previousTime),trialLength) == 1750 || t == runTime)
        figure(2);
        subplot(3,1,1);
        plot(spkPerStep);
        title('Spikes Per Timestep (simulated milliseconds)');
        subplot(3,1,2);
        plot(spkPerSec);
        title('Spikes Per Second (real life seconds)');
        subplot(3,1,3);
        plot(timestepPerSecond);
        title('Simulated Milliseconds per Real Life Second  //  Timesteps per Second');
        print(gcf,['data/',myIdentifier,'/Performance.png'],'-dpng');
        close(gcf);
    end
    
%     %% JJ to get new pfc noise every time
    if t > 10 && (mod(t-(previousTime),trialLength) == 0 || t == runTime)
        save(['data/',myIdentifier,'/t',myFileIdentifier,'.mat'],'wholeNetwork')
        save(['data/',myIdentifier,'/t',myFileIdentifier,'.mat'],'tmpDA', '-append')
        save(['data/',myIdentifier,'/t',myFileIdentifier,'.mat'],'onset', '-append')


        %% Cut down on the size of WholeNetwork.firings
        % thus it will only hold info from one trial from now in the
        % past

        intMemory = 1;

        wholeNetwork.firings = ...
            wholeNetwork.firings(...
                wholeNetwork.firings(:,1)...
                    > (max(1,t-(intMemory*trialLength)))...
                ...
            ,:);

        %% Cut down on the size of conductingpotentials
        % thus it will only hold info from one trial from now in the
        % past

        for tmpii = 1:wholeNetwork.numNuclei
            for tmpjj = 1:wholeNetwork.numNuclei
                if isfield(wholeNetwork.nuclei{tmpii}{tmpjj}, 'conductingPotentials')
                if length(wholeNetwork.nuclei{tmpii}{tmpjj}.conductingPotentials) > intMemory*trialLength
                    wholeNetwork.nuclei{tmpii}{tmpjj}.conductingPotentials = ...
                        wholeNetwork.nuclei{tmpii}{tmpjj}.conductingPotentials(...
                            wholeNetwork.nuclei{tmpii}{tmpjj}.conductingPotentials(:,1)...
                                > (max(1,t-(intMemory*trialLength)))...
                            ...
                        ,:);
                end
                end
            end
        end


        %%
        spkPerStep = [spkPerStep; length(wholeNetwork.firings)/(intMemory*trialLength)];
        spkPerSec = [spkPerSec; length(wholeNetwork.firings)/(toc(treset))];
        timestepPerSecond = [timestepPerSecond; trialLength/(toc(treset))];
        treset = tic;
        if length(spkPerStep) == 1
            spkPerStep = spkPerStep * intMemory;
            spkPerSec = spkPerSec * intMemory;
        end
        a = whos('wholeNetwork');
        
        disp(['     -- Spikes/Timestep:  ',num2str(spkPerStep(end))]);
        disp(['     -- Spikes/Second:    ',num2str(spkPerSec(end))]);
        disp(['     -- Timesteps/Second: ',num2str(timestepPerSecond(end))]);
        disp(['     -- Size of wholeNetwork: ', num2str(a.bytes.*9.53674e-7),' MB']);

        clear intMemory

        % new stimulus onset time
        onset = ceil(rand*1000) + t;%%%%%%%%%%%%
        % reset the random streams
        s1=RandStream.create('mt19937ar','seed',1);
        s2=RandStream.create('mt19937ar','seed',2);
        s3=RandStream.create('mt19937ar','seed',3);

       
%% Reset Eligibility Traces (optional)
%         for ii = 1:wholeNetwork.numNuclei
%             for jj = 1:wholeNetwork.numNuclei
%                 if ~isempty(wholeNetwork.nuclei{ii}{jj}.lastFire)
%                     wholeNetwork.nuclei{ii}{jj}.eligibilityTrace = zeros(size(wholeNetwork.nuclei{ii}{jj}.eligibilityTrace));
%                 end
%             end
%         end
%         disp(stimType);



        %% %% Set Stimulus parameters:
        %% Basic Conditioning
        stimType = myStimType.aLF;
        
        %% Extinction:
        %    Half of the time extinguish in other context
        
        if t >= timeBLFStart && t < timeBorLStart
            if wholeNetwork.Extinction_Trials_in_A
                stimType = myStimType.aL;
                %% this repeated code resets food since last reward to 1 when the context changes
                if currentContext == 0
                    HippStim.vars.T = 1;
                end
                currentContext = 1;
            else
                stimType = myStimType.bL;
                if currentContext == 1
                    HippStim.vars.T = 1;
                end
                currentContext = 0;
            end
        elseif t >= timeBorLStart %&& t < timePostExtinctionStart
            if wholeNetwork.Renewal_Trials_in_A
                if wholeNetwork.BorLbool
                    stimType = myStimType.aB;
                else
                    stimType = myStimType.aL;
                end
                if currentContext == 0
                    HippStim.vars.T = 1;
                end
                currentContext = 1;
            else
                if wholeNetwork.BorLbool
                    stimType = myStimType.bB;
                else
                    stimType = myStimType.bL;
                end
                if currentContext == 1
                    HippStim.vars.T = 1;
                end
                currentContext = 0;
            end
%         elseif t >= timePostExtinctionStart
%             if wholeNetwork.Post_Ext_Trials_in_A
%                 stimType = myStimType.aL;
%             else
%                 stimType = myStimType.bL;
%             end
        end

    end
   
%   stimType 0 = Both Bell and Food
%   stimType 1 = Bell only
%   stimType 2 = Food only
%   stimType 3 = Nothing
      
    networkIteratePARALLEL

    %% Time remainging estimation:
    if mod(t,runTime/1000) == 0 && t ~= runTime
        if t < runTime/100 || mod(t,runTime/100) == 0
            fprintf('% 5.1f', (100*t/runTime));
            disp(['% -- ', 'Time Passed (HH:MM:SS):  ', ...
                datestr((toc(tstart))/24/60/60, 'HH:MM:SS')] );
            disp(['       Time Remaining (HH:MM:SS):  ', ...
                datestr((toc(tstart)/(t/runTime)-toc(tstart))/24/60/60, 'HH:MM:SS')] );
            disp(['     Estimated Finish (HH:MM:SS):  ', ...
                datestr((toc(tstart)/(t/runTime)-toc(tstart))/24/60/60+now, 'HH:MM:SSPM')] );
        end
    elseif t == runTime
        disp([num2str(100*t/runTime), ...
            '% -- ', 'Total Time    (HH:MM:SS):  ', ...
            datestr((toc(tstart))/24/60/60, 'HH:MM:SS')] );
    end
    
end

%%
% %% how to reset myMat:
myMat = []
save(['data/data.mat'],'myMat')
%%

load(['data/data.mat'],'myMat')

myMat = [myMat;wholeNetwork.tmpUpp,...
    wholeNetwork.tmpInitHipp,...
    wholeNetwork.tmpMult,...
    wholeNetwork.tmpSENINTMult,...
    wholeNetwork.tmpPFCSTR,0];

save(['data/data.mat'],'myMat', '-append');

clear nucleiNumber previousTime runTime
clear t jj ii runTime inputNoise
clear tstart treset stimType
clear s1 s2 s3
clear all
close all
pngNucleiDriverJeff
