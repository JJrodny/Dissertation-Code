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
%% open
% if matlabpool('size') == 0 % checking to see if my pool is already open
%     matlabpool(6)
% end
%% parallel loops for each iterate loop
% pmode close
% pmode open 6

%% Actual code
% clear all;
close all;
wholeNetwork.expId = datestr(now,30);
wholeNetwork.t = 0;
wholeNetwork.firings = [];


% Each nucleus is managed within wholeNetwork.nuclei. Each nucleus has a nucleus
% index. globalNuclei is a cell array that is arranged like a connectivity
% matrix, but the elements are connectivity matrices between the nuclei.
% Since each connectivity matrix goes in two directions, wholeNetwork.nuclei is an
% upper triangular array. This means that for any two nuclei indices, the
% lower index must be the first dimension, and the higher index must be the
% second dimension. Along the diagonal, there are self-connection
% connectivity matrices, as well as a parameter list for that nucleus.

% %% PRE PFC
% % network parameters.
% params.Name = 'prePFC';
% %Number of excitatory neurons
% params.Ne = 1000;
% %Number of inhibitory neurons
% params.Ni = 0;
% params.N = params.Ne + params.Ni;
% % params.percentConnectivity = 0.1;
% params.percentConnectivity = 0;
% % neuron model parameters for the excitatory neurons
% params.exA = 0.02;
% params.exB = 0.2;
% params.exC = -65;
% params.exD = 8;
% % neuron model parameters for the inhibitory neurons.
% params.inA = 0.1;
% params.inB = 0.2;
% params.inC = -65;
% params.inD = 2;
% % Number of outgoing connections per neuron.
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% % params.numConnectionsPerNeuron = 100;
% % a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
% % Izhi&Szat have 20ms max, uniform.
% params.delayRange = [1 20];
% params.initExWeight = 5;
% params.initInWeight = -4;
% % weight bounds, applied every millisecond.
% params.weightUpperbound = 8;
% params.weightLowerbound = -8;
% % Izhikevich iters two .5 ms to add to 1ms. This is not really a parameter.
% params.timeStep = .5; % 0.5 ms.
%
% wholeNetwork.nuclei{8}{8} = networkBuild(params);
% wholeNetwork.nuclei{8}{8}.nucleiIndex = 8;
% clear params

nucIndex.PFC = 1;
nucIndex.STRstr = 2;
nucIndex.SNc = 3;
nucIndex.SENa = 4;
nucIndex.SENb = 5;
nucIndex.INTa = 6;
nucIndex.INTb = 7;

% nucIndex.SUB = 8;
% nucIndex.EC = 9;
% nucIndex.DG = 10;
% nucIndex.CA3 = 11;
% nucIndex.CA1 = 12;

nucIndex.HIPP = 8;
nucIndex.HIPP2 = 9;
nucIndex.STRdi = 14;
nucIndex.STRddi = 15;
nucIndex.GPi = 16;
nucIndex.GPe = 17;
nucIndex.STN = 18;
nucIndex.SNr = 19;
nucIndex.THAL = 20;

wholeNetwork.numNuclei = 7;
wholeNetwork.nucIndex = nucIndex;

wholeNetwork.plasticConnections = zeros(wholeNetwork.numNuclei);

%% JJ EDIT
for ii = 1:wholeNetwork.numNuclei
    for jj = 1:wholeNetwork.numNuclei
        wholeNetwork.nuclei{ii}{jj}.lastFire = [];
    end
end
%% JJ end
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0

clear tmp;

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
params.Ne = 1000;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.PFC;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% STR striosomes
% params.Name = 'STRstr';
params.Name = 'STR';
params.Ne = 100;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.STRstr;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% SNc
% params.Name = 'SNc';
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
% params.Name = 'SENa';
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
% params.Name = 'SENb';
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

%% INTa
% params.Name = 'INTa';
params.Name = 'INT-Bell';
params.Ne = 50;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.INTa;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% INTb
% params.Name = 'INTb';
params.Name = 'INT-Food';
params.Ne = 50;
params.Ni = 0;
params.N = params.Ne + params.Ni;
params.percentConnectivity = 0;
params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
params.noiseAmplitude = 6.5;
tmp = nucIndex.INTb;
wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

% %% SUB
% params.Name = 'SUB';
% params.Ne = 100;
% params.Ni = 0;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.noiseAmplitude = 6.5;
% tmp = nucIndex.SUB;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% 
% %% EC
% params.Name = 'EC';
% params.Ne = 100;
% params.Ni = 0;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.noiseAmplitude = 6.5;
% tmp = nucIndex.EC;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% 
% %% Hipp
% params.Name = 'Hipp';
% params.Ne = 400;
% params.Ni = 100;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0.0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.numConnectionsPerInhibNeuron = 400;
% params.numConnectionsPerExciteNeuron = 100;
% params.noiseAmplitude = 6.5;%% JJ Fix from 11 | 1/27
% params.initExWeight = 5;
% params.initInWeight = -600;
% params.weightUpperbound = 8;
% params.weightLowerbound = -608;
% tmp = nucIndex.HIPP;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% % reset values
% params.numConnectionsPerInhibNeuron = 0;
% params.numConnectionsPerExciteNeuron = 0;
% params.initExWeight = 5;
% params.initInWeight = -4;
% params.weightUpperbound = 8;
% params.weightLowerbound = -8;
% 
% %% Hipp2
% params.Name = 'Hipp2';
% params.Ne = 400;
% params.Ni = 0;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0.0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.numConnectionsPerInhibNeuron = 0;
% params.numConnectionsPerExciteNeuron = 0;
% tmp = nucIndex.HIPP2;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% % reset values
% params.numConnectionsPerInhibNeuron = 0;
% params.numConnectionsPerExciteNeuron = 0;
% params.initExWeight = 5;
% params.initInWeight = -4;
% params.weightUpperbound = 8;
% params.weightLowerbound = -8;
% 
% %% DG
% params.Name = 'DG';
% params.Ne = 400;
% params.Ni = 100;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.numConnectionsPerInhibNeuron = 400;
% params.numConnectionsPerExciteNeuron = 100;
% params.noiseAmplitude = 6.5;%% JJ Fix from 11 | 1/27
% params.initInWeight = -10;
% tmp = nucIndex.DG;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% % reset values
% params.initInWeight = -4;
% params.numConnectionsPerInhibNeuron = 0;
% params.numConnectionsPerExciteNeuron = 0;
% 
% %% CA3
% params.Name = 'CA3';
% params.Ne = 75;
% params.Ni = 25;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.numConnectionsPerInhibNeuron = 75;
% params.numConnectionsPerExciteNeuron = 25;
% params.noiseAmplitude = 6.5;
% tmp = nucIndex.CA3;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;
% % reset values
% params.numConnectionsPerInhibNeuron = 0;
% params.numConnectionsPerExciteNeuron = 0;
% 
% %% CA1
% params.Name = 'CA1';
% params.Ne = 100;
% params.Ni = 0;
% params.N = params.Ne + params.Ni;
% params.percentConnectivity = 0;
% params.numConnectionsPerNeuron = floor(params.N * params.percentConnectivity);
% params.noiseAmplitude = 6.5;
% tmp = nucIndex.CA1;
% wholeNetwork.nuclei{tmp}{tmp} = networkBuild(params);
% wholeNetwork.nuclei{tmp}{tmp}.nucleiIndex = tmp;

%% JJ Comment
% for alpha, used for modifying STDP and STR excitability
wholeNetwork.DA = 0;
%% End JJ

% %% JJ Comment addition
% % init a lastFire = [] for each empty region.
% for n1 = 1:wholeNetwork.numNuclei
%     for n2 = 1:wholeNetwork.numNuclei
%         wholeNetwork.nuclei{n1}{n2}.lastFire = [];
%     end
% end
% clear n1 n2
% %% end JJ

% a vector of inclusive lower bound and inclusive upper bound of the range of conductance delays. in ms
% delayRange = [1 10];
% weightUpperbound = 8;
% weightLowerbound = 0;

% Feed forward connect nuclei index to another nuclei index.
% ffConnect = [1 2];

% JJ: prePFC-PFC
% numConnectionsPerNeuron = length(wholeNetwork.nuclei{1}{1}.S)*.2;
% ffConnect = [8 1];
% weightMultiplicand = .8;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, [46,48], ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand);

%Rand parameters:
% .1-.9
% A = ceil(rand * 3)*.3 - .2
% A = ceil(rand * 3)*.1+.1;
% A = 0.3
% % 10, 100
% B = 10 + floor(rand * 2) * 90;
% B = 10
% % 4,6,8,10
% C = ceil(rand * 3)*1+6;
% C = 10
% % -0.2, -0.7, -1.2
% D = -0.8 + -0.2 * ceil(rand * 3)
% D = -1;

%% JJ: PFC-STR

% tmpPFC = nucIndex.PFC;
% tmpSTRdi = nucIndex.STRdi;
% tmpSTN = nucIndex.STN;
% tmpSTRddi = nucIndex.STRddi;
% tmpSTRstr = nucIndex.STRstr;
% tmpTHAL = nucIndex.THAL;

% % JJ: PFC-STRdi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRdi}{tmpSTRdi}.S)*0.2);
% ffConnect = [tmpPFC tmpSTRdi];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: PFC-STRddi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRddi}{tmpSTRddi}.S)*0.2);
% ffConnect = [tmpPFC tmpSTRddi];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

tmpPFC = nucIndex.PFC;
tmpSTRstr = nucIndex.STRstr;

% JJ: PFC-STRstr
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRstr}{tmpSTRstr}.S)*0.3);
ffConnect = [tmpPFC tmpSTRstr];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
weightMultiplicand = 0.12;
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);


% tmpPFC = nucIndex.PFC;
% tmpEC = nucIndex.EC;
% 
% % JJ: PFC-EC
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpEC}{tmpEC}.S)*0.2);
% ffConnect = [tmpPFC tmpEC];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);


% % JJ: PFC-STN
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTN}{tmpSTN}.S)*0.2);
% ffConnect = [tmpPFC tmpSTN];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: PFC-THAL
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpTHAL}{tmpTHAL}.S)*0.3);
% ffConnect = [tmpPFC tmpTHAL];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

clear tmpPFC;
% clear tmpSTRdi;
% clear tmpSTRddi;
clear tmpSTRstr;
% clear tmpTHAL;
%% JJ: HIPP-STR

% tmpHIPP = nucIndex.HIPP;
% tmpSTRdi = nucIndex.STRdi;
% tmpSTRddi = nucIndex.STRddi;

% % JJ: HIPP-STRdi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRdi}{tmpSTRdi}.S)*0.1);
% ffConnect = [tmpHIPP tmpSTRdi];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% % Chorley Seth: Weights initialized to Min and Max respectively
% weightMultiplicand = 0;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: HIPP-STRddi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRddi}{tmpSTRddi}.S)*0.1);
% ffConnect = [tmpHIPP tmpSTRddi];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% % Chorley Seth: Weights initialized to Min and Max respectively
% weightMultiplicand = 0;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpHIPP;
% clear tmpSTRdi;
% clear tmpSTRddi;

%% JJ: STRdi

% tmpSTRdi = nucIndex.STRdi;
% tmpGPi = nucIndex.GPi;
% tmpSNr = nucIndex.SNr;

% % JJ: STRdi-GPe
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpGPi}{tmpGPi}.S)*1);
% ffConnect = [tmpSTRdi tmpGPi];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: STRdi-SNr
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNr}{tmpSNr}.S)*1);
% ffConnect = [tmpSTRdi tmpSNr];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPi;
% clear tmpSTRdi;
% clear tmpSNr;

%% JJ: STRddi

% tmpSTRddi = nucIndex.STRddi;
% tmpGPe = nucIndex.GPe;

% % JJ: STRddi-GPe
% numConnectionsPerNeuron = length(wholeNetwork.nuclei{tmpGPe}{tmpGPe}.S)*1;
% ffConnect = [tmpSTRddi tmpGPe];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPe;
% clear tmpSTRddi;

%% JJ: STRstr

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

%% JJ: GPi

% tmpTHAL = nucIndex.THAL;
% tmpGPi = nucIndex.GPi;

% % JJ: GPi-THAL
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpTHAL}{tmpTHAL}.S)*1);
% ffConnect = [tmpGPi tmpTHAL];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPi;
% clear tmpTHAL;

%% JJ: GPe

% tmpSTN = nucIndex.STN;
% tmpGPe = nucIndex.GPe;
% tmpGPi = nucIndex.GPi;

% % JJ: GPe-STN
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTN}{tmpSTN}.S)*1);
% ffConnect = [tmpGPe tmpSTN];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: GPe-GPi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpGPi}{tmpGPi}.S)*1);
% ffConnect = [tmpGPe tmpGPi];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPe;
% clear tmpGPi;
% clear tmpSTN;

%% JJ: SNr, SNc

% tmpSNc = nucIndex.SNc;
% tmpSNr = nucIndex.SNr;
% tmpTHAL = nucIndex.THAL;

% % JJ: SNr-SNc
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNc}{tmpSNc}.S)*1);
% ffConnect = [tmpSNr tmpSNc];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: SNr-THAL
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpTHAL}{tmpTHAL}.S)*1);
% ffConnect = [tmpSNr tmpTHAL];
% delayRange = [1 10];
% weightUpperbound = 0;
% weightLowerbound = -10;
% weightMultiplicand = -2;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPe;
% clear tmpSTN;

%% JJ: STN

% tmpSTN = nucIndex.STN;
% tmpGPi = nucIndex.GPi;
% tmpGPe = nucIndex.GPe;
% tmpSNr = nucIndex.SNr;

% % JJ: STN-GPi
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpGPi}{tmpGPi}.S)*1);
% ffConnect = [tmpSTN tmpGPi];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: STN-GPe
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpGPe}{tmpGPe}.S)*1);
% ffConnect = [tmpSTN tmpGPe];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
%
% % JJ: STN-SNr
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNr}{tmpSNr}.S)*1);
% ffConnect = [tmpSTN tmpSNr];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);

% clear tmpGPi;
% clear tmpGPe;
% clear tmpSTN;
% clear tmpSNr;

% %% Hippocampus
% 
% %% Hipp-STR
% 
% tmpHIPP = nucIndex.HIPP;
% tmpHIPP2 = nucIndex.HIPP2;
% tmpSTRstr = nucIndex.STRstr;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSTRstr}{tmpSTRstr}.S)*1);
% ffConnect = [tmpHIPP tmpSTRstr];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.08;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% % SenA and SenB to Hipp
% 
% tmpSENa = nucIndex.SENa;
% tmpSENb = nucIndex.SENb;
% tmpHIPP = nucIndex.HIPP;
% tmpPFC = nucIndex.PFC;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpHIPP}{tmpHIPP}.S)*.5);
% ffConnect = [tmpPFC tmpHIPP];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.08;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpHIPP2}{tmpHIPP2}.S)*.5);
% ffConnect = [tmpPFC tmpHIPP2];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.08;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);


% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpHIPP}{tmpHIPP}.S)*1);
% ffConnect = [tmpSENa tmpHIPP];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.08;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpHIPP}{tmpHIPP}.S)*1);
% ffConnect = [tmpSENb tmpHIPP];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.08;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);


% %% SNc
% 
% tmpSNc = nucIndex.SNc;
% tmpSUB = nucIndex.SUB;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSUB}{tmpSUB}.S)*1);
% ffConnect = [tmpSNc tmpSUB];
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.06;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% %% SUB
% 
% tmpSUB = nucIndex.SUB;
% tmpEC = nucIndex.EC;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpEC}{tmpEC}.S)*1);
% ffConnect = [tmpSUB tmpEC];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpSUB = nucIndex.SUB;
% tmpCA1 = nucIndex.CA1;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpCA1}{tmpCA1}.S)*1);
% ffConnect = [tmpSUB tmpCA1];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% %% EC
% 
% tmpEC = nucIndex.EC;
% tmpSUB = nucIndex.SUB;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSUB}{tmpSUB}.S)*1);
% ffConnect = [tmpEC tmpSUB];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpEC = nucIndex.EC;
% tmpCA1 = nucIndex.CA1;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpCA1}{tmpCA1}.S)*1);
% ffConnect = [tmpEC tmpCA1];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpEC = nucIndex.EC;
% tmpCA3 = nucIndex.CA3;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpCA3}{tmpCA3}.S)*1);
% ffConnect = [tmpEC tmpCA3];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpEC = nucIndex.EC;
% tmpDG = nucIndex.DG;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpDG}{tmpDG}.S)*1);
% ffConnect = [tmpEC tmpDG];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpEC = nucIndex.EC;
% tmpPFC = nucIndex.PFC;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpPFC}{tmpPFC}.S)*1);
% ffConnect = [tmpEC tmpPFC];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% %% CA3
% 
% tmpCA3 = nucIndex.CA3;
% tmpCA1 = nucIndex.CA1;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpCA1}{tmpCA1}.S)*1);
% ffConnect = [tmpCA3 tmpCA1];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% %% CA1
% 
% tmpCA1 = nucIndex.CA1;
% tmpEC = nucIndex.EC;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpEC}{tmpEC}.S)*1);
% ffConnect = [tmpCA1 tmpEC];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% 
% tmpCA1 = nucIndex.CA1;
% tmpSUB = nucIndex.SUB;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSUB}{tmpSUB}.S)*1);
% ffConnect = [tmpCA1 tmpSUB];
% % wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% %% DG
% 
% tmpDG = nucIndex.DG;
% tmpCA3 = nucIndex.CA3;
% 
% numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpCA3}{tmpCA3}.S)*.2);
% ffConnect = [tmpDG tmpCA3];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
% delayRange = [1 10];
% weightUpperbound = 10;
% weightLowerbound = 0;
% weightMultiplicand = 0.12;
% wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
%     numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
%     weightMultiplicand, 0);
% 
% clear tmpSNc;
% clear tmpDG;
% clear tmpCA1;
% clear tmpCA3;
% clear tmpPFC;
% clear tmpSUB;

%% Sensory Pathway

tmpSENa = nucIndex.SENa;
tmpSENb = nucIndex.SENb;
tmpINTa = nucIndex.INTa;
tmpINTb = nucIndex.INTb;
tmpSNc = nucIndex.SNc;

% JJ: SENA-INTA
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTa}{tmpINTa}.S)*1);
ffConnect = [tmpSENa tmpINTa];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = 0; %%%%%%% set to zero
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

% JJ: SENB-INTB
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpINTb}{tmpINTb}.S)*1);
ffConnect = [tmpSENb tmpINTb];
wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
% Chorley Seth: Weights initialized to Min and Max respectively
weightMultiplicand = 2; % 5 * 2 = 10
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

% JJ: INTA-DA
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNc}{tmpSNc}.S)*1);
ffConnect = [tmpINTa tmpSNc];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
weightMultiplicand = 0.12;%0.12 % 5 * 0.12 = 0.6
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

% JJ: INTB-DA
numConnectionsPerNeuron = ceil(length(wholeNetwork.nuclei{tmpSNc}{tmpSNc}.S)*1);
ffConnect = [tmpINTb tmpSNc];
% wholeNetwork.plasticConnections(ffConnect(1),ffConnect(2)) = 1;
delayRange = [1 10];
weightUpperbound = 10;
weightLowerbound = 0;
weightMultiplicand = 0.12;%0.12
wholeNetwork = buildConnections(wholeNetwork, ffConnect, delayRange, ...
    numConnectionsPerNeuron, weightLowerbound, weightUpperbound, ...
    weightMultiplicand, 0);

clear tmpSENa;
clear tmpSENb;
clear tmpINTa;
clear tmpINTb;
clear tmpSNc;


clear ffConnect delayRange numConnectionsPerNeuron
clear weightUpperbound weightLowerbound weightMultiplicand

%Init stimuli. It will yield a cell stim with stimNumber of elements, each
%referring to a 1000ms period. Each spike stimuli is a row element of
%stim{stimNumber}, consisting of [neuron timeOfSpike stimmedNuclei]
onset = 400;

%% JJ test same numbers generated
% aaa = [];
% bbb = [];
% %(then plot(aaa(1:1000),aaa(1001:2000)) to test same numbers generated)
%% JJ
% JJ Comment out:
% pngSimStimuliInitJeff2;
% a vector pattern of stimulation, where 0 implies no stimulus.
%% stimPattern = [1 0 0 2 0 0 1 0 0 2 0 0];
% stimPattern = [1 1 2 2 3 3 4 4 5 5];
% stimPattern = [1,2,1,2,2,2,1,2,2,1,1,2,2,1,2,2,1];% 2 3 1 2 3 1 2 3];
% stimPattern = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];% 2 3 1 2 3 1 2 3];
%
% % Determine background noise. Must be between 0 and 1.
% % This corresponds to the likelihood that a neuron will fire due to noise
% % per second.
% bgNoise = 100;
% % The amplitude of the noise, in mV.
% noiseAmplitude = 6.5;

% firings = [];

% JJ:
tstart = tic;
treset = tic;

% JJ
stimType = 0;
% bParams = [];

% JJ
clear tmpDAWeightChanges tmpDA;
% % tmpDAWeightChanges = [0,0,0];
% tmpDA = zeros(size(runTime),1); % []
spkPerStep = [];
spkPerSec = [];
timestepPerSecond = [];
s1=RandStream.create('mt19937ar','seed',1);
s2=RandStream.create('mt19937ar','seed',2);

previousTime = wholeNetwork.t;
%% Experimenter Input
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numLearningTrials = 100;
numExtinctionTrials = 500;
numTestTrialsAfterExtinction = 100;

trialLength = 5000; %10000

numTrials = numLearningTrials + numExtinctionTrials + numTestTrialsAfterExtinction; %100 % #1 to #(numTrials-numNonLearningTrials=Learning Trials)
numNonLearningTrials = numExtinctionTrials + numTestTrialsAfterExtinction; % #(numTrials-numNonLearningTrials=Learning Trials) to #(200-11)
numTestAfterExtTrials = numTestTrialsAfterExtinction; % #(200-11) to #(200-1)
runTime = numTrials*trialLength;%2000*120;%length(stimPattern)*2000;%900


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

myIdentifier = datestr(now,'mm-dd HH.MM.SS.FFF')
mkdir('data',myIdentifier);

% Write a file in the folder with info on this run
% fclose(fopen([...
%     'data/',myIdentifier,'/'...
%     'A=',num2str(A),...
%     ', B=',num2str(B),...
%     ', C=',num2str(C),...
%     ', D=',num2str(D),...
%     ', offset=',num2str(offset),'.txt'], 'w'));

% wholeNetwork.firings = zeros(runTime*2,3);
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
       
    end
           
%     %% JJ to get new pfc noise every time
    if t > 10 && (mod(t-(previousTime),trialLength) == 0 || t == runTime)
%         t
        save(['data/',myIdentifier,'/t',int2str(t),'.mat'],'wholeNetwork')
        save(['data/',myIdentifier,'/t',int2str(t),'.mat'],'tmpDA', '-append')
        save(['data/',myIdentifier,'/t',int2str(t),'.mat'],'onset', '-append')


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
% %             b = whos('tmpDA');
%         c = wholeNetwork.nuclei;
%         d = whos('c');
%         e = wholeNetwork.nuclei{1}{1};
%         f = whos('e');
%         g = wholeNetwork.nuclei{1}{2};
%         h = whos('g');
%         i = wholeNetwork.firings;
%         j = whos('i');
%         kk = wholeNetwork.nuclei{2}{3};
%         k = whos('kk');
%         ll = wholeNetwork.nuclei{2}{2};
%         l = whos('ll');
%         mm = wholeNetwork.nuclei{4}{6};
%         m = whos('mm');
%         oo = wholeNetwork.nuclei{4}{4};
%         o = whos('oo');
%         pp = wholeNetwork.nuclei{5}{7};
%         p = whos('pp');
%         qq = wholeNetwork.nuclei{5}{5};
%         q = whos('qq');
%         rr = wholeNetwork.nuclei{6}{3};
%         r = whos('rr');
%         ss = wholeNetwork.nuclei{6}{6};
%         s = whos('ss');
%         tt = wholeNetwork.nuclei{7}{3};
%         ttt = whos('tt');
%         uu = wholeNetwork.nuclei{7}{7};
%         u = whos('uu');
        disp(['     -- Spikes/Timestep:  ',num2str(spkPerStep(end))]);
        disp(['     -- Spikes/Second:    ',num2str(spkPerSec(end))]);
        disp(['     -- Timesteps/Second: ',num2str(timestepPerSecond(end))]);
        disp(['     -- Size of wholeNetwork: ', num2str(a.bytes.*9.53674e-7),' MB']);
%             disp(['     -- Size of tmpDA:        ', num2str(b.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei:       ', num2str(d.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{1}{1}: ', num2str(f.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{1}{2}: ', num2str(h.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.firings:      ', num2str(j.bytes.*9.53674e-7),'MB']);
% 
%         disp(['     -- Size of wholeNetwork.nuclei{2}{3}: ', num2str(k.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{2}{2}: ', num2str(l.bytes.*9.53674e-7),'MB']);
% 
%         disp(['     -- Size of wholeNetwork.nuclei{4}{6}: ', num2str(m.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{4}{4}: ', num2str(o.bytes.*9.53674e-7),'MB']);
% 
%         disp(['     -- Size of wholeNetwork.nuclei{5}{7}: ', num2str(p.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{5}{5}: ', num2str(q.bytes.*9.53674e-7),'MB']);
% 
%         disp(['     -- Size of wholeNetwork.nuclei{6}{3}: ', num2str(r.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{6}{6}: ', num2str(s.bytes.*9.53674e-7),'MB']);
% 
%         disp(['     -- Size of wholeNetwork.nuclei{7}{3}: ', num2str(ttt.bytes.*9.53674e-7),'MB']);
%         disp(['     -- Size of wholeNetwork.nuclei{7}{7}: ', num2str(u.bytes.*9.53674e-7),'MB']);

        %             save(['test',int2str(t),'.mat'],'firings','-appen
        %             d')


        clear intMemory
%         clear a b c d e f g h i j;
%         clear k kk l ll m mm o oo p pp q qq r rr s ss tt ttt u uu;

        % new stimulus onset time
        onset = ceil(rand*1000) + t;%%%%%%%%%%%%
        % reset the random streams
        s1=RandStream.create('mt19937ar','seed',1);
        s2=RandStream.create('mt19937ar','seed',2);

       
%% Reset Eligibility Traces
%         for ii = 1:wholeNetwork.numNuclei
%             for jj = 1:wholeNetwork.numNuclei
%                 if ~isempty(wholeNetwork.nuclei{ii}{jj}.lastFire)
%                     wholeNetwork.nuclei{ii}{jj}.eligibilityTrace = zeros(size(wholeNetwork.nuclei{ii}{jj}.eligibilityTrace));
%                 end
%             end
%         end
    end
   
%   stimType 0 = Both Bell and Food
%   stimType 1 = Bell only
%   stimType 2 = Food only
%   stimType 3 = Nothing
   
%     % 2 trials of nothing (not sure why 4 gives 2)
%     % will need to fix when s1 and s2 and onset are initialized if we dont
%     % have these trials
%     if t < (4*trialLength)
%         stimType = 3;
%     else
        stimType = 0;
%     end

%     numLearningTrials = 3;
%     numExtinctionTrials = 3;
%     numTestTrialsAfterExtinction = 6;
%     numTrials = numLearningTrials + numExtinctionTrials + numTestTrialsAfterExtinction; %100 % #1 to #(numTrials-numNonLearningTrials=Learning Trials)
%     numNonLearningTrials = numExtinctionTrials + numTestTrialsAfterExtinction; % #(numTrials-numNonLearningTrials=Learning Trials) to #(200-11)
%     numTestAfterExtTrials = numTestTrialsAfterExtinction; % #(200-11) to #(200-1)

    % Bell only
    if t > trialLength * (numLearningTrials) % 100-201
        stimType = 1;
    end
    % Food only
    if t > trialLength * (numLearningTrials + numExtinctionTrials) % 201-300
        stimType = 2;
    end
    % Final Bell-Food trials
    if t > trialLength * (numLearningTrials + numExtinctionTrials + numTestTrialsAfterExtinction/2) % 301-400
        stimType = 0;
    end

%     % two final Bell-Food
%     if t > runTime-((round(numTestAfterExtTrials/2))*trialLength)
%         stimType = 0;
%     % Bell only
%     elseif t > runTime-((numNonLearningTrials+numTestAfterExtTrials)*trialLength)
%         stimType = 1;
%     % Food only
%     elseif t > runTime-(numTestAfterExtTrials*trialLength)
%         stimType = 2;
%     end
   
    networkIterate
%     networkIterateParallel

    %% JJ:
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
    % End JJ Edit
end
% ddd
% pmode close;
% clear willFire neuronFireTimes stimulusPattern
clear nucleiNumber previousTime runTime
clear t jj ii runTime inputNoise
clear tstart treset stimType
clear s1 s2
pngNucleiDriverJeff