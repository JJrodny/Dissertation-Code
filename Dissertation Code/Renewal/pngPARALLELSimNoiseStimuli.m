function output = pngPARALLELSimNoiseStimuli(stimType, s1, s2, s3, wholeNetwork, t, onset, offset, nucleiNumber, HippStim, myStimType)

% (taken from Ben St. Clair's code):
% Determine background noise. Must be between 0 and 1.
% This corresponds to the likelihood that a neuron will fire due to noise
% per second.
% The amplitude of the noise, in mV.
% noiseAmplitude = 6.5;

% onset = the time the trial officiall starts - the onset of the trial
% offset = the time inbetween when the trial starts and the first stimuli
% is presented - the stimuli offset

%% Visualized times of stimulus onsets and changes in stimuli:
% 0-500-1000
%   500-1000-1500
%       1000-1500-2000
% 
% 0-750-1000
%   750-1000-1500-1750
%            1500-1750-2500
%
% 0-250-500-750-1000
%           750-1000-1250-1500-1750
%                         1500-1750-2000-2250-2500
%
% 0-500-1000
% 0-500-1000
%   500-1000-1500
%            

%% Combinations of Bell, Light, Food, etc for reference
% myStimType.aBLF = 0;
% myStimType.aBF = 1;
% myStimType.aBL = 2;
% myStimType.aLF = 3;
% myStimType.aB = 4;
% myStimType.aL = 5;
% myStimType.aF = 6;
% myStimType.bBLF = 7;
% myStimType.bBF = 8;
% myStimType.bBL = 9;
% myStimType.bLF = 10;
% myStimType.bB = 11;
% myStimType.bL = 12;
% myStimType.bF = 13;
% myStimType.Nothing = 14;

noiseAmplitude = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.noiseAmplitude;

numNeurons = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.N;
inputNoise = zeros(numNeurons,1);

switch nucleiNumber
    case wholeNetwork.nucIndex.PFC
        % PFC
        mynoiseAmplitude = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.noiseAmplitude;
        
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        iN1 = inputNoise(1:floor(numNeurons/3));
%         iN2 = inputNoise(floor(numNeurons/3):2*floor(numNeurons/3));
        iN3 = inputNoise(floor(2*numNeurons/3+1):end);
        iN2to1 = inputNoise(1:2*floor(numNeurons/3));
        iN3to2 = inputNoise(floor(numNeurons/3)+1:end);

        % s1 = tmp1 = Bell
        % s2 = tmp2 = light
        % s3 = tmp3 = food
        
        % B ONLY, during t=0-500
        if t >= onset + offset ...
                && t < onset + 500 + offset ...
                && (stimType == myStimType.aBF || ...
                    stimType == myStimType.aB || ...
                    stimType == myStimType.bBF || ...
                    stimType == myStimType.bB)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN3to2;tmp1];
        % L ONLY, during t=0-500 
        elseif t >= onset + offset ...
                && t < onset + 500 + offset ...
                && (stimType == myStimType.aLF || ...
                    stimType == myStimType.aL || ...
                    stimType == myStimType.bLF || ...
                    stimType == myStimType.bL)
                
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN1;tmp2;iN3];
        % B AND L during t=0-500
        elseif t >= onset+offset ...
                && t < onset + offset + 500 ...
                && (stimType == myStimType.aBLF || ...
                    stimType == myStimType.aBL || ...
                    stimType == myStimType.bBLF || ...
                    stimType == myStimType.bBL)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN1;tmp2;tmp1];
        
       
        % F ONLY during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aF || ...
                    stimType == myStimType.bF)
                
            tmp3 = rand(s3,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [tmp3;iN2to1];
            
        % L ONLY during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aL || ...
                    stimType == myStimType.bL)
                
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN1;tmp2;iN3];
            
        % B ONLY during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aB || ...
                    stimType == myStimType.bB)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN3to2;tmp1];
            
        % L and F during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aLF || ...
                    stimType == myStimType.bLF)
                
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp3 = rand(s3,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [tmp3;tmp2;iN3];
            
        % B and F during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aBF || ...
                    stimType == myStimType.bBF)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp3 = rand(s3,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [tmp3;iN2;tmp1];
            
        % B and L during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aBL || ...
                    stimType == myStimType.bBL)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [iN1;tmp2;tmp1];
            
        % B L and F during t=500-1000
        elseif t >= onset+offset + 500 ...
                && t < onset + offset + 1000 ...
                && (stimType == myStimType.aBLF || ...
                    stimType == myStimType.bBLF)
                
            tmp1 = rand(s1,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp2 = rand(s2,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            tmp3 = rand(s3,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [tmp3;tmp2;tmp1];
            
            
        % F during t=1000-1500
        elseif t >= onset+offset + 1000 ...
                && t < onset + offset + 1500 ...
                && (stimType == myStimType.aBLF || ...
                    stimType == myStimType.aLF || ...
                    stimType == myStimType.aBF || ...
                    stimType == myStimType.aF || ...
                    stimType == myStimType.bBLF || ...
                    stimType == myStimType.bLF || ...
                    stimType == myStimType.bBF || ...
                    stimType == myStimType.bF)
                
            tmp3 = rand(s3,floor(numNeurons/3),1)*13-(13-noiseAmplitude);
            inputNoise = [tmp3;iN2to1];
            
        else
            % neither
        end
        clear numNeurons tmp1 tmp2 tmp3;
        
    case wholeNetwork.nucIndex.SENa
        % SEN A
        noiseAmplitude = 6;% from 6
        if t >= onset && t <= onset + 15 ...
                && (stimType == myStimType.aB || ...
                    stimType == myStimType.aBL || ...
                    stimType == myStimType.aBLF || ...
                    stimType == myStimType.aBF || ...
                    stimType == myStimType.bB || ...
                    stimType == myStimType.bBL || ...
                    stimType == myStimType.bBLF || ...
                    stimType == myStimType.bBF)
            noiseAmplitude = noiseAmplitude + 2.5;% * (t - onset + 1);
        end
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        
    case wholeNetwork.nucIndex.SENc
        % SEN c
        noiseAmplitude = 6;% from 6
        if t >= onset && t <= onset + 15 ...
                && (stimType == myStimType.aL || ...
                    stimType == myStimType.aBL || ...
                    stimType == myStimType.aBLF || ...
                    stimType == myStimType.aLF || ...
                    stimType == myStimType.bL || ...
                    stimType == myStimType.bBL || ...
                    stimType == myStimType.bBLF || ...
                    stimType == myStimType.bLF)
            noiseAmplitude = noiseAmplitude + 2.5;% * (t - (onset + 500) + 1);
        end
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        
    case wholeNetwork.nucIndex.SENb
        % SEN B
        noiseAmplitude = 6;% from 6
        if t >= onset + 500 && t <= onset + 515 ...
                && (stimType == myStimType.aF || ...
                    stimType == myStimType.aBF || ...
                    stimType == myStimType.aBLF || ...
                    stimType == myStimType.aLF || ...
                    stimType == myStimType.bF || ...
                    stimType == myStimType.bBF || ...
                    stimType == myStimType.bBLF || ...
                    stimType == myStimType.bLF)
            noiseAmplitude = noiseAmplitude + 2.5;% * (t - (onset + 500) + 1);
        end
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        
    case wholeNetwork.nucIndex.STRstr
        % STR
%         noiseAmplitude = 6;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.INTa
        % INT A
        noiseAmplitude = 2;% from 6
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.INTb
        % INT B
        noiseAmplitude = 2;% from 6
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.INTc
        % INT B
        noiseAmplitude = 2;% from 6
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.SNc
        % DA
%         noiseAmplitude = 4;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.HIPP
        % HIPP
%         noiseAmplitude = 6;

    %% refernce:
        % Bell
        % Light
        % Food
        % RoomA
        % RoomB
        % TimeSinceReward

        %BL,BF,BRa,BRb,BT
        %LF,LRa,LRb,LT
        %FRa,FRb,FT
        %RaRb,RaT
        %RbT


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

        HippStim.vars.T = HippStim.vars.T * 0.99999; % all the time

        
        % reference:
        % myStimType.aBLF = 0;
        % myStimType.aBF = 1;
        % myStimType.aBL = 2;
        % myStimType.aLF = 3;
        % myStimType.aB = 4;
        % myStimType.aL = 5;
        % myStimType.aF = 6;
        % myStimType.bBLF = 7;
        % myStimType.bBF = 8;
        % myStimType.bBL = 9;
        % myStimType.bLF = 10;
        % myStimType.bB = 11;
        % myStimType.bL = 12;
        % myStimType.bF = 13;
        % myStimType.Nothing = 14;
        
        if (stimType == myStimType.aB || ...
               stimType == myStimType.aL || ...
               stimType == myStimType.aF || ...
               stimType == myStimType.aBL || ...
               stimType == myStimType.aBF || ...
               stimType == myStimType.aLF || ...
               stimType == myStimType.aBLF)
            HippStim.TRa.b = HippStim.vars.T;
        elseif (stimType == myStimType.bB || ...
               stimType == myStimType.bL || ...
               stimType == myStimType.bF || ...
               stimType == myStimType.bBL || ...
               stimType == myStimType.bBF || ...
               stimType == myStimType.bLF || ...
               stimType == myStimType.bBLF)
            HippStim.TRb.b = HippStim.vars.T;
        end
            
        
        if t >= onset && t <= onset + 15
        
            if (stimType == myStimType.aB || ...
                    stimType == myStimType.aBF)
                HippStim.TB.b = HippStim.vars.T;
                HippStim.BRa.b = 1;
                
            elseif (stimType == myStimType.bB || ...
                    stimType == myStimType.bBF)
                HippStim.TB.b = HippStim.vars.T;
                HippStim.BRb.b = 1;
            
                
            elseif (stimType == myStimType.aL || ...
                    stimType == myStimType.aLF)
                HippStim.TL.b = HippStim.vars.T;
                HippStim.LRa.b = 1;
                
            elseif (stimType == myStimType.bL || ...
                    stimType == myStimType.bLF)
                HippStim.TL.b = HippStim.vars.T;
                HippStim.LRb.b = 1;
                
                
            elseif (stimType == myStimType.aBL || ...
                    stimType == myStimType.aBLF)
                HippStim.BL.b = 1;
                HippStim.TB.b = HippStim.vars.T;
                HippStim.TL.b = HippStim.vars.T;
                HippStim.BRa.b = 1;
                HippStim.LRa.b = 1;
                
            elseif (stimType == myStimType.bBL || ...
                    stimType == myStimType.bBLF)
                HippStim.BL.b = 1;
                HippStim.TB.b = HippStim.vars.T;
                HippStim.TL.b = HippStim.vars.T;
                HippStim.BRb.b = 1;
                HippStim.LRb.b = 1;
            end
            
        elseif t >= onset + 500 && t <= onset + 515
                
            if (stimType == myStimType.aBLF)
                HippStim.vars.T = 1;
%                 HippStim.BL.b = 1;
                HippStim.BF.b = 1;
                HippStim.LF.b = 1;
%                 HippStim.BRa.b = 1;
%                 HippStim.LRa.b = 1;
                HippStim.FRa.b = 1;
%                 HippStim.TB.b = HippStim.vars.T;
%                 HippStim.TL.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
            elseif (stimType == myStimType.bBLF)
                HippStim.vars.T = 1;
%                 HippStim.BL.b = 1;
                HippStim.BF.b = 1;
                HippStim.LF.b = 1;
%                 HippStim.BRb.b = 1;
%                 HippStim.LRb.b = 1;
                HippStim.FRb.b = 1;
%                 HippStim.TB.b = HippStim.vars.T;
%                 HippStim.TL.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
                
            elseif (stimType == myStimType.aBF)
                HippStim.vars.T = 1;
                HippStim.BF.b = 1;
%                 HippStim.BRa.b = 1;
                HippStim.FRa.b = 1;
%                 HippStim.TB.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
            elseif (stimType == myStimType.bBF)
                HippStim.vars.T = 1;
                HippStim.BF.b = 1;
%                 HippStim.BRb.b = 1;
                HippStim.FRb.b = 1;
%                 HippStim.TB.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
                
            elseif (stimType == myStimType.aLF)
                HippStim.vars.T = 1;
                HippStim.LF.b = 1;
%                 HippStim.LRa.b = 1;
                HippStim.FRa.b = 1;
%                 HippStim.TL.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
            elseif (stimType == myStimType.bLF)
                HippStim.vars.T = 1;
                HippStim.LF.b = 1;
%                 HippStim.LRb.b = 1;
                HippStim.FRb.b = 1;
%                 HippStim.TL.b = HippStim.vars.T;
                HippStim.TF.b = HippStim.vars.T;
                
                
            elseif (stimType == myStimType.aF)
                HippStim.FRa.b = 1;
                HippStim.TF.b = HippStim.vars.T;
                HippStim.vars.T = 1;
                
            elseif (stimType == myStimType.bF)
                HippStim.FRb.b = 1;
                HippStim.TF.b = HippStim.vars.T;
                HippStim.vars.T = 1;
                
            end
            
        end
        
        if stimType == myStimType.aNothing
            HippStim.TRa.b = HippStim.vars.T;
        elseif stimType == myStimType.bNothing
            HippStim.TRb.b = HippStim.vars.T;
        end
        

        % myStimType.BLF = 0;
        % myStimType.BF = 1;
        % myStimType.BL = 2;
        % myStimType.LF = 3;
        % myStimType.B = 4;
        % myStimType.L = 5;
        % myStimType.F = 6;
        % myStimType.Nothing = 7;


        %%
        
        tmpfields = fieldnames(HippStim);
        inputNoise = [];
        for tf = 1:numel(tmpfields)-1
            inputNoise = [inputNoise; (rand(floor((numNeurons/(numel(fieldnames(HippStim))-1))),1)*13-(13-(noiseAmplitude + (4 * HippStim.(tmpfields{tf}).b))))];
        end
        clear tf
    otherwise
        if noiseAmplitude > 0
            inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        end
end
clear numNeurons;
output.inputNoise = inputNoise;
output.HippStim = HippStim;
    
end