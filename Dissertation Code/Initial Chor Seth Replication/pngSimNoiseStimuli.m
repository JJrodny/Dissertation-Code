function output = pngSimNoiseStimuli(stimType, s1, s2, wholeNetwork, t, onset, offset, nucleiNumber)

% (taken from Ben St. Clair's code):
% Determine background noise. Must be between 0 and 1.
% This corresponds to the likelihood that a neuron will fire due to noise
% per second.
% bgNoise = [130,130,110,110,90,90,100];
% bgNoise = [130,0,110,110,0,0,0];
% The amplitude of the noise, in mV.
% noiseAmplitude = 6.5;
noiseAmplitude = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.noiseAmplitude;

numNeurons = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.N;
inputNoise = zeros(numNeurons,1);

switch nucleiNumber
    case wholeNetwork.nucIndex.PFC
        % PFC
        mynoiseAmplitude = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.noiseAmplitude;
        
%         inputNoise1 = mynoiseAmplitude * .8 + 4 * mod(round(t),1000)/1000;
%         inputNoise2 = mynoiseAmplitude * .8 + 4 * (1000 - mod(round(t),1000))/1000;
%         noiseAmplitude = min(inputNoise1,inputNoise2);
           
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        if t >= onset+offset && t < onset + 500+offset && stimType ~= 2 && stimType ~=3
            % S1 
            tmp = rand(s1,numNeurons/2,1)*13-(13-noiseAmplitude);
            inputNoise = [inputNoise(numNeurons/2+1:end);tmp];
            clear tmp;
            
        elseif t >= onset + 500+offset && t < onset + 1000+offset
            % S1 + S2
            tmp1 = rand(numNeurons/2,1)*13-(13-noiseAmplitude);
            tmp2 = rand(numNeurons/2,1)*13-(13-noiseAmplitude);
            if stimType ~= 1 && stimType ~= 4
                tmp1 = rand(s2,numNeurons/2,1)*13-(13-noiseAmplitude);
            end
            if stimType ~= 2 && stimType ~= 3
                tmp2 = rand(s1,numNeurons/2,1)*13-(13-noiseAmplitude);
            end
            
            inputNoise = [tmp1;tmp2];
            clear tmp1 tmp2;
            
        elseif t >= onset + 1000+offset && t < onset + 1500+offset && stimType ~= 1 && stimType ~= 4
            % S2
            tmp = rand(s2,numNeurons/2,1)*13-(13-noiseAmplitude);
            inputNoise = [tmp;inputNoise(numNeurons/2+1:end)];
            clear tmp;
        else
            % neither
        end
        clear numNeurons;
        
    case wholeNetwork.nucIndex.SENa
        % SEN A
        if t >= onset && t <= onset + 15 && stimType ~= 2 && stimType < 3
            noiseAmplitude = noiseAmplitude + 2.5;% * (t - onset + 1);
        end
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        
    case wholeNetwork.nucIndex.SENb
        % SEN B
        if t >= onset + 500 && t <= onset + 515 && stimType ~= 1 && stimType < 3
            noiseAmplitude = noiseAmplitude + 2.5;% * (t - (onset + 500) + 1);
        end
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        
    case wholeNetwork.nucIndex.STRstr
        % STR
%         noiseAmplitude = 6;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.INTa
        % INT A
%         noiseAmplitude = 6;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.INTb
        % INT B
%         noiseAmplitude = 6;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    case wholeNetwork.nucIndex.SNc
        % DA
%         noiseAmplitude = 6;
        inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
%     case wholeNetwork.nucIndex.HIPP
%         % HIPP
% %         noiseAmplitude = 6;
%         inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
    otherwise
        if noiseAmplitude > 0
            inputNoise = rand(numNeurons,1)*13-(13-noiseAmplitude);
        end
end
clear numNeurons;
output.inputNoise = inputNoise;
    
end