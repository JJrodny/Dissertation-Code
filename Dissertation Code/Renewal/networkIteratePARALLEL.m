
for nucleiNumber = 1:wholeNetwork.numNuclei
    % assume network is called: wholeNetwork.
% the nuclei to iterate within wholeNetwork is: nucleiNumber

%membranePotential =>
%       wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential
%
% these scripts and functions were written by:
% William Benjamin St. Clair wst.clair@ucmerced.edu
% over a period from 02/2010-11/2012
% Further modified by Jeffrey Rodny jrodny@ucmerced.edu 2013-2017
% Here, Jeff Rodny modified mostly those labeled: 'JJ CODE' among others


%%%% Find action potentials.
fired=find(wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential>=30);    % indices of spikes

% For observation, we record the occurrence of the firings.
wholeNetwork.firings = [wholeNetwork.firings; t+0*fired,fired, nucleiNumber+0*fired];

% Instantiate threshold effects, setting membrane potential to c and
% adding d to the equalibrium force.
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential(fired)=wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.c(fired);
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce(fired)=wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce(fired)+wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.d(fired);

% Since the Polychronization paper uses nonplastic inhibitory
% connections, we only should change weight values for excitatory
% neurons that fired. Hence, we obtain the excitatory subset of the
% fired neurons for this purpose.
excitatorySubset = fired(fired <= wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.Ne);

% Calculate STDP for connections in the nuclei based on new firings.
% Calculate STDP for all outgoing connections.
for receivingNuclei = 1:wholeNetwork.numNuclei
    % Only attempt to calculate STDP when there exists some connectivity!
    %% JJ CODE
    % Chorley and Seth 2012 implementation
    if wholeNetwork.plasticConnections(nucleiNumber,receivingNuclei) == 1 ||...
            wholeNetwork.plasticConnections(receivingNuclei,nucleiNumber) == 1

        if ~isempty(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire) ...
         || ~isempty(wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.lastFire)
            % STDP applies spike time dependent plasticity
            [stdpEffectsPos, stdpEffectsNeg] = ...
                STDP(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire, ...
                excitatorySubset, ...
                t, ...
                wholeNetwork.nuclei{receivingNuclei}{receivingNuclei}.params.Ne, ...
                wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.lastFire);

            if ~isempty(wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.lastFire)
                % Apply STDP for excitatory connections in the positive
                % time case.
                % Maximum weight size is an indication of the hypothetical minimum
                % number of incoming spikes required to cause a postsynaptic spike.
                % Clip weights to between upper bounds ...
                wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.S(wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.S ...
                    > wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.params.weightUpperbound) ...
                    = wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.params.weightUpperbound;
               
                %% JJ
               
                wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.eligibilityTrace(:,excitatorySubset) = ...
                    wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.eligibilityTrace(:,excitatorySubset) + ...
                    stdpEffectsPos(:,excitatorySubset);
               
                % JJ end
            end
            if ~isempty(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire)
                % Apply STDP for excitatory connections in the negative
                % time case.
                
                % ... and lower bounds.
                wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S ...
                    < wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.params.weightLowerbound) = ...
                    wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.params.weightLowerbound;
               
                %% JJ
               
                wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace(excitatorySubset,:) = ...
                    wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace(excitatorySubset,:) + ...
                    stdpEffectsNeg(excitatorySubset,:);
               
                % JJ end
            end

            if ~isempty(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire)
                % Don't let excitatory synapses become inhibitory synapses by
                % STDP. This process does not need to be repeated for inhibitory
                % synapses because STDP should not be applied to inhibitory.
                % It does not need to applied to the nuclei{receivingNuclei}{nucleiNumber} because it is not receiving any negative STDP, so it could not have just become inhibitory.
                wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S(1:wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.Ne,1:wholeNetwork.nuclei{receivingNuclei}{receivingNuclei}.params.Ne) < 0) = 0;
            end

            %% JJ CODE
            % calculate eligibility traces
            if ~isempty(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire)
                % CHORLEY SETH
                if (nucleiNumber == wholeNetwork.nucIndex.STRstr && receivingNuclei == wholeNetwork.nucIndex.PFC) || ...
                        (nucleiNumber == wholeNetwork.nucIndex.PFC && receivingNuclei == wholeNetwork.nucIndex.STRstr)
                    wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace = ...
                        wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace ...
                        * (1 - (1/200)); %0.2s = 200ms
                else
                    wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace = ...
                        wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace ...
                        * (1 - (1/1000)); % 1s = 1000ms
                end
                   
                % make the weight changes
                if wholeNetwork.plasticConnections(nucleiNumber,receivingNuclei) == 1
                        
                    wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S = ...
                        wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.S ...
                        + 0.05*((wholeNetwork.DA^2 ...
                        .* wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.eligibilityTrace) ...
                        / 5);
                end
            end
            %% end JJ
        end
    end
end

% This is code remaining from a gating value implementation.
%   It would need modification to be used,
%   so it remains to make it easier to implement gating values
%   in the future.
%   This is part of the model conditions required to exhibit
%   the izhi-szat working memory conditions.
% Apply the STDP effects to the short-term STDP mechanism.
%          gatingValues(:,excitatorySubset) = gatingValues(:,excitatorySubset) ...
%              + stdpEffects(:,excitatorySubset)*100;
% ASSUMPTION: Gating values should not go below zero.
%          gatingValues(gatingValues < 0) = 0;



if ~isempty(fired)
  % Every spike causes N potentials. Each iteration of this for
  % loop adds N spike seperate potentials to conductingPotentials.
    for k = 1:length(fired)
        for receivingNuclei = 1:wholeNetwork.numNuclei
            % When a neuron fires, all of its dendrites immediately know, which is relevant for STDP.
            % Hence, we must update all of their lastFires.
            % We loop through each possible sending nuclei.
            % Only update lastFires where there actually exist connections to the nuclei of interest! (nucleiNumber)
            % For this statement, we use receivingNuclei as if it is the "sendingNuclei", so the self-documenting variable name is misleading for this statement.
            if ~isempty(wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.lastFire)
                % There exist connections, so update the relevant lastFire entries.
                wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.lastFire(wholeNetwork.nuclei{receivingNuclei}{nucleiNumber}.connectivityMatrix(:,fired(k)) == 1,fired(k)) = t;
            end
           
            if ~isempty(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.lastFire)
                % Determine the relevant outgoing connections.
                effectiveConnections = find(wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.connectivityMatrix(fired(k),:));
                % Assemble a list of events for all nonzero outgoing
                % connections.
                newPotentials = [t*ones(length(effectiveConnections),1)  ... %the current time + how long it takes to send is
                  + wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.conductanceDelays(fired(k),effectiveConnections)', ...%% BEN+JJ EDIT - remove the '-1' - 1, ... %when they will arrive
                ones(length(effectiveConnections),1)*fired(k), ... %where they came from (neuron number)
                transpose(effectiveConnections)]; % where they are going (neuron number)

                % Combine the new list with the old list.
                wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.conductingPotentials = [wholeNetwork.nuclei{nucleiNumber}{receivingNuclei}.conductingPotentials; newPotentials];
            end
        end
    end
end


% If no potentials have arrived, input will stay zero. Otherwise,
% it will grow.
input = zeros(wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.N,1);
% Find all received potentials for each sending nuclei
for sendingNuclei = 1:wholeNetwork.numNuclei
    if ~isempty(wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.lastFire)
        if ~isempty(wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials)
            % row indices of all potentials that have arrived at their
            % destination.
            arrivedPotentials = wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(:,1) == t;
            if ~isempty(arrivedPotentials)
                % We create a potentialList to minimize the size of the index.
                % Since conductingPotentials can grow quite large, we want to
                % only reference it the minimum number of times. Particularly, I
                % have noticed that using syntax like
                % conductingPotentials(arrivedPotentials,:) is very
                % computationally expensive and has enhanced huge slow downs.
                potentialList = wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(arrivedPotentials,:);
                % for each arrived potential
                for potIndex = 1:size(potentialList,1)
                    % The potential is being received by a neuron, so update the lastFire term.
                    wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.lastFire(potentialList(potIndex,2), ...
                           potentialList(potIndex,3)) = t;

                    % Calculate the input from the potential.
                    % In izhi&szat, I = S*(1+sd), where sd is a term that
                    % grows and shrinks by a form of STDP, and has constant
                    % decay.
                    % without gatingValues
                    input(potentialList(potIndex,3)) = ...
                      + input(potentialList(potIndex,3)) ...
                      + wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.S(potentialList(potIndex,2), potentialList(potIndex,3));
                end
                % Clip out arrived action potentials.
                wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials = wholeNetwork.nuclei{sendingNuclei}{nucleiNumber}.conductingPotentials(~arrivedPotentials, :);
            end
        end
    end
end

% include thalamic input
tmp = pngPARALLELSimNoiseStimuli(stimType, s1, s2, s3, wholeNetwork, t, onset, offset, nucleiNumber, HippStim, myStimType);
HippStim = tmp.HippStim;
input=input + tmp.inputNoise;
clear tmp

% Membrane Potential Update.
% step 0.5 ms for numerical stability
%  v=v+0.5*((0.04*v+5).*v+140-u+I);
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential + ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.timeStep * ...
    ((0.04*wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential+5).*...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential + 140 - ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + input);
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential = wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential + ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.timeStep * ...
    ((0.04*wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential+5).*...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential + 140 - ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + input);
%% JJ EDIT
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential = min(wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential,140);
%% JJ END

% u=u+a.*(0.2*v-u);
wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce = ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce + ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.a.* ...
    (wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.params.b.*...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.membranePotential - ...
    wholeNetwork.nuclei{nucleiNumber}{nucleiNumber}.equalibriumForce);

% This accomplishes an exponential decay which nears 0.0183 (close to
% zero) after 5000 iterations. 0.0183 = exp(5000 * log(.9992)). I
% chose the factor .9992 visually after plotting many points, so as
% to best fit the description in izhi&szat: "it decays back to 0 with
% a time constant 5 seconds."
%          gatingValues = gatingValues*.9992;
% Gating Values have a maximum of 100% increase of excitation.
% "A maximum of 100% temporary increase relative to baseline"
%          gatingValues(gatingValues > 1) = 1;

end

%% JJ CODE DOPAMINE
% get the alpha for the DA nuclei (nuclei 7)

% alpha is step-increased by 0.05 uM for each spike of a DA neuron
% while diffusing with exponential time constant t=0.1s
% side note, a baseline DA concentration of between 0.5 and 1uM is
% maintained by the background DA activity
% when it does increase, the concentration will go up to 2 or 3 uM

tmp = wholeNetwork.firings(wholeNetwork.firings(:,1)==t,:);
firingsNuclei = tmp(tmp(:,3)==wholeNetwork.nucIndex.SNc,1:2);

if ~isempty(firingsNuclei)
    wholeNetwork.DA = wholeNetwork.DA + 0.05*(length(firingsNuclei));
    if wholeNetwork.DA > 10
        wholeNetwork.DA = 10;
        disp('DA Maxed out');
    end
end

wholeNetwork.DA = wholeNetwork.DA .* (1 - 1/(100)); % exponential time constant 0.1s (Chorley, Seth)

%% now use this value to modulate excitability of STR neuorns
% wholeNetwork.nuclei{wholeNetwork.nucIndex.STRddi}{wholeNetwork.nucIndex.STRddi}.params.b = ...
%     (0.21 - 0.0025*wholeNetwork.DA^2);
% wholeNetwork.nuclei{wholeNetwork.nucIndex.STRdi}{wholeNetwork.nucIndex.STRdi}.params.b = ...
%     (0.195 + 0.005*wholeNetwork.DA^2);
wholeNetwork.nuclei{wholeNetwork.nucIndex.STRstr}{wholeNetwork.nucIndex.STRstr}.params.b = ...
    (0.19 + 0.01*wholeNetwork.DA^2);

tmpDA(t) = wholeNetwork.DA;

%% GRAPH IT! (JJ CODE)
if t == onset + 1900 %% || (t < 15000 && t > 1999 && mod(t,100) == 0) %use this to make a gif
    XaxisLength = 2000;
    % This figure line opens a maximized figure
    figureID = figure('units','normalized','outerposition',[0 0 1 1]);
    firingsNuclei = cell(1,3);
   

    for ii = 1:wholeNetwork.numNuclei
        tmp = wholeNetwork.firings(wholeNetwork.firings(:,1)>t-XaxisLength,:);
        firingsNuclei{ii} = tmp(tmp(:,3)==ii,1:2);
               
        switch ii
            case wholeNetwork.nucIndex.PFC
                %PFC 2
                subplot(3,3,1);
            case wholeNetwork.nucIndex.STRstr
                % STRstr 15
                subplot(3,3,2);
            case wholeNetwork.nucIndex.SENa
                subplot(6,3,12);
                % SENa 8
            case wholeNetwork.nucIndex.SENb
                subplot(6,3,18);
                % SENb 9
            case wholeNetwork.nucIndex.SENc
                subplot(6,3,15);
            case wholeNetwork.nucIndex.INTc
                continue;
            case wholeNetwork.nucIndex.INTa
                continue;
                % INTa 10
            case wholeNetwork.nucIndex.INTb
                continue;
                subplot(9,3,17);
                % INTb 11
            case wholeNetwork.nucIndex.SNc
                % SNc 13
                subplot(2,3,3)
            case wholeNetwork.nucIndex.HIPP
                subplot(3,3,[4,7]);
            case wholeNetwork.nucIndex.INTh
                subplot(3,3,[5,8]);
            otherwise
                continue;
        end
        
        numUnits = length(wholeNetwork.nuclei{ii}{ii}.S);
        % (-0.048*numUnits+52) is a formula just to give good sized dots
        
        % Color points in Scatter plot differently for when thinking of
        % Bell or Food
        % get the points where the X value is after onset, and the ones
        % where x is before end of onset+1000, == will get the points
        % between.
        % Then == that to the Y values between 500 and 1000
        % same for Food, but X t is moved over by 500 and down by 500
        myBellScatterPoints = ((firingsNuclei{ii}(:,1)>onset+offset) .* (firingsNuclei{ii}(:,1)<onset+offset+1000) == 1) .* ((firingsNuclei{ii}(:,2)>2*floor(numUnits/3)) .* (firingsNuclei{ii}(:,2)<=numUnits) == 1) == 1;
        myFoodScatterPoints = ((firingsNuclei{ii}(:,1)>onset+offset+500) .* (firingsNuclei{ii}(:,1)<onset+offset+1500) == 1) .* ((firingsNuclei{ii}(:,2)>0) .* (firingsNuclei{ii}(:,2)<=floor(numUnits/3)) == 1) == 1;
        myLightScatterPoints = ((firingsNuclei{ii}(:,1)>onset+offset) .* (firingsNuclei{ii}(:,1)<onset+offset+1000) == 1) .* ((firingsNuclei{ii}(:,2)>floor(numUnits/3)) .* (firingsNuclei{ii}(:,2)<=2*floor(numUnits/3)) == 1) == 1;
        
        if stimType == myStimType.aBLF ...
                || stimType == myStimType.aBL ...
                || stimType == myStimType.aBF ...
                || stimType == myStimType.aB ...
                || stimType == myStimType.bBLF ...
                || stimType == myStimType.bBL ...
                || stimType == myStimType.bBF ...
                || stimType == myStimType.bB 
            myBellColor = [0,1,0]; % RGB = Blue
        else
            myBellColor = [0,0,0];
        end
            
        myBellScatterPointsColor = ones(length(myBellScatterPoints),3).*[myBellScatterPoints,myBellScatterPoints,myBellScatterPoints];
        myBellScatterPointsColor = myBellScatterPointsColor.*[(...
            ones(length(myBellScatterPoints),1).*myBellColor(1)),...
            ones(length(myBellScatterPoints),1).*myBellColor(2),...
            ones(length(myBellScatterPoints),1).*myBellColor(3)];
        
        if stimType == myStimType.aBLF ...
                || stimType == myStimType.aBF ...
                || stimType == myStimType.aLF ...
                || stimType == myStimType.aF ...
                || stimType == myStimType.bBLF ...
                || stimType == myStimType.bBF ...
                || stimType == myStimType.bLF ...
                || stimType == myStimType.bF
            myFoodColor = [1,0,0]; % RGB = Green
        else
            myFoodColor = [0,0,0];
        end
        myFoodScatterPointsColor = ones(length(myFoodScatterPoints),3).*[myFoodScatterPoints,myFoodScatterPoints,myFoodScatterPoints];
        myFoodScatterPointsColor = myFoodScatterPointsColor.*[(...
            ones(length(myFoodScatterPoints),1).*myFoodColor(1)),...
            ones(length(myFoodScatterPoints),1).*myFoodColor(2),...
            ones(length(myFoodScatterPoints),1).*myFoodColor(3)];
        
        
        if stimType == myStimType.aBLF ...
                || stimType == myStimType.aBL ...
                || stimType == myStimType.aLF ...
                || stimType == myStimType.aL ...
                || stimType == myStimType.bBLF ...
                || stimType == myStimType.bBL ...
                || stimType == myStimType.bLF ...
                || stimType == myStimType.bL
            myLightColor = [1,0,1]; % RGB = Green
        else
            myLightColor = [0,0,0];
        end
        myLightScatterPointsColor = ones(length(myLightScatterPoints),3).*[myLightScatterPoints,myLightScatterPoints,myLightScatterPoints];
        myLightScatterPointsColor = myLightScatterPointsColor.*[(...
            ones(length(myLightScatterPoints),1).*myLightColor(1)),...
            ones(length(myLightScatterPoints),1).*myLightColor(2),...
            ones(length(myLightScatterPoints),1).*myLightColor(3)];
        
        
        myScatterPointsColor = myFoodScatterPointsColor + myBellScatterPointsColor + myLightScatterPointsColor;
        
        if ii == wholeNetwork.nucIndex.PFC
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,myScatterPointsColor);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.SENa
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[0 1 0],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.INTa
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[.75 1 .75],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.SENb
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.INTb
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[1 .75 .75],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.SENc
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[1 0 1],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.INTc
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[1 .75 1],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        elseif ii == wholeNetwork.nucIndex.SNc
            % % Data to plot.
            x1 = firingsNuclei{ii}(:,1);
            x2 = max(t-XaxisLength,1):t;
            y1 = firingsNuclei{ii}(:,2);
            y2 = tmpDA(max(t-XaxisLength,1):t)';
            AX = subplot(2,3,3);
            
            P = get(AX,'pos');    % Get the position.
            delete(AX)            % Delete the subplot axes
            subplot(2,3,6);
            [AX,H1,H2] = plotyy(x1,y1,x2,y2,'scatter','plot');
            set(AX,'pos',P)       % Recover the position.
            axis(AX(1),[t-XaxisLength,t,0,numUnits]);
            set(AX(1),'ytickmode','auto','FontSize',12);
            set(H1,'SizeData',1000/numUnits,'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor',[0 0 0]);
            set(H2,'LineWidth',2);
            
            ylabel(AX(2),'DA Concentration (uM)','FontSize',20);
            axis(AX(2),[t-XaxisLength,t,0,5]);
            set(AX(2),'ytickmode','auto','xtick',[],'FontSize',20);
            
        elseif ii == wholeNetwork.nucIndex.HIPP
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        else
            scatter(firingsNuclei{ii}(:,1),firingsNuclei{ii}(:,2),1000/numUnits,'MarkerFaceColor',[.9 .9 .9],'MarkerEdgeColor',[0 0 0]);%max(4,floor(-0.048*numUnits+52)));
        end
        
        exciteNumber = find(wholeNetwork.nuclei{ii}{ii}.S<0,1)-1;
        if exciteNumber > 0
            line([t-XaxisLength,t],[exciteNumber,exciteNumber],'Color','k','LineStyle',':');
        end
        if ii == wholeNetwork.nucIndex.PFC
            if stimType == myStimType.aBLF ...
                    || stimType == myStimType.aBL ...
                    || stimType == myStimType.aBF ...
                    || stimType == myStimType.aB ...
                    || stimType == myStimType.bBLF ...
                    || stimType == myStimType.bBL ...
                    || stimType == myStimType.bBF ...
                    || stimType == myStimType.bB 
                line([onset+offset,     onset+offset],      [2*floor(numUnits/3),   numUnits],              'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+1000+offset,onset+1000+offset], [2*floor(numUnits/3),   numUnits],              'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+offset,     onset+1000+offset], [2*floor(numUnits/3),   2*floor(numUnits/3)],   'Color',[.1,.1,.1],'LineStyle','-.');
            end
            if stimType == myStimType.aBLF ...
                    || stimType == myStimType.aBL ...
                    || stimType == myStimType.aLF ...
                    || stimType == myStimType.aL ...
                    || stimType == myStimType.bBLF ...
                    || stimType == myStimType.bBL ...
                    || stimType == myStimType.bLF ...
                    || stimType == myStimType.bL
                line([onset+offset, onset+offset],          [floor(numUnits/3),     2*floor(numUnits/3)],   'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+1000+offset,onset+1000+offset], [floor(numUnits/3),     2*floor(numUnits/3)],   'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+offset, onset+1000+offset],     [floor(numUnits/3),     floor(numUnits/3)],     'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+offset, onset+1000+offset],     [2*floor(numUnits/3),   2*floor(numUnits/3)],   'Color',[.1,.1,.1],'LineStyle','-.');
            end
            if stimType == myStimType.aBLF ...
                    || stimType == myStimType.aBF ...
                    || stimType == myStimType.aLF ...
                    || stimType == myStimType.aF ...
                    || stimType == myStimType.bBLF ...
                    || stimType == myStimType.bBF ...
                    || stimType == myStimType.bLF ...
                    || stimType == myStimType.bF
                line([onset+500+offset,    onset+500+offset], [0,                 floor(numUnits/3)], 'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+1500+offset,    onset+1500+offset], [0,                 floor(numUnits/3)], 'Color',[.1,.1,.1],'LineStyle','-.');
                line([onset+offset+500,    onset+1500+offset], [floor(numUnits/3), floor(numUnits/3)], 'Color',[.1,.1,.1],'LineStyle','-.');
            end
        elseif ii == wholeNetwork.nucIndex.STRddi || ii == wholeNetwork.nucIndex.STRdi || ii == wholeNetwork.nucIndex.STRstr
            line([onset+5+offset,onset+5+offset],[0,numUnits],'Color','k','LineStyle','--');
            line([onset+10,onset+10],[0,numUnits],'Color','k','LineStyle','--');
            line([onset+505+offset,onset+505+offset],[0,numUnits],'Color','k','LineStyle','--');
            line([onset+510,onset+510],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.SENa
            line([onset,onset],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.INTa
            line([onset+5,onset+5],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.SENb
            line([onset+500,onset+500],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.INTb
            line([onset+505,onset+505],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.SENc
            line([onset+0,onset+0],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.INTc
            line([onset+5,onset+5],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.SNc
            line([onset+10,onset+10],[0,numUnits],'Color','k','LineStyle','--');
            line([onset+510,onset+510],[0,numUnits],'Color','k','LineStyle','--');
        elseif ii == wholeNetwork.nucIndex.HIPP
            
            for hippIndex = 1:14
                line([t-XaxisLength,t],[hippIndex*floor(numUnits/(numel(fieldnames(HippStim))-1)),hippIndex*floor(numUnits/(numel(fieldnames(HippStim))-1))],'Color',[.1,.1,.1],'LineStyle','-');
            end
            
        elseif ii == wholeNetwork.nucIndex.INTh
            
            
            for hippIndex = 1:14
                line([t-XaxisLength,t],[hippIndex*floor(numUnits/(numel(fieldnames(HippStim))-1)),hippIndex*floor(numUnits/(numel(fieldnames(HippStim))-1))],'Color',[.1,.1,.1],'LineStyle','-');
            end
            
        end
        H = gca;
        set(H,'xtick', [floor((t-XaxisLength)/1000)*1000:1000:ceil((t)/1000)*1000]);
        set(H,'xticklabel', [floor((t-XaxisLength)/1000)*1000:1000:ceil((t)/1000)*1000]./1000);
        if ii ~= wholeNetwork.nucIndex.SENb && ...
                ii ~= wholeNetwork.nucIndex.INTb && ...
                ii ~= wholeNetwork.nucIndex.PFC && ...
                ii ~= wholeNetwork.nucIndex.STRstr && ...
                ii ~= wholeNetwork.nucIndex.HIPP && ...
                ii ~= wholeNetwork.nucIndex.INTh% && ...
                %ii ~= wholeNetwork.nucIndex.SNc
            set(gca,'XTickLabel',[]);
        elseif ii ~= wholeNetwork.nucIndex.PFC && ...
                ii ~= wholeNetwork.nucIndex.STRstr% && ...
                %ii ~= wholeNetwork.nucIndex.SENb && ...
                %ii ~= wholeNetwork.nucIndex.INTb
            xlabel('Time (seconds)','FontSize',12);
        end

        if ii == wholeNetwork.nucIndex.PFC %|| ...
                %ii == wholeNetwork.nucIndex.STRstr || ...
                %ii == wholeNetwork.nucIndex.SNc
            ylabel('Neuron Number','FontSize',15);
        elseif ii == wholeNetwork.nucIndex.SENb || ...
                ii == wholeNetwork.nucIndex.SENb|| ...
                ii == wholeNetwork.nucIndex.SENb %|| ...
                %ii == wholeNetwork.nucIndex.INTb || ...
                %ii == wholeNetwork.nucIndex.INTa
%                 ylabel('Neuron No.','FontSize',15);
        end
        if ii == wholeNetwork.nucIndex.HIPP
            ylabel('Neuron Group','FontSize',15);
        end

        if ii == wholeNetwork.nucIndex.HIPP
            H = gca;
            set(H,'ytick', [(numUnits/(numel(fieldnames(HippStim))-1))/2:(numUnits/(numel(fieldnames(HippStim))-1)):numUnits]);
            set(gca,'FontSize',8)
            set(H,'yticklabel', {'BL','BF','BRa','BRb','LF','LRa','LRb','FRa','FRb','RaRb','TB','TL','TF','TRa','TRb'});
            set(gca,'FontSize',8)
            axis([t-XaxisLength,t,1,numUnits]);
        end

        axis([t-XaxisLength,t,0,max(0,numUnits)]);
        title(wholeNetwork.nuclei{ii}{ii}.params.Name,'FontSize',7);
        set(gca,'FontSize',20); % Axis number sizes
    end
    try
        print(figureID,['data/',myIdentifier,'/tb',num2str(t),'.png'],'-dpng');
    catch
        disp('something went wrongfdsfsdfdsfsdfsd');
        saveas(figureID,['data/',myIdentifier,'/tb',num2str(t),'2.png'],'png');
    end
    close(figureID);
    
    
    
    %% Find the MAx DA spike level and graph it
    tmpMaxDA1 = max(tmpDA(onset+20:onset+70));
    tmpMaxDA2 = max(tmpDA(onset+520:onset+570));
    wholeNetwork.DA1avg = [wholeNetwork.DA1avg,tmpMaxDA1];
    wholeNetwork.DA2avg = [wholeNetwork.DA2avg,tmpMaxDA2];
    clear tmpMaxDA1 tmpMaxDA2;
    % 'convolve' the data to a fractional dataset that adds up to one, in
    % order to smooth out the curve
    NNN = 20;
    figID = figure(28);
    subplot(2,1,1);
    plot([1:length(wholeNetwork.DA1avg)],wholeNetwork.DA1avg,[NNN/2:length(wholeNetwork.DA1avg)-NNN/2],conv(wholeNetwork.DA1avg,[zeros(NNN/2,1)/NNN;ones(NNN/2,1)/NNN*2],'valid'))
    title('Average Max DA spike value for Light');
    
    if t >= timeBLFStart
        tmpT = floor(timeBLFStart/3000);
        line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
        clear tmpT;
    end
    if t >= timeBorLStart
        tmpT = floor(timeBorLStart/3000);
        line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
        clear tmpT;
    end
    
%     if t >= timePostExtinctionStart
%         tmpT = floor(timePostExtinctionStart/3000);
%         line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
%         clear tmpT;
%     end
    
    subplot(2,1,2);
    plot([1:length(wholeNetwork.DA2avg)],wholeNetwork.DA2avg,[NNN/2:length(wholeNetwork.DA2avg)-NNN/2],conv(wholeNetwork.DA2avg,[zeros(NNN/2,1)/NNN;ones(NNN/2,1)/NNN*2],'valid'))
    title('Average Max DA spike value for Food');
    xlabel('trial');
    
    if t >= timeBLFStart
        tmpT = floor(timeBLFStart/3000);
        line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
        clear tmpT;
    end
    if t >= timeBorLStart
        tmpT = floor(timeBorLStart/3000);
        line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
        clear tmpT;
    end
    
%     if t >= timePostExtinctionStart
%         tmpT = floor(timePostExtinctionStart/3000);
%         line([tmpT,tmpT],[0,5],'Color',[.1,.1,.1],'LineStyle','--');
%         clear tmpT;
%     end
    
    saveas(figID,['data/',myIdentifier,'/DAoverTime.png'],'png');
    close(figID);
    clear NNN figID;
    
end
%% End JJ CODE

clear tmp firingsNuclei
   
clear receivingNuclei

clear firingsNuclei tmp numUnits XaxisLength exciteNumber

clear arrivedPotentials effectiveConnections excitatorySubset fired
clear k newPotentials potIndex potentialList sendingNuclei
clear stdpEffectsNeg stdpEffectsPos input