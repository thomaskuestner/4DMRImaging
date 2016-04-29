function [dKSpace, dRespPhases] = fRetroGating(sFilename, iNPhases, dTolerance, dTime, sGatingMode, sGatingGather, sPrecision, iEcho, lDebug)
%FRETROGATING
%
%   Copyright 2014-2016 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   and Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

% inputs
% sFilename         path to measurement file
% iNPhases          number of respiratory gates
% dTolerance        view sharing blending factor b of gates (1 (non-overlapping) <= b <= 2 (completely overlapping))
% dTime             crop scan time (0 := no cropping | scan time [s])
% sGatingMode       select gating procedure
%                   'percentile'    select gate centroids according to 5% and 95% percentile
%                   'kmeans'        perform k-means clustering to select gate centroids
% sGatingGather     sample gathering mode: 
%                   'closest'   take closest sample to gate centroid
%                   'average'   take all possible samples and weight them
%                               according to their distance to the gate centroid
%                   'collect'   gather all samples and let the
%                               reconstruction decide for the most favourable one
%                   'none'      no k-Space population, show simulations of the k-Space sampling density
% iEcho             current to be processed echo
% lDebug            debugging flag for additional plots
% 
% outputs
% dKSpace           gated kSpace
% dRespPhases       respiratory phases (TR level)

if nargin < 9, lDebug = false; end
if nargin < 8, iEcho = 1; end
if nargin < 7, sPrecision = 'single'; end
if nargin < 6, sGatingGather = 'closest'; end
if nargin < 5, sGatingMode = 'percentile'; end
if nargin < 4, dTime = 0; end
if nargin < 3, dTolerance = 1; end
if nargin < 2, iNPhases = 1; end

% parse inputs 
[sPath, sName, sExt] = fileparts(sFilename);
if isempty(sExt), sFilename = [sFilename, '.dat']; end
% sFilename = which(sFilename);
[sPath, sName] = fileparts(sFilename);
sFilename = [sPath, filesep, sName];
iNPhasesLoop = 1:iNPhases;

% -------------------------------------------------------------------------
% Get some relevant data from the drecksMDH
drecksMDH = fMeas_mainLUTDrecks( sFilename );
iNPartitions    = drecksMDH.Geo.MatrixSize(3);
iNLines         = drecksMDH.Geo.MatrixSize(2);
iBaseRes        = drecksMDH.Geo.MatrixSize(1);
dTR             = drecksMDH.Contrast.TR./1E6;
dNavPeriod      = drecksMDH.Wip.NavPeriod;
if(isfield(drecksMDH.Geo,'EchoPosition'))
    iEchoLine       = drecksMDH.Geo.EchoPosition(1);
    iEchoPartition  = drecksMDH.Geo.EchoPosition(2);
else
    iEchoLine       = drecksMDH.Geo.MatrixSize(2)/2;
    iEchoPartition  = drecksMDH.Geo.MatrixSize(3)/2;
end   
if(isfield(drecksMDH,'Wip'))
    iNavPERes       = drecksMDH.Wip.NavRes(1);
    iNav3DRes       = drecksMDH.Wip.NavRes(2);
else
    iNavPERes       = 8;
    iNav3DRes       = 1;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Get loop counters from file and reduce to first channel and first echo
load([sFilename, '.mat'], 'iLC');
load([sFilename, '.mat'], 'iSP');
load([sFilename, '.mat'], 'dNavInt');
iNChannels = max(unique(double(iLC(:,2))));

lMask = iLC(:, 1) == 2*iBaseRes & iLC(:, 11) < 65000 & iLC(:,7) == iEcho-1; % throw out dummy measurements (which are needed to capture surrogate signals and keep steady state)
iSP = iSP(lMask);
iLCReduced = iLC(lMask, :); % only "real" readouts
iSP = iSP(1:iNChannels:end,:);
iLCReduced = iLCReduced(1:iNChannels:end, :);

fprintf(1, 'Scantime was %f3 s\n', length(iLCReduced)*dTR);

% -------------------------------------------------------------------------
% check cropping parameter dTime
dPRECROP = 5.0; % Crop first 5 s
fprintf(1, 'Cropped initial %.3f s for steady-state!\n', dPRECROP);
if(isscalar(dTime))
    dTimeStart = dPRECROP;
    dTimeEnd = dTime;
    dTime = num2str(dTime,'%03d');
else
    dTimeStart = dTime(:,1);   
    dTimeEnd = dTime(:,2);
    dTime = dTime.'; % for backward compatibility
    dTime = strrep(num2str(dTime(:).','%03d'),' ','');
end

if(dTimeStart(1) < dPRECROP)
    dTimeStart(1) = dPRECROP;
end
if(any(dTimeEnd == 0))
   dTimeEnd(dTimeEnd == 0) = length(iLCReduced)*dTR; 
end
if(dTimeEnd(end) > length(iLCReduced)*dTR)
    dTimeEnd(end) = length(iLCReduced)*dTR;
end

iCropRange = [];
for iI = 1:length(dTimeStart)
    iCropRange = [iCropRange,round([dTimeStart(iI)/dTR:dTimeEnd(iI)/dTR])];
    fprintf(1, 'Scantime cutout from %.3f s to %.3f s\n', dTimeStart(iI),dTimeEnd(iI));
end
iCropRange = unique(iCropRange);

iLCReduced  = iLCReduced(iCropRange,:);
iSP         = iSP       (iCropRange,:);
dNavInt     = dNavInt   (iCropRange);

iLine = iLCReduced(:,3);
iPartition = iLCReduced(:,6);

iADCLength = 128/4 + iBaseRes*4; % MDH + ADC

if(usejava('jvm') && ~feature('ShowFigureWindows'))
    % use text-based alternative
    lflagDisp = false;
else
    lflagDisp = true;
end
  
% -------------------------------------------------------------------------
% respiratory gating
% -------------------------------------------------------------------------
if(all(iNPhasesLoop > 0) && exist('dNavInt', 'var'))
    % Calculate the center positions of the phases
    % ---------------------------------------------------------------------
    % distinguish between inhale and exhale
    iNPhasesOld = iNPhases;
   
    % gating method
    dNavInt = dNavInt(1:length(iLine));
    dNavMin = min(dNavInt);
    dNavMax = max(dNavInt);
    dToleranceOld = dTolerance;
    if(dNavMin == dNavMax) % BH acquisition
        dGatePos = dNavMin;
        iNPhases = 1;
    else    
        switch sGatingMode
            case 'percentile'
                % 5% and 95% percentile gating
                dX = linspace(dNavMin, dNavMax, 256);
                iHist = hist(dNavInt, dX);
                iSum = sum(iHist);

                iLowerLim = 1;
                while sum(iHist(iLowerLim:end))/iSum > 0.9
                    iLowerLim = iLowerLim + 1;
                end
                iUpperLim = 256;
                while sum(iHist(1:iUpperLim))/iSum > 0.9
                    iUpperLim = iUpperLim - 1;
                end

                % equally spaced gate centroids
                dMin = dX(iLowerLim);
                dMax = dX(iUpperLim);
                dGatePos = linspace(dMax, dMin, iNPhasesOld)';
                dTolerance = abs(diff(dGatePos)).*dTolerance./2;
                if(isempty(dTolerance)) % breathhold
                    dTolerance = dToleranceOld;
                end
                dTolerance = repmat(dTolerance(1),[iNPhases, 1]);
                
            case 'kmeans'
                % kmeans clustering of gate centroids
                dGrad = gradient(dNavInt);
                lExhale = dGrad >= 0;
                lInhale = dGrad < 0;
                iInd = find(lExhale,1,'first'):find(lInhale,1,'last'); % without outliers (due to dNav extraction)
                lInd = false(1,length(dNavInt));
                lInd(iInd) = true;
                [iClusterEx, dGatePosEx] = kmeans(dNavInt(lInd).',iNPhases); % without outliers (due to dNav extraction)
                dGatePos = sort(dGatePosEx,'descend');

                % get tolerance window width
                dMinDist = zeros(length(dGatePos),2); % 1st col: upper dist, 2nd col: lower dist
                for iI = 1:length(dGatePos)
                    lIndCl = lInd;
                    lIndCl(lInd) = iClusterEx == find(dGatePos(iI) == dGatePosEx);
                    dMinDist(iI,1) = max(dNavInt(lIndCl));
                    dMinDist(iI,2) = min(dNavInt(lIndCl));
                end
                dMinDist = abs(dMinDist - repmat(dGatePos, [1 2]));
                dMaxDist = abs([[dMinDist(1,1); diff(dGatePos(1:iNPhasesOld))], [diff(dGatePos(1:iNPhasesOld)); dMinDist(iNPhasesOld,2)]]);
                dTolerance = dMinDist + (dTolerance-1) .* dMaxDist;                
        end
    end
    % get respiratory phase mask
    dRespPhases = zeros(length(dNavInt),1);
    for iPh = iNPhasesLoop
        if(strcmp(sGatingMode,'kmeans'))
            iTol = [1 2];
        else
            iTol = [1 1];
        end
        lMask = dNavInt <= dGatePos(iPh) + dTolerance(iPh,iTol(1)) & dNavInt > dGatePos(iPh) - dTolerance(iPh,iTol(2));
        dRespPhases(lMask) = iPh;
    end
else % ignore respiratory gating
    dNavInt = ones(length(iLCReduced), 1)';
    iNPhases = 1; % needed in for loop
    iNPhasesLoop = 1;
    dGatePos = 1;
    dTolerance = 1;
    dRespPhases = ones(length(iLCReduced),1);
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Show the navigator signal, phase centers and tolerances
if lDebug
    dColors = [0 1 0;0 0 1;1 0 0;0 0.965517241379310 1;1 0.482758620689655 0.862068965517241;1 0.862068965517241 0.482758620689655;0 0.551724137931035 1;0 0.413793103448276 0.0344827586206897;0.620689655172414 0.275862068965517 0.206896551724138;0.241379310344828 0 0.413793103448276;1 0.931034482758621 1];
    hfig_respSignal = figure('name',sprintf('Gating - Ph%02d, Tol%03.0f, t%s, Esp%s, InExh%s, %s, %s', iNPhases, 100*dToleranceOld, dTime, sGatingMode, sGatingGather));
    plot(dNavInt, 'Color', 'k');
    hold on;
    dXMax = max(get(gca, 'XLim'));
    for iI = 1:length(dGatePos)
        if(strcmp(sGatingMode,'kmeans'))
            line([0 dXMax], [dGatePos(iI) dGatePos(iI)], 'Color', dColors(iI,:), 'LineWidth', 1.5); 
            line([0 dXMax], [dGatePos(iI) + dTolerance(iI,1) dGatePos(iI) + dTolerance(iI,1)], 'Color', dColors(iI,:), 'LineStyle', '--');
            line([0 dXMax], [dGatePos(iI) - dTolerance(iI,2) dGatePos(iI) - dTolerance(iI,2)], 'Color', dColors(iI,:), 'LineStyle', '--');
            hold on;               
            lIndCl = lInd;
            lIndCl(lInd) = iClusterEx == find(dGatePos(iI) == dGatePosEx);
            if(iI < iNPhasesOld)
                lIndVS = (iClusterEx == find(dGatePos(iI) == dGatePosEx) & iClusterEx == find(dGatePos(iI+1) == dGatePosEx)).'; % view sharing samples      
            else
                lIndVS = false(1,length(lIndCl));
            end 
            lIndCl(lIndCl & lIndVS) = false;

            plot(find(lIndCl),dNavInt(lIndCl),'Color', dColors(iI,:), 'Marker', 'x', 'LineStyle', 'none');
            hold on;
            plot(find(lIndVS),dNavInt(lIndVS),'Color', (dColors(iI,:)+dColors(iI+1,:))/2, 'Marker', 'x', 'LineStyle', 'none');
            hold on;
        else % percentile
            line([0 dXMax], [dGatePos(iI) dGatePos(iI)], 'Color', dColors(iI,:), 'LineWidth', 1.5); 
            line([0 dXMax], [dGatePos(iI) + dTolerance(iI) dGatePos(iI) + dTolerance(iI)], 'Color', dColors(iI,:), 'LineStyle', '--');
            line([0 dXMax], [dGatePos(iI) - dTolerance(iI) dGatePos(iI) - dTolerance(iI)], 'Color', dColors(iI,:), 'LineStyle', '--');
            hold on;
            lIndCl = dNavInt <= dGatePos(iI) + dTolerance(iI) & dNavInt > dGatePos(iI) - dTolerance(iI); 
            if(iI == 1)
                lIndVS = (dNavInt <= dGatePos(iI+1) + dTolerance(iI+1) & dNavInt > dGatePos(iI) - dTolerance(iI));
            elseif(iI == iNPhasesOld)
                lIndVS = (dNavInt > dGatePos(iI-1) - dTolerance(iI-1) & dNavInt <= dGatePos(iI) + dTolerance(iI));
            else
                lIndVS = (dNavInt <= dGatePos(iI+1) + dTolerance(iI+1) & dNavInt > dGatePos(iI) - dTolerance(iI)) | ...
                    (dNavInt > dGatePos(iI-1) - dTolerance(iI-1) & dNavInt <= dGatePos(iI) + dTolerance(iI)); 
            end

            lIndCl(lIndCl & lIndVS) = false;

            plot(find(lIndCl),dNavInt(lIndCl),'Color', dColors(iI,:), 'Marker', 'x', 'LineStyle', 'none');
            hold on;
            plot(find(lIndVS),dNavInt(lIndVS),'Color', (dColors(iI,:)+dColors(iI+1,:))/2, 'Marker', 'x', 'LineStyle', 'none');
            hold on;
        end
    end
    plot(dNavInt, 'Color', 'k');
end
% -------------------------------------------------------------------------


% -------------------------------------------------------------------------
% Populate the k-Space
fprintf('*** k-Space statistics: ***\n');
switch sGatingGather
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    case 'closest'
        if(lflagDisp), hF = waitbar(0, 'Populating k-Space'); end;
        fid = fopen([sPath, filesep, sName, '.dat'], 'r');
        
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        % Pre-allocate k-Space memory
        dKSpace = complex(zeros(iNLines, iBaseRes.*2, iNPartitions, length(iNPhasesLoop), iNChannels,sPrecision), ...
                      zeros(iNLines, iBaseRes.*2, iNPartitions, length(iNPhasesLoop), iNChannels,sPrecision));
        dKSpaceImg = zeros(iNLines, iNPartitions, length(iNPhasesLoop), 'uint8');       
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -        
        
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        % RESP PHASES loop
        for iPh = iNPhasesLoop
            lRespMask = dRespPhases == iPh;
            dWeight = abs(dNavInt - dGatePos(iPh))';

            % -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
            % LINES loop
            for iL = 1:iNLines
                lLineMask = iLine == (iL - 1);

                % -    -    -    -    -    -    -    -    -    -    -  
                % PARTITIONS loop
                for iP = 1:iNPartitions
                    lMask = lRespMask & lLineMask & iPartition == (iP - 1);
                    dThisDist = dWeight(lMask);
                    if isempty(dThisDist), continue, end;

                    [dVal, iInd] = min(dThisDist);
                    if(strcmp(sGatingMode,'kmeans'))
                        if(dNavInt(iInd) >= dGatePos(iPh))
                            iTol = 1;
                        else
                            iTol = 2;
                        end;
                    else
                        iTol = 1;
                    end
                    if(dVal > dTolerance(iPh,iTol))
                        continue; % just for backup
                    end

                    iThisSP = iSP(lMask);
                    fseek(fid, double(iThisSP(iInd)), 'bof');
                    dData = fread(fid, iADCLength.*iNChannels, 'float');
                    dData = reshape(dData, [iADCLength, iNChannels]);
                    dData = dData(128/4 + 1:end,:);
                    dData = complex(dData(1:2:end,:), dData(2:2:end,:));
                    dKSpace(iL, :, iP, iPh==iNPhasesLoop, :) = permute(dData, [3 1 4 2]);
                    dKSpaceImg(iL, iP, iPh==iNPhasesLoop) = 1;
                end
                % end of PARTITIONS loop
                % -    -    -    -    -    -    -    -    -    -    -  
                if(lflagDisp), waitbar((double(iPh - 1).*double(iNLines) + double(iL))./(double(iNPhases).*double(iNLines)), hF); end;

            end
            % end of LINES loop
            % -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
            fprintf(1, '  Phase %u acceleration is %1.3f\n', iPh, double(numel(dKSpaceImg(:,:,iPh==iNPhasesLoop)))./double(sum(sum(dKSpaceImg(:,:,iPh==iNPhasesLoop) > 0))));

        end
        % end of PHASES loop
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -             
        fclose(fid);
        if(lflagDisp), close(hF); end;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -   
    case 'collect'
        if(lflagDisp), hF = waitbar(0, 'Populating k-Space'); end;
        fid = fopen([sPath, filesep, sName, '.dat'], 'r');
                
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        % find maximal amount of samples per gate
        dKSpaceImg = cell(iNLines, iNPartitions, iNPhases);
        
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        % PHASES loop
        for iPh = 1:iNPhases
            if(strcmp(sGatingMode,'kmeans'))
                iTol = [1 2];
            else
                iTol = [1 1];
            end
            dWeight = dNavInt <= dGatePos(iPh) + dTolerance(iPh,iTol(1)) & dNavInt > dGatePos(iPh) - dTolerance(iPh,iTol(2));
            dDist = abs(dNavInt - dGatePos(iPh))';
            
            % -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
            % LINES loop
            for iL = 1:iNLines
                lLineMask = iLine == (iL - 1);
                
                % -    -    -    -    -    -    -    -    -    -    -    -
                % PARTITIONS loop
                for iP = 1:iNPartitions
                    lPartitionMask = iPartition == (iP - 1);
                    iInd = lLineMask & lPartitionMask & dWeight.';
                                        
                    if isempty(iInd), continue, end
                    [~,dThisDist] = sort(dDist(iInd));
                    iInd = find(iInd);
                    iInd = iInd(dThisDist);
                    
                    dKSpaceImg{iL, iP, iPh} = iInd;
                end
                if(lflagDisp), waitbar((double(iPh - 1).*double(iNLines) + double(iL))./(double(iNPhases).*double(iNLines)), hF); end;
            end            
            fprintf(1, '  Phase %u acceleration is %1.3f\n', iPh, (iNPartitions*iNLines)./(double(nnz(cellfun(@isempty,dKSpaceImg(:,:,iPh))))));
        end
        if(lflagDisp), close(hF); end;
        
        % decide for the most occuring gate filling situations
        % if any gate has more samples (k-space center), take the iNGate 
        % nearest samples
        iSamples = cellfun(@length, dKSpaceImg);
        n_elements = hist(iSamples(:),1:max(iSamples(:)));
        iNGate = 1;
        while(sum(n_elements(1:iNGate)) < 0.99*sum(n_elements(:))) 
            iNGate = iNGate + 1;
        end        
        
        lMask = iSamples <= iNGate & iSamples > 0;
        lMask(iEchoLine-round(iNavPERes/2)+1:iEchoLine+round(iNavPERes/2), iEchoPartition-round(iNav3DRes/2)+1:iEchoPartition+round(iNav3DRes/2), :) = true; % keep navigator lines
        lMask = find(~lMask);
        for iMask = 1:length(lMask)
            [cLine,cPar,cPha] = ind2sub(size(dKSpaceImg),lMask(iMask));
            dKSpaceImg{cLine,cPar,cPha} = [];
        end
        
        % Pre-allocate k-Space memory
        dKSpace = complex(zeros(iNLines, iBaseRes.*2, iNPartitions, iNPhases, iNGate, iNChannels, sPrecision), ...
                          zeros(iNLines, iBaseRes.*2, iNPartitions, iNPhases, iNGate, iNChannels, sPrecision));
                      
        if(lflagDisp), hF = waitbar(0, 'Populating k-Space'); end;
        
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
        % PHASES loop
        for iPh = 1:iNPhases
            
            % -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
            % LINES loop
            for iL = 1:iNLines
                
                % -    -    -    -    -    -    -    -    -    -    -    -
                % PARTITIONS loop
                for iP = 1:iNPartitions  
                    if(isempty(dKSpaceImg{iL,iP,iPh})), continue; end;
                    
                    % -     -     -     -     -     -     -     -     -    
                    % GATING loop
                    for iG = 1:length(dKSpaceImg{iL,iP,iPh})
                        if(iG > iNGate), break, end;
                        iThisSP = iSP(dKSpaceImg{iL,iP,iPh}(iG));
                        fseek(fid, double(iThisSP), 'bof');
                        dData = fread(fid, iADCLength.*iNChannels, 'float');
                        dData = reshape(dData, [iADCLength, iNChannels]);
                        dData = dData(128/4 + 1:end,:);
                        dData = complex(dData(1:2:end,:), dData(2:2:end,:));
                        dKSpace(iL, :, iP, iPh, iG, :) = permute(dData, [3 1 4 2]);
                    end
                    % end of GATING loop
                    % -     -     -     -     -     -     -     -     -                       
                end
                % end of PARTITIONS loop
                % -    -    -    -    -    -    -    -    -    -    -    -
                if(lflagDisp), waitbar((double(iPh - 1).*double(iNLines) + double(iL))./(double(iNPhases).*double(iNLines)), hF); end;
            end
            % end of LINES loop
            % -   -   -   -   -   -   -   -   -   -   -   -   -   -   -   -
        end
        % end of PHASES loop
        % -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -       
        fclose(fid);
        if(lflagDisp), close(hF); end;
        
        
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % No k-Space population, show simulations of the k-Space sampling
    % density
    case 'none'
        dKSpaceImg = zeros(iNLines, iNPartitions, iNPhases, 'uint8');
        for iPh = 1:iNPhases
            lMask = abs(dNavInt - dGatePos(iPh)) <= dTolerance;
            
            iThisLine = iLine(lMask);
            iThisPartition = iPartition(lMask);
            
            for iI = 1:length(iThisLine)
                dKSpaceImg(iThisLine(iI) + 1, iThisPartition(iI) + 1, iPh) = ...
                dKSpaceImg(iThisLine(iI) + 1, iThisPartition(iI) + 1, iPh) + 1;
            end
            fprintf(1, '  Phase %u acceleration is %1.3f\n', iPh, double(numel(dKSpaceImg(:,:,iPh)))./double(sum(sum(dKSpaceImg(:,:,iPh) > 0))));
        end
        dKSpace = dKSpaceImg;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -       
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% show subsampling masks
if lDebug
    dColMap = OptimalColor(256 + 20);
    dColMap = dColMap(21:end, :);
    dColMap(1,:) = [0, 0, 0];
    
    for iPh = 1:iNPhases
        figure, imshow(dKSpaceImg(:,:,iPh)', 'InitialMagnification', 400);
        set(gcf, 'Name', sprintf('k-Space density for phase %u', iPh));
        colormap(dColMap);
        set(gca, 'CLim', [0 255]);
    end
end
% -------------------------------------------------------------------------