function dImg = fRetroRecon4D_main(sFilename, iNPhases, dTolerance, dTime, sGating, sGatingGather, lDebug)
%
% main function for CS Retro reconstruction
%
%   Copyright 2014-2016 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   and Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

% inputs
% sFilename         string of TWIX measurement file
% iNPhases          number of gates
% dTolerance        view sharing blending factor(s) b of 
%                   respiratory gates: 1 (non-overlapping) <= b <= 2 (completely overlapping) 
% dTime             scalar/vector/array of crop scan time(s) (0 := no cropping | scan time [s])
% sGating           string to select gating procedure
%                   'percentile'    select gate centroids according to 5% and 95% percentile
%                   'kmeans'        perform k-means clustering to select gate centroids
% sGatingGather     string of sample gathering mode: 
%                   'closest'   take closest sample to gate centroid
%                   'collect'   gather all samples and let the
%                               reconstruction decide for the most
%                               favourable one; not working for cardiac
%                   'none'      no k-Space population, show simulations of
%                               the k-Space sampling density
% 
% output
% dImg              reconstructed image

% -------------------------------------------------------------------------
% parse inputs
if nargin < 7, lDebug = false; end
if nargin < 6, sGatingGather = 'closest'; end
if nargin < 5, sGating = 'percentile'; end
if nargin < 4, dTime = 0; end
if nargin < 3, dTolerance = 1; end
if nargin < 2, iNPhases = 4; end
if(nargin < 1 || isempty(sFilename))
    [sFile, sPath] = uigetfile({'*.dat', 'ADC measdata (*.dat)'} ,'Select ADC measdata file', 'MultiSelect', 'off', pwd);
    if(~isempty(sFile))
        [~, sName] = fileparts(sFile);
    else
        return;
    end
    sFilename = [sPath,sFile];
end

if(~exist('CS_reconstruction','file'))
    error('CS_reconstruction code not found! Please download it from "https://github.com/thomaskuestner/CS_LAB" or put it on Matlab path');
end
currpath = fileparts(mfilename('fullpath'));
addpath(genpath([currpath,filesep,'utils']));

% calculation precision
sPrecision = 'single';
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% load raw measurment file
[sPath, sName] = fileparts(sFilename);
[~, ~, ~, para.drecksMDH, para.iLCPositions] = fMeas_main([sPath, filesep, sName, '.dat'], false);
load([sPath, filesep, sName, '.mat'], 'iLC');
para.iLC = iLC;
clear 'iLC';
para.drecksMDH.Wip.NSamples = 0;
para.drecksMDH.Wip.Phases = iNPhases;
para.measPara.sequenceName = 'CS_Retro';
para.measPara.dim = [para.drecksMDH.Geo.MatrixSize(2), 2*para.drecksMDH.Geo.MatrixSize(1), para.drecksMDH.Geo.MatrixSize(3), iNPhases, double(para.iLC(1,2))];
para.measPara.LCall = ones(1,4);
dNavPERes = para.drecksMDH.Wip.NavRes(1);
dNavSLRes = para.drecksMDH.Wip.NavRes(2);
if(para.measPara.dim(3) == 1)
    para.measPara.dimension = '2D';
else
    if(para.measPara.dim(4) == 1)
        para.measPara.dimension = '3D';
    else
        para.measPara.dimension = '4D';
    end    
end
iNEchos = length(unique(para.iLC(:,7)));
% -------------------------------------------------------------------------
        
% -------------------------------------------------------------------------
% prepare output file
if(isscalar(dTime))
    dTShow = num2str(dTime,'%03d');
else
    dTShow = dTime.';
    dTShow = strrep(num2str(dTShow(:).','%03d'),' ','');
end     
sSaveFilename = [sPath, filesep, sName, sprintf('_Ph%1d_Tol%03.0f_t%s_%s_%s', iNPhases, 100*dTolerance, dTShow, sGating, sGatingGather)];                                   
% -------------------------------------------------------------------------
        
for iEcho = 1:iNEchos
    if(iNEchos > 1)
        sSaveFilename = [sSaveFilename,'_',num2str(iEcho,'E%02d')];
    end
    [~, sFilenameShow, ~] = fileparts(sSaveFilename);

    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % extract self-navigation signal
    if(~any(arrayfun(@(x) strcmp(x.name,'dNavInt'), whos('-file', [sPath, filesep, sName, '.mat']))))
        if(dNavPERes == 1 && dNavSLRes == 1) % 1D navigator signal
            dNav = fRetroGetNav(sFilename,lDebug);
        elseif(dNavPERes > 1 && dNavSLRes == 1) % 2D navigator signal
            dNav = fRetroGetNav2D(sFilename,lDebug);
        else % 3D navigator signal
            dNav = fRetroGetNav3D(sFilename,lDebug);
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % perform retrospective gating
    if(~exist([sSaveFilename,'_kSpace.mat'],'file'))
        kSpace = fRetroGating(sFilename, iNPhases, dTolerance, dTime, sGating, sGatingGather, sPrecision, iEcho, false);
        save([sSaveFilename,'_kSpace.mat'], 'kSpace', '-v7.3');
    else
        fprintf('Loading k-space: %s\n',sFilenameShow);
        load([sSaveFilename,'_kSpace.mat']);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Compressed Sensing reconstruction                             
    if(strcmp(sGatingGather,'collect'))
        kSpace = mat2cell(kSpace, size(kSpace, 1), size(kSpace, 2), size(kSpace, 3), size(kSpace, 4), size(kSpace,5), ones(1, size(kSpace, 6)));
        kSpace = shiftdim(kSpace, 4);
    else
        kSpace = mat2cell(kSpace(:,:,:,:,:), size(kSpace, 1), size(kSpace, 2), size(kSpace, 3), size(kSpace, 4), ones(1, size(kSpace, 5)));
        kSpace = shiftdim(kSpace, 3);
    end

    nSize = size(kSpace{1});
    para = fPrepareRecon(para,nSize,sPath,sGatingGather);                                    
    [imgAll,imgCha] = CS_reconstruction(kSpace,para);
    if(strcmp(sGatingGather,'collect') || ndims(imgAll) == 4)
        dImg = imgAll;
    else
        dImg = zeros([size(imgAll,1), size(imgAll,2), nSize(3), nSize(4)]);
        for i=1:iNPhases
            dImg(:,:,:,i) = imgAll(:,:,(i-1)*nSize(3)+1:i*nSize(3));
        end
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    
    
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % save results
    save([sSaveFilename,'_channel.mat'], 'imgCha', '-v7.3');
    save([sSaveFilename,'_recon.mat'], 'dImg', '-v7.3');                                
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end
end


function para = fPrepareRecon(para,nSize,sPath,sGatingGather)
% set CS reconstruction parameters
    nPha = nSize(1);
    nFreq = nSize(2);
    nZ = nSize(3);
    nTime = nSize(4);
    if(length(nSize) > 4)
        nSmp = nSize(5);
    else
        nSmp = 1;
    end
    
    % additional length of gating dimension
    % drecksMDH v2.0
    para.drecksMDH.Wip.Phases = nTime;
    if(strcmp(sGatingGather,'collect'))
        para.drecksMDH.Wip.NSamples = nSmp;
    end                   
    
    % CS recon parameter
    if(strcmp(para.measPara.dimension,'3D'))
        para.trafo.trafodim = [1 2 1 0; 1 0 1 0]; % 3D
    elseif(strcmp(para.measPara.dimension,'4D'))
        para.trafo.trafodim = [1 2 1 0; 0 0 0 1]; % 4D
    end
    para.cstype = 'FOCUSS';
    para.transformation = 'pca';
    para.trafo.scrambledim = [1 1 1 0];
    para.lambda = 0.01;
    para.lambdaCalib = 0; % GRAPPA off
    para.lambdaTV = 0;
    para.lambdaESPReSSo = 0;
    para.lambdaMC = 0;
    para.sScaling = 'self';
    para.flagOversampling = logical([1 0 0]); % oversampling correction
    para.prop.flagPlot = false; % Imagine off
    para.prop.flagOptim = false; % no FFT optimization
    para.savePath = sPath;
    
end