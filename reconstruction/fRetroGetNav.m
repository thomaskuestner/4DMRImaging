function dNav = fRetroGetNav(sFilename, lDebug)
%FFLASHRETROGETNAV Extracts navigator data from CS_FLASH_retro sequence
%   DNAV = FFLASHRETROGETNAV(SFILENAME) Extracts navigator data DNAV from
%   the siemens meas-data file SFILENAME. Must have been parsed with
%   FMeasCreateLUT before.
%
% See also: FMEASCREATELUT
%
%   Copyright 2014-2016 Christian Wuerslin, University of Tuebingen, Germany
%   christian.wuerslin@med.uni-tuebingen.de
%   and Thomas Kuestner, University of Tuebingen, Germany
%   thomas.kuestner@med.uni-tuebingen.de

if(nargin < 1)
    lDebug = false;
end
[sPath, sName, sExt] = fileparts(sFilename);

% -------------------------------------------------------------------------
% Get some relevant data from the drecksMDH
try
    load(fullfile(sPath,[sName,'.mat']), 'SDrecksMDH');
    dBaseRes        = SDrecksMDH.Geo.MatrixSize(1);
    dTR             = SDrecksMDH.Contrast.TR./1000; % in ms
    if(isfield(SDrecksMDH.Geo,'EchoPosition'))
        dEchoLine       = SDrecksMDH.Geo.EchoPosition(1);
        dEchoPartition  = SDrecksMDH.Geo.EchoPosition(2);
    else
        dEchoLine       = SDrecksMDH.Geo.MatrixSize(2)/2;
        dEchoPartition  = SDrecksMDH.Geo.MatrixSize(3)/2;
    end            
    dNavPeriod      = SDrecksMDH.Wip.NavPeriod;
    dNavPERes       = SDrecksMDH.Wip.NavRes(1);
catch
    % Load the loop counters
    load([sPath, filesep, sName, '.mat']);
    dBaseRes        = 256;
    dEchoLine       = 128;
    dEchoPartition  = 36;
    dNavPeriod      = 200;
    dNavPERes       = 1;
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Get ktSpace
[dKtSpace, iLC] = fMeasRead(sFilename, 'Set', [1 7], 'Eco', 0, 'Lin', dEchoLine, 'Par', dEchoPartition);
if(length(unique(iLC(:,10))) > 1) % in the case of switched on surrogate monitoring
    lMask = true(size(iLC, 1), 1);
    lMask = lMask & (iLC(:,10) == 1 | iLC(:,10) == 3 | iLC(:,10) == 5 | iLC(:,10) == 7);
    iLC = iLC(lMask,:);
    dKtSpace = dKtSpace(lMask,:);
end    
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Determine data size
iNChannels = double(iLC(1, 2));
iNSamples = size(dKtSpace, 2);
iNMeasurements = size(dKtSpace, 1)./(iNChannels);
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Build up the kSpaces of all channels
dKtSpace = reshape(dKtSpace, [iNChannels, iNMeasurements, iNSamples]);
dKtSpace = permute(dKtSpace, [3, 2, 1]); % -> RO x t x CH
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Reconstruct the 1-D projections for all measuremetns and all channels
dImg = fftshift(ifft(ifftshift(dKtSpace(:,:,:))));
dImg = flipdim(dImg, 1); % Invert the RO direction: 1-N -> H-F
dImg = dImg(iNSamples/4:iNSamples.*3/4 - 1, :, :, :); % RO x t x CH

% Get x range of respiratory motion and sort out Channels with no relevant
% information in target area
dFreq = fft(dImg, [], 2);
dPower = dFreq.*conj(dFreq);
dIMGres = 1./(double(dNavPeriod)./1000.*double(iNMeasurements)); % The frequency resolution of dIMG in Hz
dPower = squeeze(sum(dPower(:, round(1./(5.*dIMGres)):round(1./(3.*dIMGres)), :), 2)); % RO x CH
dPowerInChan = dPower; % RO x CH
dPowerAcrossChannels = sum(dPowerInChan, 2); % RO x 1
dPowerAcrossChannels(length(dPowerAcrossChannels).*3/4:end) = 0; % Prevent detection of regions in the abdomen
dPowerAcrossChannels = conv(dPowerAcrossChannels, fGaussianLP(20), 'same');
[dMax, dX] = max(dPowerAcrossChannels);
dMax = max(dPowerInChan(dX - 20: dX + 20, :));
lGoodChannels = dMax > 0.4.*max(dMax);

dSOSImg = sqrt(sum(dImg(:,:,lGoodChannels).*conj(dImg(:,:,lGoodChannels)), 3));
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Get the navigator signal
iDisplacement = 80; % [mm] diaphragm displacement (max. +/- 80mm)
dDisplacement = round(iDisplacement/(SDrecksMDH.Geo.FOV(1)/SDrecksMDH.Geo.MatrixSize(1))); % [px]

dRefImg = zeros(size(dSOSImg));
dRefImg(:,end) = dSOSImg(:,end);
idx = 1:size(dRefImg,2);
dNav = zeros(1,size(dRefImg,2));

hF = waitbar(0, 'Getting Navigator');
for i=size(dRefImg,2)-1:-1:2
    dRMSImg = zeros(size(dRefImg,1),2*dDisplacement+1);
    tmp = dRefImg(:,idx(i+1:end));
    for iD = -dDisplacement:dDisplacement
        dRMSImg(:,dDisplacement-iD+1) = sum((tmp - repmat(circshift(dSOSImg(:,idx(i)), iD),[1 size(tmp,2)])).^2,2); % RMS
        waitbar((abs(i-size(dRefImg,2)-2)*(2*dDisplacement+1) + (iD + dDisplacement + 1))/((size(dRefImg,2)-2)*(2*dDisplacement+1)),hF);
    end
    [~, dNav(i)] = min(sum(dRMSImg(dX-round(dDisplacement/2):dX+round(dDisplacement/2),:))); % RMS
    dNav(i) = dDisplacement + 1 - dNav(i);
    dRefImg(:,idx(i)) = circshift(dSOSImg(:,idx(i)), dNav(i));
    waitbar(abs(i-size(dRefImg,2)-2)/(size(dRefImg,2)-2),hF);
end
close(hF);
dNav = -dNav;
dNav = conv(dNav, fGaussianLP(5), 'same');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Visualize the result
if(lDebug)
    figure('Name','Navi'), imagesc(dSOSImg(:, 3:end));
    colormap gray;
    hold all;
    plot(dNav(3:end)-mean(dNav)+dX);
end
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Interpolate the navigator data to TR intervals
load([sPath, filesep, sName, '.mat'], 'iLC');
iLC = iLC(1:iNChannels:end, :);
iLC = iLC(iLC(:,7) == 0, :); % first echo
iLC = iLC(iLC(:, 1) == 2*dBaseRes, :); % only "real" readouts

dNav = - dNav;
dNav = dNav - min(dNav);

iNavInd = find((iLC(:,3) == dEchoLine) & (iLC(:,6) == dEchoPartition));
if(length(iNavInd) ~= length(dNav)) % happens at the end when navi not fully acquired
    iNavInd = iNavInd(1:length(dNav));
end
dNavInt = interp1(iNavInd, dNav, 1:length(iLC), 'pchip');
dNavInt_ms = interp1(iNavInd.*dTR, dNav, 1:length(iLC).*dTR, 'pchip');
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% Append the navigator data to the LUT file
save([sPath, filesep, sName, '.mat'], 'dNavInt', '-append');
save([sPath, filesep, sName, '.mat'], 'dNavInt_ms', '-append');
save([sPath, filesep, sName, '.mat'], 'dSOSImg', '-append');
% -------------------------------------------------------------------------