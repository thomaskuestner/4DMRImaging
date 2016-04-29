function dDataOut = fGaussianLP(dData, dW)
%GAUSSIANLP performs a 1D gaussian lowpass filter.
%
% HOUT = GAUSSIANLP(HIN, L, MODE) performs a gaussian lowpass of length L
% to the array in HIN. The gaussian lowpass is defined in the range
% [-floor(L/2), floor(L/2)] and the variance sigma of the Gaussinan is
% chosen such, that the definition boundaries of the Gaussian are equal to
% 3*sigma. MODE defines the method for padding the input array HIN and can
% be either 'circular', 'replicate' or 'symmetroc'. See also the
% documentation of 'padarray' for more information on the padding methods.

if nargin == 1, dW = dData; end
    
% -------------------------------------------------------------------------
% Create Gaussian kernel.
dX = (-floor(dW/2):floor(dW/2))';
dX = dX.*6./dW;
dG = exp(-0.5.*dX.^2);
dG = dG./sum(dG);
if nargin == 1
    dDataOut = dG;
    return
end
% -------------------------------------------------------------------------

dDataOut = zeros(size(dData));
for iI = 1:size(dData, 2), dDataOut(:,iI) = conv(dData(:,iI), dG, 'same'); end