function BW = binarizeHyphae(filterResponse, thetaMax)
% Binarizes the filter response.
% Input:
% filterResponse [n x m], the maximal response magnitude.
% thetaMax [n x m], the angle (deg) yielding the maximal filter response.
% Output:
% BW [n x m], binary matrix where 1 indicates presence of hyphae.

% Normalize ridge filter response
responseMap = mat2gray(filterResponse);

% Local adaptive treshholding using Bradley's method.
BW1 = imbinarize(responseMap, 'adaptive','Sensitivity', 0.5);

% Identify pixels repressenting local maxima.
lmax = nonMaxSuppression(filterResponse, thetaMax);

% Create ridge filter response map without non maxima.
lmaxRespMap = responseMap;
lmaxRespMap(~lmax) = 0;

% Create histogram of filter response values which are also local maxima in
% the maximal response direction. Identify threshold using Otsu's method.
lmaxRespList = lmaxRespMap(:);
lmaxRespList(lmaxRespList==0) = [];
T = otsuthresh(imhist(lmaxRespList));

% Identify the seed index map by thresholding the locally maximal filter
% responses.
BW2 = (lmaxRespMap>T);
seed_indices=sub2ind(size(BW1),find(BW2));

% Perform 8-connected hysteresis 
hys=imfill(~BW1,seed_indices,8);
BW = hys & BW1;

end