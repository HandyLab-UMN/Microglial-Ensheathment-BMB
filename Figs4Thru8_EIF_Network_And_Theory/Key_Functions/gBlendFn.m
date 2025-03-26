function [gColorScheme] = gBlendFn(cMin,cMax)

blueStart = [0 0.4470 0.7410]; % 'top' color 
% redStart = [0.85 0.325 0.098]; % 'bottom' color
redStart = [1 1 1]; % 'bottom' color
middle = [1 1 1]; % middle color
numColors = 256;

% Find where the white color is located
cVec = cMin:(cMax-cMin)/(numColors-1):cMax;
middleIndex = find(cVec>=0,1);

%% Blend the red color to white
rColorSlope = middle-redStart;
t = 0:1/(middleIndex-1):1;
rBlend = zeros(length(t),3);
for ii = 1:length(t)
    rBlend(ii,:) = redStart + rColorSlope.*t(ii);
end

%% Blend the white color to blue 
bColorSlope = blueStart-middle;
t = 0:1/(length(cVec)-middleIndex-1):1;
bBlend = zeros(length(t),3);
for ii = 1:length(t)
    bBlend(ii,:) = middle + bColorSlope.*t(ii);
end

gColorScheme = [rBlend; bBlend];

end