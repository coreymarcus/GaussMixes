clear
close all
clc

% Setup terrain parameters
rng(3)
width = 80;
height = 80;
xmin = -40;
ymin = -40;
xmax = xmin + width;
ymax = ymin + height;
gsd = .1;
xsample = xmin:gsd:xmax;
ysample = ymin:gsd:ymax;
Nx = length(xsample);
Ny = length(ysample);
groundA = 1E-4;
groundFlatRadius = 200;
rockRadiusMean = 2.5;
rockRadiusStdDev = 0.25;
rockAbundance = 0.08;
craterDepthToWidthMean = 0.20;
craterDepthToWidthStdDev = 0.05;
craterRadiusMean = 10;
craterRadiusStdDev = 0.5;
craterAbundance = 0.08;

% Use ALHAT approximate data?
useALHAT = true;

if useALHAT
    % Read rocks
    rocks = xlsread('ApproxALHATRocks.xlsx');
    
    % Extract data
    numRocks = size(rocks,1);
    rockX = rocks(:,1);
    rockY = rocks(:,2);
    rockRadius = rocks(:,3);
    
    % read craters
    craters = xlsread('ApproxALHATCraters.xlsx');
    
    % Extract data
    numCraters = size(craters,1);
    craterX = craters(:,1);
    craterY = craters(:,2);
    craterRadius = craters(:,3);
    craterDepth = craters(:,4);
    
    % Calculate surface normal parameters
    rise = 1.0;
    run = 60;
    nhat = [-rise 0 run]';
    nhat = nhat/norm(nhat);
    offset = 0;
    
else
    % Find the number of rocks and craters
    terrainArea = width*height;
    numRocks = round(rockAbundance*terrainArea/(pi*rockRadiusMean^2));
    numCraters = round(craterAbundance*terrainArea/(pi*craterRadiusMean^2));
    
    % Sample rock and crater paremeters
    rockX = width*rand([numRocks,1]);
    rockY = height*rand([numRocks,1]);
    rockRadius = abs(mvnrnd(rockRadiusMean,rockRadiusStdDev^2,numRocks));
    
    craterX = width*rand([numCraters,1]);
    craterY = height*rand([numCraters,1]);
    craterRadius = abs(mvnrnd(craterRadiusMean,craterRadiusStdDev^2,numCraters));
    craterDepth = abs(mvnrnd(craterDepthToWidthMean,craterDepthToWidthStdDev^2,numCraters).*2.0.*craterRadius);
    
    % Use a flat ground plane
    nhat = [0 0 1]';
    offset = 0;
end

% Generate rock plots
rockDEM = zeros(Ny,Nx);
for ii = 1:numRocks
    x = round(rockX(ii));
    y = round(rockY(ii));
    rad = rockRadius(ii);
    cellRad = floor(rad/gsd);
    rockSample = gsd*(-cellRad:cellRad);
    Nrock = length(rockSample);
    rockDEMiter = zeros(Nrock);
    
    % Create rock height
    for jj = 1:Nrock
        for kk = 1:Nrock
            dist = sqrt(rockSample(jj)^2 + rockSample(kk)^2);
            if(dist >= rad)
                continue;
            end
            rockDEMiter(jj,kk) = sqrt(rad^2 - rockSample(jj)^2 - rockSample(kk)^2);
        end
    end
    
    % Add rock to large DEM
    rockSampleX = round((rockSample + x - xmin)/gsd);
    rockSampleY = round((rockSample + y - ymin)/gsd);
    targX = rockSampleX >= 0 & rockSampleX < Nx;
    targY = rockSampleY >= 0 & rockSampleY < Ny;
    rockDEMiter = rockDEMiter(targY, targX);
    rockDEM(rockSampleY(targY) + 1,rockSampleX(targX) + 1) = ...
        max(rockDEM(rockSampleY(targY) + 1,rockSampleX(targX) + 1), ...
        rockDEMiter);
    
end

% Generate crater plots
craterDEM = zeros(Ny,Nx);
for ii = 1:numCraters
    x = round(craterX(ii));
    y = round(craterY(ii));
    rad = craterRadius(ii);
    depth = craterDepth(ii);
    cellRad = floor(rad/gsd);
    craterSample = gsd*(-cellRad:cellRad);
    Ncrater = length(craterSample);
    craterDEMiter = zeros(Ncrater);
    
    % Find z location of crater sphere center
    R = (depth^2 + 0.25*(2*rad)^2)/(2*depth);
    z = R - depth;
    
    % Create crater height
    for jj = 1:Ncrater
        for kk = 1:Ncrater
            dist = sqrt(craterSample(jj)^2 + craterSample(kk)^2);
            if(dist >= R)
                continue;
            end
            craterDEMiter(jj,kk) = z - sqrt(R^2 - craterSample(jj)^2 - craterSample(kk)^2);
        end
    end
    
    % Add crater to large DEM
    craterSampleX = round((craterSample + x - xmin)/gsd);
    craterSampleY = round((craterSample + y - ymin)/gsd);
    targX = craterSampleX >= 0 & craterSampleX < Nx;
    targY = craterSampleY >= 0 & craterSampleY < Ny;
    craterDEMiter = craterDEMiter(targY, targX);
    craterDEM(craterSampleY(targY) + 1,craterSampleX(targX) + 1) = ...
        min(craterDEM(craterSampleY(targY) + 1,craterSampleX(targX) + 1), ...
        craterDEMiter);
    
end

% Generate ground plane DEM
groundDEM = zeros(Ny,Nx);
for ii = 1:Nx
    for jj = 1:Ny
        groundDEM(jj,ii) = (offset - ii*nhat(1)*gsd - jj*nhat(2)*gsd)/nhat(3);
    end
end

% Create final DEM by adding rocks and craters together
rockAndCraterDEM = craterDEM + rockDEM + groundDEM;

% Save DEM to disk
save('RockAndCraterDEM','rockAndCraterDEM',"gsd");

% Plot
figure
scatter(rockX,rockY)
hold on
scatter(craterX,craterY)
title('Rock / Crater Locations')
legend('Rocks','Craters')

% figure
% mesh(rockDEM)
% title('Rock DEM')
%
% figure
% mesh(craterDEM)
% title('Crater DEM')

figure
mesh(xsample,ysample,rockAndCraterDEM)
title('Combined DEM')
colorbar
view(2)
