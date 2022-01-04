function [h, slope] = TruthEval(x, y, shape)
%TRUTHEVAL Evaluates the true height of the model at each z = f(x(i),y(j))

%locals
nx = length(x);
ny = length(y);
h = zeros(ny,nx);
slope = 0;

% Persistent data for crater terrain model
persistent gsdPersist terrainPersist;
if strcmp(shape, 'RocksAndCraters')
    if isempty(gsdPersist)
        load('RockAndCraterDEM.mat','gsd','rockAndCraterDEM');
        gsdPersist = gsd;
        terrainPersist = rockAndCraterDEM;
    end
end

% run loop
for ii = 1:nx
    for jj = 1:ny
        
        switch shape
            %         case 'DoubleRamp'
            %
            %             if(r(ii) < 0)
            %                 h(ii) = r(ii);
            %                 slope(ii) = 1;
            %             else
            %                 h(ii) = 0.5*r(ii);
            %                 slope(ii) = 0.5;
            %             end
            %         case 'Parabola'
            %
            %             h(ii) = 0.05*r(ii)^2;
            %             slope(ii) = 0.10*r(ii);
            
            case 'GentleParabola'
                a = 1E-4;
                h(jj,ii) = a*(x(ii)^2 + y(jj)^2);
                %             slope(ii) = 2*a*r(ii);
                
                %         case 'Rock'
                %             if(abs(r(ii)) > 0.5)
                %                 h(ii) = 0;
                %                 slope(ii) = 0;
                %             else
                %                 h(ii) = cos(pi*r(ii));
                %                 slope(ii) = -pi*sin(pi*r(ii));
                %             end
                %
                %         case 'Line'
                %             m = .5;
                %             b = 2;
                %             h(ii) = m*r(ii) + b;
                %             slope(ii) = m;
                
            case 'FlatBottomParabola'
                a = 1E-3;
                r = 200; % radius of flat bottom
                
                rad = x(ii)^2 + y(jj)^2;
                if(sqrt(rad) < r)
                    h(jj,ii) = 0;
                else
                    h(jj,ii) = a*(rad - r^2);
                end

            case 'RocksAndCraters'

%                 % define the x coordinates of the terrain
%                 [NterrainY, NterrainX] = size(terrainPersist);
%                 spanX = NterrainX*gsdPersist;
%                 spanY = NterrainY*gsdPersist;
%                 xterrain = linspace(-spanX/2,spanX/2,NterrainX);
%                 yterrain = linspace(-spanY/2,spanY/2,NterrainY);
% 
%                 % Interpolate terrain at query points
%                 xeval = repmat(x,ny,1);
%                 yeval = repmat(y',1,nx);
%                 h = interp2(xterrain,yterrain,terrainPersist,xeval,yeval);
% 
%                 % we are able to do this all in one call so break
%                 return

                  % Convert x and y position to indicies in terrain
                  [NterrainY, NterrainX] = size(terrainPersist);
                  idxX = round(x(ii)/gsdPersist + NterrainX/2);
                  idxY = round(y(jj)/gsdPersist + NterrainY/2);

                  % Constrain idxes
                  idxX = max(1,idxX);
                  idxX = min(idxX,NterrainX);
                  idxY = max(1,idxY);
                  idxY = min(idxY,NterrainY);

                  h(jj,ii) = terrainPersist(idxY,idxX);

            otherwise
                disp("Error: Invalid Truth Shape!")
                return
        end
        
    end
    
end

end

