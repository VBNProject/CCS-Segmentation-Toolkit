function [xmin,xmax,ymin,ymax,rangePoints] = splineRange(basicFunctions,opt)

n = size(basicFunctions,1);
xmin = Inf;
ymin = Inf;
xmax = -Inf;
ymax = -Inf;
rangePoints = double(zeros(4,2));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

for i = 1:n
    bonex = basicFunctions(i,1,1);
    boney = basicFunctions(i,1,2);
    btwox = basicFunctions(i,2,1);
    btwoy = basicFunctions(i,2,2);
    bthrx = basicFunctions(i,3,1);
    bthry = basicFunctions(i,3,2);
    bfoux = basicFunctions(i,4,1);
    bfouy = basicFunctions(i,4,2);

    xminfun = @(t)bonex+btwox*t+bthrx*(t^2)+bfoux*(t^3);
    yminfun = @(t)boney+btwoy*t+bthry*(t^2)+bfouy*(t^3);
    xmaxfun = @(t)(-1)*(bonex+btwox*t+bthrx*(t^2)+bfoux*(t^3));
    ymaxfun = @(t)(-1)*(boney+btwoy*t+bthry*(t^2)+bfouy*(t^3));

    [mintx,minfx] = fminbnd(xminfun,0,1,optimset('TolX',1e-6));
    [minty,minfy] = fminbnd(yminfun,0,1,optimset('TolX',1e-6));
    [maxtx,maxfx] = fminbnd(xmaxfun,0,1,optimset('TolX',1e-6));
    [maxty,maxfy] = fminbnd(ymaxfun,0,1,optimset('TolX',1e-6));
    maxfx = -maxfx;
    maxfy = -maxfy;
    
    if ~isnan(mintx)
        if minfx<xmin
            xmin = minfx;
            rangePoints(1,1) = bonex+btwox*mintx+bthrx*(mintx^2)+bfoux*(mintx^3);
            rangePoints(1,2) = boney+btwoy*mintx+bthry*(mintx^2)+bfouy*(mintx^3);            
        end
    end

    if ~isnan(minty)
        if minfy<ymin
            ymin = minfy;
            rangePoints(3,1) = bonex+btwox*minty+bthrx*(minty^2)+bfoux*(minty^3);
            rangePoints(3,2) = boney+btwoy*minty+bthry*(minty^2)+bfouy*(minty^3); 
        end
    end

    if ~isnan(maxtx)
        if maxfx>xmax
            xmax = maxfx;
            rangePoints(2,1) = bonex+btwox*maxtx+bthrx*(maxtx^2)+bfoux*(maxtx^3);
            rangePoints(2,2) = boney+btwoy*maxtx+bthry*(maxtx^2)+bfouy*(maxtx^3);  
        end
    end

    if ~isnan(maxty)
        if maxfy>ymax
            ymax = maxfy;
            rangePoints(4,1) = bonex+btwox*maxty+bthrx*(maxty^2)+bfoux*(maxty^3);
            rangePoints(4,2) = boney+btwoy*maxty+bthry*(maxty^2)+bfouy*(maxty^3);
        end
    end
    
end

if plotFlag == true
    for i = 1:4
        hold on;
        xcoord = rangePoints(i,1);
        ycoord = rangePoints(i,2);
        plot(xcoord,ycoord,'b+');
    end
end

end