function [tgUnitVec, tgStrength] = pointTangent(basicFunctions, targetPoints, opt)

pointNum = size(targetPoints,1);
tgUnitVec = double(zeros(pointNum,2));
tgStrength = double(zeros(pointNum,1));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

[onSplineFlag, segNum, segT] = isOnSpline(basicFunctions, targetPoints, opt);

for i = 1:pointNum
    if onSplineFlag(i,1)== false
        tgUnitVec(i,:) = [NaN NaN];
        tgStrength(i,1) = NaN;
    else

        b2x = basicFunctions(segNum(i),2,1);
        b2y = basicFunctions(segNum(i),2,2);
        b3x = basicFunctions(segNum(i),3,1);
        b3y = basicFunctions(segNum(i),3,2);
        b4x = basicFunctions(segNum(i),4,1);
        b4y = basicFunctions(segNum(i),4,2);
        t = segT(i);
        partialX = b2x + 2*b3x*t + 3*b4x*t*t;
        partialY = b2y + 2*b3y*t + 3*b4y*t*t;
        tgVec = [partialX partialY];
        tgStrength(i,1) = norm(tgVec);
        tgUnitVec(i,:) = tgVec./tgStrength(i,1);
        
    end
end

if plotFlag == true
    figure;
    hold on;
    xcoord = targetPoints(:,1);
    ycoord = targetPoints(:,2);
    xvec = tgUnitVec(:,1);
    yvec = tgUnitVec(:,2);
    
    xcoord = xcoord(~isnan(tgStrength));
    ycoord = ycoord(~isnan(tgStrength));
    xvec = xvec(~isnan(tgStrength));
    yvec = yvec(~isnan(tgStrength));
    quiver(xcoord, ycoord, xvec, yvec);
end