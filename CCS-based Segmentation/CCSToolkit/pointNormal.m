function [normUnitVec, normStrength] = pointNormal(basicFunctions, targetPoints, opt)

pointNum = size(targetPoints,1);
normUnitVec = double(zeros(pointNum,2));
normStrength = double(zeros(pointNum,1));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

try
    outFlag = opt.plotFlag;
catch
    outFlag = true;
end

opt.plotFlag = false;
[tgUnitVec, tgStrength] = pointTangent(basicFunctions, targetPoints, opt);

n = size(basicFunctions,1);
xcoord = double(zeros(n,1));
ycoord = double(zeros(n,1));

B1 = double(zeros(2,1));
for i = 1:n
    B1(1,1) = basicFunctions(i,1,1);
    B1(2,1) = basicFunctions(i,1,2);
    if i == 1
        xcoord(n,1) = B1(1,1);
        ycoord(n,1) = B1(2,1);
    else
        xcoord((i-1),1) = B1(1,1);
        ycoord((i-1),1) = B1(2,1);
    end
end

xcent = mean(xcoord);
ycent = mean(ycoord);

for i = 1:pointNum
    if ~isnan(tgStrength(i))
        radVec = [(targetPoints(i)-xcent);(targetPoints(i)-ycent)];
        tmpNormVec = [tgUnitVec(i,2); -tgUnitVec(i,1)];
        tmpDot = dot(radVec, tmpNormVec);
        if tmpDot~=0
            tmpNormVec2 = -tmpNormVec;
            if tmpDot > 0
                if outFlag == true
                    normUnitVec(i,:) = tmpNormVec;
                else
                    normUnitVec(i,:) = tmpNormVec2;
                end
            else
                if outFlag == true
                    normUnitVec(i,:) = tmpNormVec2;
                else
                    normUnitVec(i,:) = tmpNormVec;
                end
            end
        else
            normUnitVec(i,:) = tmpNormVec;
        end
        
    else
        normUnitVec(i,:) = [NaN NaN];
        normStrength(i) = NaN;
    end
end

normStrength(:,:) = tgStrength(:,:);

if plotFlag == true
    hold on;
    xcoord = targetPoints(:,1);
    ycoord = targetPoints(:,2);
    xvec = normUnitVec(:,1);
    yvec = normUnitVec(:,2);
    
    xcoord = xcoord(~isnan(normStrength));
    ycoord = ycoord(~isnan(normStrength));
    xvec = xvec(~isnan(normStrength));
    yvec = yvec(~isnan(normStrength));
    quiver(xcoord, ycoord, xvec, yvec);
end