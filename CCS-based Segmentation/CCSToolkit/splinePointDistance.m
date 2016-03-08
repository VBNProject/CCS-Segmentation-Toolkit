function [minPointList, distanceList] = splinePointDistance(basicFunctions, targetPoints, opt)

n = size(basicFunctions,1);
len = size(targetPoints,1);
distanceList = double(zeros(len,1));
minPointList = double(zeros(len,2));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

try
    fastFlag = opt.fastFlag;
catch
    fastFlag = false;
end

try
    secondFlag = opt.secondFlag;
catch
    secondFlag = false;
end

exOpt.plotFlag = false;
pointList = extractNodePoints(basicFunctions, exOpt);

for i = 1:len
    px = targetPoints(i,1);
    py = targetPoints(i,2);
    if fastFlag 
        pxArray = repmat(px,n,1);
        pyArray = repmat(py,n,1);
        pointDist = ((pxArray-pointList(:,1)).^2+(pyArray-pointList(:,2)).^2).^0.5;
        [~,minInd] = min(pointDist);
        B = double(zeros(2,4));
        C = double(zeros(2,4));
        for j = 1:4
            B(1,j) = basicFunctions(minInd,j,1);
            B(2,j) = basicFunctions(minInd,j,2);
        end
        if minInd == n
            for j = 1:4
                C(1,j) = basicFunctions(1,j,1);
                C(2,j) = basicFunctions(1,j,2);
            end                     
        else
            for j = 1:4
                C(1,j) = basicFunctions((minInd+1),j,1);
                C(2,j) = basicFunctions((minInd+1),j,2);
            end                
        end

        fLeft = @(t)(B(1,1) + t.*B(1,2) + (t*t).*B(1,3) + (t^3).*B(1,4)-px).^2 + (B(2,1) + t.*B(2,2) + (t*t).*B(2,3) + (t^3).*B(2,4)-py).^2;
        [tLeftVal, fLeftVal] = fminbnd(fLeft,0,1);
        fRight = @(t)(C(1,1) + t.*C(1,2) + (t*t).*C(1,3) + (t^3).*C(1,4)-px).^2 + (C(2,1) + t.*C(2,2) + (t*t).*C(2,3) + (t^3).*C(2,4)-py).^2;
        [tRightVal, fRightVal] = fminbnd(fRight,0,1);
        if fLeftVal < fRightVal
            distanceList(i,1) = sqrt(fLeftVal);
            minPointList(i,1) = B(1,1) + tLeftVal.*B(1,2) + (tLeftVal*tLeftVal).*B(1,3) + (tLeftVal^3).*B(1,4);
            minPointList(i,2) = B(2,1) + tLeftVal.*B(2,2) + (tLeftVal*tLeftVal).*B(2,3) + (tLeftVal^3).*B(2,4);
        else
            distanceList(i,1) = sqrt(fRightVal);
            minPointList(i,1) = C(1,1) + tRightVal.*C(1,2) + (tRightVal*tRightVal).*C(1,3) + (tRightVal^3).*C(1,4);
            minPointList(i,2) = C(2,1) + tRightVal.*C(2,2) + (tRightVal*tRightVal).*C(2,3) + (tRightVal^3).*C(2,4);
        end 
    else
        localDistanceList = double(zeros(n,1));
        localMinPointList = double(zeros(n,2));
        for j = 1:n
            B = double(zeros(2,4));
            for k = 1:4 
                B(1,k) = basicFunctions(j,k,1);
                B(2,k) = basicFunctions(j,k,2);
            end
            f = @(t)(B(1,1) + t.*B(1,2) + (t*t).*B(1,3) + (t^3).*B(1,4)-px).^2 + (B(2,1) + t.*B(2,2) + (t*t).*B(2,3) + (t^3).*B(2,4)-py).^2;
            [tVal, fVal] = fminbnd(f,0,1);
            localDistanceList(j,1) = sqrt(fVal);
            localMinPointList(j,1) = B(1,1) + tVal.*B(1,2) + (tVal*tVal).*B(1,3) + (tVal^3).*B(1,4);
            localMinPointList(j,2) = B(2,1) + tVal.*B(2,2) + (tVal*tVal).*B(2,3) + (tVal^3).*B(2,4);
        end
        if secondFlag == false
            [~,minInd] = min(localDistanceList);
            distanceList(i,1) = localDistanceList(minInd);
            minPointList(i,:) = localMinPointList(minInd,:);
        else
            [sortDist,sortInd] = sort(localDistanceList);
            sortInd = sortInd(sortDist>1);
            secondInd = sortInd(1);
            distanceList(i,1) = localDistanceList(secondInd);
            minPointList(i,:) = localMinPointList(secondInd,:);
        end
    end
end

if plotFlag
    figure;
    for i = 1:n
        hold on;
        plot(pointList(i,1),pointList(i,2), 'r+');
    end
    for i = 1:len
        hold on;
        plot(targetPoints(i,1), targetPoints(i,2), 'g+');
        hold on;
        plot(minPointList(i,1), minPointList(i,2), 'g+');
        hold on;
        plot([targetPoints(i,1), minPointList(i,1)], [targetPoints(i,2), minPointList(i,2)], 'g');
    end
end