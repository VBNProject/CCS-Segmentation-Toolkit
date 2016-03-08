function insideSplineFlag = isInsideSpline(basicFunctions, targetPoints, opt)

opt.plotFlag = false;
n = size(basicFunctions,1);
len = size(targetPoints,1);
insideSplineFlag = false(len,1);

B1 = double(zeros(2,1));
xcoord = double(zeros(n,1));
ycoord = double(zeros(n,1));

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

try
    insidePoint = opt.insidePoint;
    xcent = insidePoint(1,1);
    ycent = insidePoint(1,2);
catch
    xcent = mean(xcoord);
    ycent = mean(ycoord);
end

targetPointsx = targetPoints(:,1);
targetPointsy = targetPoints(:,2);
for i = 1:len
    segsParam = double(zeros(1,4));
    segsParam(1,3) = targetPointsx(i);
    segsParam(1,4) = targetPointsy(i);
    segsParam(1,1) = xcent;
    segsParam(1,2) = ycent;
    intersectPoints = segSplineIntersection(basicFunctions, segsParam, opt);   
    pointNum = size(intersectPoints{1,1},1);
    if ((pointNum == 0)||(pointNum == 2))
        insideSplineFlag(i,1) = true;
    elseif ((pointNum == 1)||(pointNum == 3))
        insideSplineFlag(i,1) = false;
    end
end