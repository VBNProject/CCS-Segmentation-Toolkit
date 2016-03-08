function pointList = extractPointList(basicFunctions, res, opt)

n = size(basicFunctions,1);

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

pointList = double(zeros((n*(res+1)),2));
for i = 1:n    
    B1 = double(zeros(2,1));
    B2 = double(zeros(2,1));
    B3 = double(zeros(2,1));
    B4 = double(zeros(2,1));
    B1(1,1) = basicFunctions(i,1,1);
    B1(2,1) = basicFunctions(i,1,2);
    B2(1,1) = basicFunctions(i,2,1);
    B2(2,1) = basicFunctions(i,2,2);
    B3(1,1) = basicFunctions(i,3,1);
    B3(2,1) = basicFunctions(i,3,2);
    B4(1,1) = basicFunctions(i,4,1);
    B4(2,1) = basicFunctions(i,4,2);
    
    tVec = linspace(0,1,(res+2));
    tVec = tVec';
    len = size(tVec,1);
    pointCoord = double(zeros(len,2));
    for j = 1:len
        t = tVec(j);
        pointCoord(j,:) = B1 + t.*B2 + (t*t).*B3 + (t^3).*B4;
    end
    
    if i == 1
        startIndex = (n-1)*(res+1) + 1;
        terminIndex = n*(res+1);
    else
        startIndex = (i-2)*(res+1) + 1;
        terminIndex = (i-1)*(res+1);
    end
    pointList(startIndex:terminIndex,:) = pointCoord(1:(res+1),:);
end

if plotFlag == true
    figure;
    for i = 1:size(pointList,1)
        hold on;
        plot(pointList(i,1), pointList(i,2),'g+');
    end
end