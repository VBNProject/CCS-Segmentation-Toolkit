function outputPoints = spline2Points(basicFunctions, sampleNum, opt)

n = size(basicFunctions,1);
perNum = round(sampleNum/n);
outputPoints = double(zeros((perNum*n),2));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

tVec = linspace(0,1,(perNum+1));
tVec = tVec';
B1 = double(zeros(2,1));
B2 = double(zeros(2,1));
B3 = double(zeros(2,1));
B4 = double(zeros(2,1));

for i = 1:n
    B1(1,1) = basicFunctions(i,1,1);
    B1(2,1) = basicFunctions(i,1,2);
    B2(1,1) = basicFunctions(i,2,1);
    B2(2,1) = basicFunctions(i,2,2);
    B3(1,1) = basicFunctions(i,3,1);
    B3(2,1) = basicFunctions(i,3,2);
    B4(1,1) = basicFunctions(i,4,1);
    B4(2,1) = basicFunctions(i,4,2);
    
    len = size(tVec,1);
    pointCoord = double(zeros(len,2));
    for j = 1:len
        t = tVec(j);
        pointCoord(j,:) = B1 + t.*B2 + (t*t).*B3 + (t^3).*B4;
    end

    if i==1
        pShift = (n-1)*perNum+1;
    else
        pShift = (i-2)*perNum+1;
    end
    outputPoints(pShift:(pShift+perNum-1),:) = pointCoord(1:perNum,:);
end

if plotFlag == true
    len = size(outputPoints,1);
    for i = 1:len
        hold on;
        xcoord = outputPoints(i,1);
        ycoord = outputPoints(i,2);
        plot(xcoord,ycoord,'bo');
    end
end

end