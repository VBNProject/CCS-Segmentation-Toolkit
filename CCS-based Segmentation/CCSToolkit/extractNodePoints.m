function pointList = extractNodePoints(basicFunctions, opt)

n = size(basicFunctions,1);

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

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

pointList = [xcoord ycoord];

if plotFlag == true
    hold on;
    for i = 1:n
        plot(xcoord(i),ycoord(i),'g+');
    end
end