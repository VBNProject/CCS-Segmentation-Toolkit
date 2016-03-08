function isClockwiseFlag = isClockwise(basicFunctions, opt)

n = size(basicFunctions,1);

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

try
    centerPoint = opt.centerPoint;
    xcent = centerPoint(1,1);
    ycent = centerPoint(2,1);
catch
    xcent = mean(xcoord);
    ycent = mean(ycoord);
end
pcent = [xcent ycent];

cloAcc = 0;
conAcc = 0;
for i = 1:n
    if i ~= n
        p1 = pointList(i,:);
        p2 = pointList((i+1),:);
    else
        p1 = pointList(n,:);
        p2 = pointList(1,:);        
    end
    vec1 = p1 - pcent;
    vec2 = p2 - pcent;
    evec1 = [vec1 0];
    evec2 = [vec2 0];
    crossPro = cross(evec1,evec2);
    if crossPro(1,3)<0
        conAcc = conAcc + 1;
    elseif crossPro(1,3)>0
        cloAcc = cloAcc + 1;
    end
end

if conAcc > cloAcc
    isClockwiseFlag = false;
else
    isClockwiseFlag = true;
end