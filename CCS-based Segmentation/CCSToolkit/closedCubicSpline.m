function basicFunctions = closedCubicSpline(pointList, opt)

pointList = double(pointList);
n = size(pointList,1);

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

if plotFlag == true
    try
        plotRes = opt.res;
    catch
        plotRes = 10;
    end
end

a=diag(4*ones(1,n));
b = diag(ones(1,n),-1); 
b = b(1:n,1:n);
c = diag(ones(1,n),1); 
c = c(1:n,1:n);
A = a + b+ c;
A(n,1) = 1;
A(1,n) = 1;

upList = [pointList(2:end,:); pointList(1,:)];
doList = [pointList(n,:); pointList(1:(n-1),:)];
B = 3.0*upList - 3.0*doList;

Deratives = A\B;

% b1 + b2*t + b3*t^2 + b4*t^3
basicFunctions = double(zeros(n,4,2));

basicFunctions(1,1,:) = pointList(n,:);
basicFunctions(1,2,:) = Deratives(n,:);
basicFunctions(1,3,:) = (-1).*Deratives(1,:) + (-2).*Deratives(n,:) + 3.*pointList(1,:) + (-3).*pointList(n,:);
basicFunctions(1,4,:) = Deratives(1,:) + Deratives(n,:) + (-2).*pointList(1,:) + 2.*pointList(n,:);
for i = 2:n
    basicFunctions(i,1,:) = pointList((i-1),:);
    basicFunctions(i,2,:) = Deratives((i-1),:);
    basicFunctions(i,3,:) = (-1).*Deratives(i,:) + (-2).*Deratives((i-1),:) + 3.*pointList(i,:) + (-3).*pointList((i-1),:);
    basicFunctions(i,4,:) = Deratives(i,:) + Deratives((i-1),:) + (-2).*pointList(i,:) + 2.*pointList((i-1),:);
end

if plotFlag == true
    tVec = linspace(0,1,(plotRes+2));
    tVec = tVec';
    B1 = double(zeros(2,1));
    B2 = double(zeros(2,1));
    B3 = double(zeros(2,1));
    B4 = double(zeros(2,1));
    
    figure;
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
        
        hold on;
        plot(pointCoord(:,1),pointCoord(:,2),'r');
    end
end