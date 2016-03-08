function [tgUnitVec, tgStrength] = nodeTangent(basicFunctions, opt)

n = size(basicFunctions,1);
tgUnitVec = double(zeros(n,2));
tgStrength = double(zeros(n,1));

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

B2 = double(zeros(2,1));
B1 = double(zeros(2,1)); 
xcoord = double(zeros(n,1));
ycoord = double(zeros(n,1));

for i = 1:n
    B2(1,1) = basicFunctions(i,2,1);
    B2(2,1) = basicFunctions(i,2,2);
    tmp_norm = norm(B2);
    if i == 1
        tgStrength(n,1) = tmp_norm;
        B2 = B2./tmp_norm;
        tgUnitVec(n,:) = B2(:,1);
    else
        tgStrength((i-1),1) = tmp_norm;
        B2 = B2./tmp_norm;
        tgUnitVec((i-1),:) = B2(:,1);
    end
    
    if plotFlag == true
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
end

if plotFlag == true
    hold on;
    xvec = tgUnitVec(:,1);
    yvec = tgUnitVec(:,2);
    quiver(xcoord, ycoord, xvec, yvec);
end