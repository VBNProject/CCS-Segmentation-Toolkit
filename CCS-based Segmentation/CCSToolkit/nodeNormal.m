function normUnitVec = nodeNormal(basicFunctions, outFlag, opt)

n = size(basicFunctions,1);
normUnitVec = double(zeros(n,2));
tgUnitVec = double(zeros(n,2));

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
        B2 = B2./tmp_norm;
        tgUnitVec(n,:) = B2(:,1);
    else
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

xcent = mean(xcoord);
ycent = mean(ycoord);

for i = 1:n
    radVec = [(xcoord(i)-xcent);(ycoord(i)-ycent)];
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
end

if plotFlag == true
    hold on;
    xvec = normUnitVec(:,1);
    yvec = normUnitVec(:,2);
    quiver(xcoord, ycoord, xvec, yvec);    
end

end