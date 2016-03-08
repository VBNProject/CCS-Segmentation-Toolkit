function [intersectPoints, intSegNumber, orgSegNumber, intTVal, orgTVal] = splineSplineIntersection(intBasicFunctions, orgBasicFunctions, opt)

intn = size(intBasicFunctions,1);
orgn = size(orgBasicFunctions,1);

try
    eps = opt.eps;
catch
    eps = 1e-6;
end

intersectPoints = double(zeros(0,2));
intSegNumber = uint16(zeros(0,1));
orgSegNumber = uint16(zeros(0,1));
intTVal = double(zeros(0,1));
orgTVal = double(zeros(0,1));

for i = 1:intn
    intbonex = intBasicFunctions(i,1,1);
    intboney = intBasicFunctions(i,1,2);
    intbtwox = intBasicFunctions(i,2,1);
    intbtwoy = intBasicFunctions(i,2,2);
    intbthrx = intBasicFunctions(i,3,1);
    intbthry = intBasicFunctions(i,3,2);
    intbfoux = intBasicFunctions(i,4,1);
    intbfouy = intBasicFunctions(i,4,2);
    for j = 1:orgn
        orgbonex = orgBasicFunctions(j,1,1);
        orgboney = orgBasicFunctions(j,1,2);
        orgbtwox = orgBasicFunctions(j,2,1);
        orgbtwoy = orgBasicFunctions(j,2,2);
        orgbthrx = orgBasicFunctions(j,3,1);
        orgbthry = orgBasicFunctions(j,3,2);
        orgbfoux = orgBasicFunctions(j,4,1);
        orgbfouy = orgBasicFunctions(j,4,2);
        % k11-k12 + k21*t1-k22*t2 + k31*t1^2-k32*t2^2 + k41*t^3-k42*t^3 = 0
        konex = intbonex - orgbonex;
        koney = intboney - orgboney;     
        fun = @(t)[konex + (intbtwox*t(1)-orgbtwox*t(2)) + (intbthrx*t(1)^2-orgbthrx*t(2)^2) + (intbfoux*t(1)^3-orgbfoux*t(2)^3), koney + (intbtwoy*t(1)-orgbtwoy*t(2)) + (intbthry*t(1)^2-orgbthry*t(2)^2) + (intbfouy*t(1)^3-orgbfouy*t(2)^3)];
        t0 = [0.5,0.5];
        options = optimoptions('fsolve','Display','off');
        [t, fval] = fsolve(fun,t0,options);
        if (t(1)>=0)&&(t(1)<=1)&&(t(2)>=0)&&(t(2)<=1)
            if (abs(fval(1))<eps)&&(abs(fval(2))<eps)
                xcoord = intbonex + intbtwox*t(1) + intbthrx*t(1)^2 + intbfoux*t(1)^3;
                ycoord = intboney + intbtwoy*t(1) + intbthry*t(1)^2 + intbfouy*t(1)^3;
                len = size(intersectPoints,1);
                tmpIntersectPoints = intersectPoints;
                tmpSegNumber = intSegNumber;
                tmpOrgSegNumber = orgSegNumber;
                tmpIntTVal = intTVal;
                tmpOrgTVal = orgTVal;
                intersectPoints = double(zeros((len+1),2));
                intSegNumber = uint16(zeros((len+1),1));
                orgSegNumber = uint16(zeros((len+1),1));
                intTVal = double(zeros((len+1),1));
                orgTVal = double(zeros((len+1),1));
                if len ~=0
                    intersectPoints(1:len,:) = tmpIntersectPoints(:,:);
                    intersectPoints((len+1),:) = [xcoord ycoord];
                    intSegNumber(1:len,:) = tmpSegNumber;
                    intSegNumber((len+1),:) = i;
                    orgSegNumber(1:len,:) = tmpOrgSegNumber;
                    orgSegNumber((len+1),:) = j;
                    intTVal(1:len,:) = tmpIntTVal;
                    intTVal((len+1),:) = t(1);
                    orgTVal(1:len,:) = tmpOrgTVal;
                    orgTVal((len+1),:) = t(2);
                else
                    intersectPoints(1,:) = [xcoord ycoord];
                    intSegNumber(1,:) = i;
                    orgSegNumber(1,:) = j;
                    intTVal = t(1);
                    orgTVal = t(2);
                end
            end
        end
    end
end

len = size(intSegNumber,1);
markFlag = true(len,1);
for i = 1:len
    curPoint = intersectPoints(i,:);
    if i ~= len
        nxtPoint = intersectPoints((i+1),:);
    else
        nxtPoint = intersectPoints(1,:);
    end
    if (abs(curPoint(1,1)-nxtPoint(1,1))<eps)&&(abs(curPoint(1,2)-nxtPoint(1,2))<eps)
        if i~=len
            markFlag((i+1),1) = false;
        else
            markFlag(1,1) = false;
        end       
    end
end

intersectPoints = intersectPoints(markFlag,:);
intSegNumber = intSegNumber(markFlag);
orgSegNumber = orgSegNumber(markFlag);
intTVal = intTVal(markFlag);
orgTVal = orgTVal(markFlag);