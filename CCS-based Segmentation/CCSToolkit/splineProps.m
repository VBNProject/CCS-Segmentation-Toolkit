function outputList = splineProps(basicFunctions, targetPoints, properties)

pointNum = size(targetPoints,1);

if (strcmp(properties, 'd')||strcmp(properties, 'derative'))
    calcFlag = 1;
elseif (strcmp(properties, 'c')||strcmp(properties, 'curvature'))
    calcFlag = 2;
elseif (strcmp(properties, 'r')||strcmp(properties, 'radius'))
    calcFlag = 3;
elseif (strcmp(properties, 'a')||strcmp(properties, 'arc'))
    calcFlag = 4;
elseif (strcmp(properties, 'ra')||strcmp(properties, 'radian'))
    calcFlag = 5;
elseif (strcmp(properties, 'dc')||strcmp(properties, 'derative of curvature'))
    calcFlag = 6;
elseif (strcmp(properties, 'dr')||strcmp(properties, 'derative of radius'))
    calcFlag = 7;
else
    calcFlag = 1;
end

opt.plotFlag = false;
[onSplineFlag, segNum, segT] = isOnSpline(basicFunctions, targetPoints, opt);

outputFlag = true(pointNum,1);
derativeList = double(zeros(pointNum,2));
arcList = double(zeros(pointNum,1));
radianList = double(zeros(pointNum,1));
curvatureList = double(zeros(pointNum,1));
radiusList = double(zeros(pointNum,1));
decurvaList = double(zeros(pointNum,1));
deradiusList = double(zeros(pointNum,1));

for i = 1:pointNum
    if onSplineFlag(i,1)== false
        outputFlag(i,:) = false;
        derativeList(i,:) = [NaN NaN];
        arcList(i,1) = NaN;
        radianList(i,1) = NaN;
        curvatureList(i,1) = NaN;
        radiusList(i,1) = NaN;
    else
        b2x = basicFunctions(segNum(i),2,1);
        b2y = basicFunctions(segNum(i),2,2);
        b3x = basicFunctions(segNum(i),3,1);
        b3y = basicFunctions(segNum(i),3,2);
        b4x = basicFunctions(segNum(i),4,1);
        b4y = basicFunctions(segNum(i),4,2);
        t = segT(i);
                 
        partialX = b2x + 2*b3x*t + 3*b4x*t*t;
        partialY = b2y + 2*b3y*t + 3*b4y*t*t;
        partialX2 = 2*b3x + 6*b4x*t;
        partialY2 = 2*b3y + 6*b4y*t;
        partialX3 = 6*b4x;
        partialY3 = 6*b4y;
        derativeList(i,:) = [partialX partialY];
        
        if partialY2*partialX > partialY*partialX2
            tmp1 = partialY2*partialX - partialY*partialX2;
            tmp1_derative = partialY3*partialX - partialY*partialX3;
        else
            tmp1 = partialY*partialX2 - partialY2*partialX;
            tmp1_derative = partialY*partialX3 - partialY3*partialX;
        end
        tmp2 = partialX.^2 + partialY.^2;
        tmp3 = partialX*partialX2 + partialY*partialY2;
        
        partialArc = tmp2.^0.5;
        partialArc2 = tmp3/(tmp2.^0.5);
        partialRad = double(abs(tmp1))/double(tmp2);
        partialRad2 = double(tmp1_derative*tmp2 - 2*tmp3*tmp1)/double(tmp2.^2);
        
        arcList(i,1) = partialArc;
        radianList(i,1) = partialRad;

        if calcFlag == 2 
            curvatureList(i,1) = partialRad/partialArc;
        end

        if calcFlag == 3 
            radiusList(i,1) = partialArc/partialRad;
        end

        if calcFlag == 6
            numerator = partialRad2*partialArc - partialArc2*partialRad;
            denominator = partialArc.^2;
            decurvaList(i,1) = numerator/denominator;
        end

        if calcFlag == 7
            numerator = partialArc2*partialRad - partialRad2*partialArc;
            denominator = partialRad.^2;
            deradiusList(i,1) = numerator/denominator;
        end
    end
end

if calcFlag == 1
    outputList = derativeList;
elseif calcFlag == 2
    outputList = curvatureList;
elseif calcFlag == 3
    outputList = radiusList;
elseif calcFlag == 4
    outputList = arcList;
elseif calcFlag == 5
    outputList = radianList;
elseif calcFlag == 6
    outputList = decurvaList;
elseif calcFlag == 7
    outputList = deradiusList;
end
