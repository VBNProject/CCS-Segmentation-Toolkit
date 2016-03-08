function [labelImg, labelNum] = closedCubicSpline2label(imwidth, imheight, basicFunctions, opt)

try
    plotFlag = opt.plotFlag;
catch 
    plotFlag = false;
end
inopt.plotFlag = false;

enhanceFactor = 3.0;
pointList = extractNodePoints(basicFunctions);
newPointList = double(zeros(size(pointList,1),2));
newPointList(1,:) = pointList(end,:);
newPointList(2:end,:) = pointList(1:(size(pointList,1)-1),:);
diffPointList = abs(newPointList - pointList);
maxGap = max(max(diffPointList));
sampleNum = enhanceFactor*maxGap*size(pointList,1);
outputPoints = spline2Points(basicFunctions, sampleNum);

labelImg = uint8(zeros(imheight, imwidth));
for i = 1:size(outputPoints,1)
    xcoord = uint16(outputPoints(i,1));
    ycoord = uint16(outputPoints(i,2));
    labelImg(ycoord, xcoord) = 255;
end
tmp = im2bw(labelImg);
sampleNum = size(find(tmp==1),1);
tmp = imfill(tmp,'holes');
labelImg = 255*uint8(tmp);

labelNum = size(find(labelImg==255),1);
if labelNum <= sampleNum
targetPoints = double(zeros(1,2));
opt.insidePoint = mean(pointList);
    for i = 1:size(outputPoints,1)
        xcoord = uint16(outputPoints(i,1));
        ycoord = uint16(outputPoints(i,2));
        localxmin = max(1, (ycoord-1));
        localxmax = min(imheight, (ycoord+1));
        localymin = max(1, (xcoord-1));
        localymax = min(imwidth, (xcoord+1)); 
        for j = localxmin:localxmax
            for k = localymin:localymax
                if labelImg(j,k) ~= 255
                    targetPoints(1,1) = k;
                    targetPoints(1,2) = j;
                    if isInsideSpline(basicFunctions, targetPoints, opt)
                        labelImg(j,k) = 255;
                    end                    
                end
            end
        end
    end

tmp = im2bw(labelImg);
tmp = imfill(tmp,'holes');
labelImg = 255*uint8(tmp);    
end

labelNum = size(find(labelImg==255),1);

if plotFlag == true
    plotRes = 20;
    CCSGroup = {basicFunctions};
    plotSpline(labelImg, plotRes, CCSGroup, inopt);
end
