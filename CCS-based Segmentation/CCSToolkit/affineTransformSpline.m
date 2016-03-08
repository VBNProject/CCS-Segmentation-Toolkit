function newBasicFunctions = affineTransformSpline(basicFunctions, opt, varargin)

exOpt.plotFlag = false; 
pointList = extractNodePoints(basicFunctions, exOpt);

try
    transformCenter = opt.transformCenter;
catch
    transformCenter = 0;
end

preMat = eye(3);
postMat = eye(3);
if transformCenter == 1
    srOpt.plotFlag = false;
    [xmin,~,ymin,~,~] = splineRange(basicFunctions,srOpt);
    preMat(1,3) = -xmin;
    preMat(2,3) = -ymin;
    postMat(1,3) = xmin;
    postMat(2,3) = ymin;    
elseif transformCenter == 2
    centCoord = mean(pointList);
    xcent = centCoord(1,1);
    ycent = centCoord(1,2);
    preMat(1,3) = -xcent;
    preMat(2,3) = -ycent;
    postMat(1,3) = xcent;
    postMat(2,3) = ycent;
else
end


if nargin == 3
    transMat = varargin{1};
elseif nargin == 4
    if isempty(varargin{1})
        rotA = varargin{2};
        transMat = [cos(rotA) -sin(rotA) 0; sin(rotA) cos(rotA) 0; 0 0 1];
    else
        transMat = eye(3);
        transMat(1,3) = varargin{1};
        transMat(2,3) = varargin{2};
    end
elseif nargin == 5
    if isempty(varargin{1})
        transMat = eye(3);
        transMat(1,1) = varargin{2};
        transMat(2,2) = varargin{3};
    else
        newBasicFunctions = basicFunctions;
        return;
    end
else
    newBasicFunctions = basicFunctions;
    return;
end

nodeNum = size(pointList,1);
addedOne = ones(nodeNum,1);
tmpList = [pointList addedOne];
tmpList = tmpList';
totalMat = postMat*transMat*preMat;
tmpList = totalMat*tmpList; 
tmpList = tmpList';
newPointList = tmpList(:,1:2);

inOpt.plotFlag = false;
newBasicFunctions = closedCubicSpline(newPointList, inOpt);