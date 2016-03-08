%% Configuration
configureKits();

windowSize = 20;
windowWidth = 80;
windowHeight = 40;
filledValue = 0;

opt.nodeFlag = true;
opt.plotFlag = false;
opt.cloThr = 5;

nxtStart = 41;
nxtTermi = 159;
ccsGroup = cell((nxtTermi-nxtStart+1),1);

baseSrc = ['../modelImage/cellColony1_base/cell_' num2str(nxtStart, '%03d')  '.tif'];
maskSrc = ['../modelImage/cellColony3_mask/cell_' num2str(nxtStart, '%03d')  '.tif'];
nxtDir = '../modelImage/cellColony3/cell_';

%% Calculation
newLabelImg = imread(baseSrc);
neighborSize = 5; 
basicFunctions = label2ClosedCubicSpline(newLabelImg, neighborSize, 50, opt);
if size(basicFunctions,1)<3  
    error('Not enough node points');
end
basicFunctions = RearrangeSpline(basicFunctions);
ccsGroup{1,1} = basicFunctions;


for n = (nxtStart+1):nxtTermi
    orgImg = imread([nxtDir num2str((n-1), '%03d') '.tif']);
    nxtImg = imread([nxtDir num2str(n, '%03d') '.tif']);
    orgImg = medfilt2(orgImg,[neighborSize neighborSize]);
    nxtImg = medfilt2(nxtImg,[neighborSize neighborSize]);
    
    oldBasicFunctions = ccsGroup{(n-nxtStart),1};
    nodeList = extractNodePoints(oldBasicFunctions);
    normList = nodeNormal(oldBasicFunctions, true, opt);
    
    nodeNum = size(nodeList,1);
    newNodeList = double(zeros(nodeNum,2));
    
    for i = 1:nodeNum
        curNode = nodeList(i,:);
        curNorm = normList(i,:);
        nodeInfo = [curNode, curNorm]';
        orgInfo = imfinfo([nxtDir num2str((n-1), '%03d') '.tif']);
        [croppedImage, ~] = NodeSurroundingExtraction(orgImg, nodeInfo, windowWidth, windowHeight, filledValue, orgInfo);
        [croppedImage2, invM] = NodeSurroundingExtraction(nxtImg, nodeInfo, windowWidth, (windowHeight+windowSize), filledValue, orgInfo);
        [global_pos, ~] = NormSearch(croppedImage, croppedImage2, invM);
        global_pos = global_pos(1:2,:);
        newNodeList(i,:) = global_pos';
    end
    
    newBasicFunctions = closedCubicSpline(newNodeList);
    disp(n);
    newBasicFunctions = RearrangeSpline(newBasicFunctions);
    ccsGroup{(n-nxtStart+1),1} = newBasicFunctions;
end

%% Plotting the result
plotSpline(nxtImg, 20, {newBasicFunctions}, opt);
opt.interval = 10;
plotSpline(nxtImg, 20, ccsGroup, opt);