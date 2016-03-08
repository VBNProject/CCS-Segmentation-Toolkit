%% Configuration
configureKits();

windowSize = 3;
windowWidth = 10;
windowHeight = 5;
filledValue = 0;

rearrangeLowTh = 5;
rearrangeHighTh = 20;
rearrangeDepth = 4;

opt.nodeFlag = true;
opt.plotFlag = false;
opt.cloThr = 5;

nxtStart = 8768;
nxtTermi = 8867;
srcDir = '../realImage/12458_CENTmo/';
baseSrc = [srcDir num2str(nxtStart, '%05d')  '_label.tif'];

orgSrc = [srcDir num2str(nxtStart, '%05d')  '.tif'];
nxtDir = srcDir;

%% Calculation
if matlabpool('size')<=0
    matlabpool local 12
end

ccsGroup = cell((nxtTermi-nxtStart+1),1);
labelImg = imread(baseSrc);
neighborSize = 5;
basicFunctions = label2ClosedCubicSpline(labelImg, neighborSize, 200, opt);
if size(basicFunctions,1)<3
    error('Not enough node points');
end
basicFunctions = RearrangeSpline(basicFunctions);
ccsGroup{1,1} = basicFunctions;

for n = (nxtStart+1):nxtTermi
    orgImg = imread([nxtDir num2str((n-1), '%05d') '.tif']);
    nxtImg = imread([nxtDir num2str(n, '%05d') '.tif']);
    orgImg = medfilt2(orgImg,[neighborSize neighborSize]);
    nxtImg = medfilt2(nxtImg,[neighborSize neighborSize]);
    
    oldBasicFunctions = ccsGroup{(n-nxtStart),1};
    nodeList = extractNodePoints(oldBasicFunctions);
    normList = nodeNormal(oldBasicFunctions, true, opt);
    
    nodeNum = size(nodeList,1);
    newNodeList = double(zeros(nodeNum,2));
      
    parfor i = 1:nodeNum
        curNode = nodeList(i,:);
        curNorm = normList(i,:);
        nodeInfo = [curNode, curNorm]';
        orgInfo = imfinfo(baseSrc);
        [croppedImage, ~] = NodeSurroundingExtraction(orgImg, nodeInfo, windowWidth, windowHeight, filledValue, orgInfo);
        [croppedImage2, invM] = NodeSurroundingExtraction(nxtImg, nodeInfo, windowWidth, (windowHeight+windowSize), filledValue, orgInfo);
        [global_pos, ~, local_out_img] = NormSearch(croppedImage, croppedImage2, invM);
        global_pos = global_pos(1:2,:);
        newNodeList(i,:) = global_pos';
    end
    
    newBasicFunctions = closedCubicSpline(newNodeList);    
    newBasicFunctions = RearrangeSpline(newBasicFunctions, rearrangeLowTh, rearrangeHighTh, rearrangeDepth);
    ccsGroup{(n-nxtStart+1),1} = newBasicFunctions;
    disp([num2str((n-nxtStart+1)) ': ' num2str(n)]);
end

%% Plotting the result
opt.nodeFlag = false;
opt.saveFlag = false;
opt.colormap = [1.0, 0.0, 0.0; 1.0, 0.0, 0.0];
plotSpline(nxtImg, 20, {newBasicFunctions}, opt);
opt.interval = 10;plotSpline(nxtImg, 20, ccsGroup, opt);
m = 20;plotSpline(imread([nxtDir num2str(nxtStart+m-1, '%05d') '.tif']), 20, ccsGroup(m), opt);