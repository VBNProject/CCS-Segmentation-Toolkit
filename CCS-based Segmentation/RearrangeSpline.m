function newBasicFunctions = RearrangeSpline2(basicFunctions, lowerThresh, upperThresh, rearrageDepth, opt)
% 用途：重整三次样条曲线基函数组，归并在空间中距离过近的节点，在缺少节点的区域添加节点
% 该函数需要依赖CCS Toolkit工具包
% 输入参数
% basicFunctions：n*4*2形式的数组，存储表征三次样条曲线的基函数组
% 对于每一个分段i，其所对应的基函数组适用于点[i-1 i]所构成的区间
% filterThresh：判断CCS的两个节点是否靠得过近的阈值，默认值为5个像素
% 该值越大，则需要优化和重排的节点越多，曲线越不容易出现缠结，但曲线所保留的细节越少
% opt：绘图控制参数
% opt.plotFlag：判断用户是否需要绘制重整前后的三次样条曲线
% opt.integerFlag：是否强制新添加的点的坐标为整数的标识位，默认为false
% opt.backgroundImage：用于对比重整前后的曲线的背景图片
% 输出参数
% newbasicFunctions：重整后的三次样条曲线基函数组
% 注意：对于从label2ClosedCubicSpline函数提取出的CCS，各个节点都是label的边界点，如果直接使用本函数进行重整，新添加的节点可能会偏离label的边界

%% 处理基本参数
% 控制节点对删除距离的参数，默认为2。
% 当rearrageDepth取值为k，表示节点编号距离小于等于k的节点对，如果相距距离小于lth，则以这个节点对为起止的所有节点都会被删除
% 一般而言，当图像序列的相邻图像差异不大时，取2即可。如果相邻图像差异较大，可以考虑取更大的值，防止扭曲缠结出现
try
    rd = rearrageDepth; 
    if isempty(rd)
        rd = 2;
    end
catch
    rd = 2;
end

% 相邻节点欧氏距离的上限与下限
try
    lth = lowerThresh;
    if isempty(lth)
        lth = 5;
    end
catch  % 如果用户没有输入
    lth = 5;
end

try
    uth = upperThresh;
    if isempty(uth)
        lth = 10;
    end
catch  % 如果用户没有输入
    uth = 10;
end

try
    integerFlag = opt.integerFlag;
catch  % 如果用户没有输入
    integerFlag = false; % 默认不强制改插入点的坐标为整数
end

% 绘图控制参数
try
    plotFlag = opt.plotFlag;
catch  % 如果用户没有输入
    plotFlag = false; % 默认不绘制
end

if plotFlag
    I = opt.backgroundImage;
end

%% 剔除与CCS靠得过近的相邻点
inOpt.secondFlag = true;
opt.plotFlag = false;
pointList = extractNodePoints(basicFunctions, opt);

[~, distanceList] = splinePointDistance(basicFunctions, pointList, inOpt);
pointLabel = distanceList>lth;
pointList = pointList(pointLabel,:);

%% 剔除靠得过近的i个相邻点
n = size(pointList,1);
dist = double(zeros(n,1));
for i = 1:n
    if i == n
        dist(i) = sqrt( (pointList(i,1)-pointList(1,1))^2 + (pointList(i,2)-pointList(1,2))^2 );
    else
        dist(i) = sqrt( (pointList(i,1)-pointList(i+1,1))^2 + (pointList(i,2)-pointList(i+1,2))^2 );
    end
end

pointLabel = dist>lth;
pointList = pointList(pointLabel,:);

%% 剔除靠得过近的任意点
pointList = filterCloseNodes(pointList, lth, rd);

%% 添加若干个节点，保证没有相隔过远的节点
newPointList = addNewNodes(pointList, uth, opt, integerFlag);

%% 再过滤一次靠得过近的节点，防止新添加节点的错误
newPointList = filterCloseNodes(newPointList, lth, rd);

%% 得到新的样条基函数组
newBasicFunctions = closedCubicSpline(newPointList,opt);

%% 结果绘制部分
if plotFlag
    opt.vectorFlag = false;
    plotRes = 10;
    CCSGroup = {basicFunctions; newBasicFunctions}; % 注意，构造的是2*1的单元数组，中间需要加上分号，而不是逗号
    plotSpline(I, plotRes, CCSGroup, opt);
end
end

%% 过滤靠得过近的任意节点的函数
function pointList = filterCloseNodes(pointList, lth, rd)
% 获取节点距离矩阵
n = size(pointList,1);
dist = double(zeros((n-1), n)); % 节点距离矩阵，第i行记录了每个节点左侧距该节点编号偏移为i的节点。第j列则记录了原始节点列表中的第j个节点
for i = 1:(n-1) % 第i次循环，得到每个节点与其左侧第i个节点的距离
    % 循环左移i个节点，得到新节点列表
    cmpPointList = double(zeros(n,2));
    cmpPointList(1:(n-i), :) = pointList((i+1):end,:);
    cmpPointList((n-i+1):end, :) = pointList(1:i,:);
    dist(i,:) = (((pointList(:,1)-cmpPointList(:,1)).^2 + (pointList(:,2)-cmpPointList(:,2)).^2).^0.5)';
end

% 提取所有小于距离阈值的节点对
closePPairs = uint16(zeros(0,2));
for i = 1:(n-1)
    for j = 1:n
        if dist(i,j)<=lth
            p1 = j;
            tmp = j - i;
            if tmp <=0
                tmp = tmp + n;
            end
            p2 = tmp;
            if p1>p2 % 按照节点对中的前一个点的编号比后一个点小的顺序存入
                closePPairs = [closePPairs; p2 p1];
            else
                closePPairs = [closePPairs; p1 p2];
            end
        end
    end
end

% 根据rearrageDepth的要求，删除以节点对为起止的节点段
if ~isempty(closePPairs)
    saveLabel = true(n,1); % 对待保留节点进行标记
    closePDist = abs(double(closePPairs(:,2))-double(closePPairs(:,1))); % 计算所有待删除节点对的两点之间的距离
    closePPairs = closePPairs(closePDist <= rd, :); % 根据节点的编号相差的数值，挑选出待删除的节点对
    for i = 1:size(closePPairs,1)
        p1 = closePPairs(i,1);
        p2 = closePPairs(i,2);
        dist1 = double(p2) - double(p1) -1.0;
        dist2 = double(p1) - double(p2) -1.0 + n;
        if dist1<=dist2
            saveLabel((p1+1):p2, :) = false;
        else
            saveLabel(p2:end, :) = false;
            if p1>1
                saveLabel(1:(p1-1), :) = false;
            end
        end
    end
    
    pointList = pointList(saveLabel, :);
end
end

%% 添加新节点的函数
function newPointList = addNewNodes(pointList, uth, opt, integerFlag)
% 计算删除了过近节点后，各个节点与右侧节点的距离
n = size(pointList,1);
dist = double(zeros(n,1));
for i = 1:n
    if i == 1
        dist(i) = sqrt( (pointList(i,1)-pointList(n,1))^2 + (pointList(i,2)-pointList(n,2))^2 );
    else
        dist(i) = sqrt( (pointList(i,1)-pointList(i-1,1))^2 + (pointList(i,2)-pointList(i-1,2))^2 );
    end
end
% 确认是否需要添加节点
pointLabel = dist<=uth;
addNum = length(find(pointLabel==0));
% 如果需要添加节点
if addNum~=0
    % 通过线性插值添加点
    newPointList = double(zeros((addNum+n),2));
    tmpBasicFunctions = closedCubicSpline(pointList,opt);
    cursor = 0; % 指向新节点列表中当前节点的上一个节点的游标
    for i = 1:n
        if dist(i)>uth
            % 计算新节点坐标
            B1 = double(zeros(2,1));
            B2 = double(zeros(2,1));
            B3 = double(zeros(2,1));
            B4 = double(zeros(2,1));
            B1(1,1) = tmpBasicFunctions(i,1,1);
            B1(2,1) = tmpBasicFunctions(i,1,2);
            B2(1,1) = tmpBasicFunctions(i,2,1);
            B2(2,1) = tmpBasicFunctions(i,2,2);
            B3(1,1) = tmpBasicFunctions(i,3,1);
            B3(2,1) = tmpBasicFunctions(i,3,2);
            B4(1,1) = tmpBasicFunctions(i,4,1);
            B4(2,1) = tmpBasicFunctions(i,4,2);
            t = 0.5;
            newPoint = (B1 + t.*B2 + (t*t).*B3 + (t^3).*B4)'; % 通过三次样条插值计算出位于两节点“正中”的点的坐标
            if integerFlag
                newPoint = uint16(newPoint);
            end 
            % 添加新节点
            if i == 1
                newPointList(end,:) = newPoint(:,:);
                newPointList((cursor+1),:) = pointList(i,:);
                cursor = cursor + 1;
            else
                newPointList((cursor+1), :) = newPoint(:,:);
                newPointList((cursor+2), :) = pointList(i,:);
                cursor = cursor + 2;
            end
        else
            newPointList((cursor+1),:) = pointList(i,:);
            cursor = cursor + 1;
        end
    end
else
    newPointList = pointList;    
end
end