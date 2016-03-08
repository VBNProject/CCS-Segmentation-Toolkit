function newBasicFunctions = RearrangeSpline2(basicFunctions, lowerThresh, upperThresh, rearrageDepth, opt)
% ��;�����������������߻������飬�鲢�ڿռ��о�������Ľڵ㣬��ȱ�ٽڵ��������ӽڵ�
% �ú�����Ҫ����CCS Toolkit���߰�
% �������
% basicFunctions��n*4*2��ʽ�����飬�洢���������������ߵĻ�������
% ����ÿһ���ֶ�i��������Ӧ�Ļ������������ڵ�[i-1 i]�����ɵ�����
% filterThresh���ж�CCS�������ڵ��Ƿ񿿵ù�������ֵ��Ĭ��ֵΪ5������
% ��ֵԽ������Ҫ�Ż������ŵĽڵ�Խ�࣬����Խ�����׳��ֲ��ᣬ��������������ϸ��Խ��
% opt����ͼ���Ʋ���
% opt.plotFlag���ж��û��Ƿ���Ҫ��������ǰ���������������
% opt.integerFlag���Ƿ�ǿ������ӵĵ������Ϊ�����ı�ʶλ��Ĭ��Ϊfalse
% opt.backgroundImage�����ڶԱ�����ǰ������ߵı���ͼƬ
% �������
% newbasicFunctions��������������������߻�������
% ע�⣺���ڴ�label2ClosedCubicSpline������ȡ����CCS�������ڵ㶼��label�ı߽�㣬���ֱ��ʹ�ñ�������������������ӵĽڵ���ܻ�ƫ��label�ı߽�

%% �����������
% ���ƽڵ��ɾ������Ĳ�����Ĭ��Ϊ2��
% ��rearrageDepthȡֵΪk����ʾ�ڵ��ž���С�ڵ���k�Ľڵ�ԣ����������С��lth����������ڵ��Ϊ��ֹ�����нڵ㶼�ᱻɾ��
% һ����ԣ���ͼ�����е�����ͼ����첻��ʱ��ȡ2���ɡ��������ͼ�����ϴ󣬿��Կ���ȡ�����ֵ����ֹŤ���������
try
    rd = rearrageDepth; 
    if isempty(rd)
        rd = 2;
    end
catch
    rd = 2;
end

% ���ڽڵ�ŷ�Ͼ��������������
try
    lth = lowerThresh;
    if isempty(lth)
        lth = 5;
    end
catch  % ����û�û������
    lth = 5;
end

try
    uth = upperThresh;
    if isempty(uth)
        lth = 10;
    end
catch  % ����û�û������
    uth = 10;
end

try
    integerFlag = opt.integerFlag;
catch  % ����û�û������
    integerFlag = false; % Ĭ�ϲ�ǿ�ƸĲ���������Ϊ����
end

% ��ͼ���Ʋ���
try
    plotFlag = opt.plotFlag;
catch  % ����û�û������
    plotFlag = false; % Ĭ�ϲ�����
end

if plotFlag
    I = opt.backgroundImage;
end

%% �޳���CCS���ù��������ڵ�
inOpt.secondFlag = true;
opt.plotFlag = false;
pointList = extractNodePoints(basicFunctions, opt);

[~, distanceList] = splinePointDistance(basicFunctions, pointList, inOpt);
pointLabel = distanceList>lth;
pointList = pointList(pointLabel,:);

%% �޳����ù�����i�����ڵ�
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

%% �޳����ù����������
pointList = filterCloseNodes(pointList, lth, rd);

%% ������ɸ��ڵ㣬��֤û�������Զ�Ľڵ�
newPointList = addNewNodes(pointList, uth, opt, integerFlag);

%% �ٹ���һ�ο��ù����Ľڵ㣬��ֹ����ӽڵ�Ĵ���
newPointList = filterCloseNodes(newPointList, lth, rd);

%% �õ��µ�������������
newBasicFunctions = closedCubicSpline(newPointList,opt);

%% ������Ʋ���
if plotFlag
    opt.vectorFlag = false;
    plotRes = 10;
    CCSGroup = {basicFunctions; newBasicFunctions}; % ע�⣬�������2*1�ĵ�Ԫ���飬�м���Ҫ���Ϸֺţ������Ƕ���
    plotSpline(I, plotRes, CCSGroup, opt);
end
end

%% ���˿��ù���������ڵ�ĺ���
function pointList = filterCloseNodes(pointList, lth, rd)
% ��ȡ�ڵ�������
n = size(pointList,1);
dist = double(zeros((n-1), n)); % �ڵ������󣬵�i�м�¼��ÿ���ڵ�����ýڵ���ƫ��Ϊi�Ľڵ㡣��j�����¼��ԭʼ�ڵ��б��еĵ�j���ڵ�
for i = 1:(n-1) % ��i��ѭ�����õ�ÿ���ڵ���������i���ڵ�ľ���
    % ѭ������i���ڵ㣬�õ��½ڵ��б�
    cmpPointList = double(zeros(n,2));
    cmpPointList(1:(n-i), :) = pointList((i+1):end,:);
    cmpPointList((n-i+1):end, :) = pointList(1:i,:);
    dist(i,:) = (((pointList(:,1)-cmpPointList(:,1)).^2 + (pointList(:,2)-cmpPointList(:,2)).^2).^0.5)';
end

% ��ȡ����С�ھ�����ֵ�Ľڵ��
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
            if p1>p2 % ���սڵ���е�ǰһ����ı�űȺ�һ����С��˳�����
                closePPairs = [closePPairs; p2 p1];
            else
                closePPairs = [closePPairs; p1 p2];
            end
        end
    end
end

% ����rearrageDepth��Ҫ��ɾ���Խڵ��Ϊ��ֹ�Ľڵ��
if ~isempty(closePPairs)
    saveLabel = true(n,1); % �Դ������ڵ���б��
    closePDist = abs(double(closePPairs(:,2))-double(closePPairs(:,1))); % �������д�ɾ���ڵ�Ե�����֮��ľ���
    closePPairs = closePPairs(closePDist <= rd, :); % ���ݽڵ�ı��������ֵ����ѡ����ɾ���Ľڵ��
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

%% ����½ڵ�ĺ���
function newPointList = addNewNodes(pointList, uth, opt, integerFlag)
% ����ɾ���˹����ڵ�󣬸����ڵ����Ҳ�ڵ�ľ���
n = size(pointList,1);
dist = double(zeros(n,1));
for i = 1:n
    if i == 1
        dist(i) = sqrt( (pointList(i,1)-pointList(n,1))^2 + (pointList(i,2)-pointList(n,2))^2 );
    else
        dist(i) = sqrt( (pointList(i,1)-pointList(i-1,1))^2 + (pointList(i,2)-pointList(i-1,2))^2 );
    end
end
% ȷ���Ƿ���Ҫ��ӽڵ�
pointLabel = dist<=uth;
addNum = length(find(pointLabel==0));
% �����Ҫ��ӽڵ�
if addNum~=0
    % ͨ�����Բ�ֵ��ӵ�
    newPointList = double(zeros((addNum+n),2));
    tmpBasicFunctions = closedCubicSpline(pointList,opt);
    cursor = 0; % ָ���½ڵ��б��е�ǰ�ڵ����һ���ڵ���α�
    for i = 1:n
        if dist(i)>uth
            % �����½ڵ�����
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
            newPoint = (B1 + t.*B2 + (t*t).*B3 + (t^3).*B4)'; % ͨ������������ֵ�����λ�����ڵ㡰���С��ĵ������
            if integerFlag
                newPoint = uint16(newPoint);
            end 
            % ����½ڵ�
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