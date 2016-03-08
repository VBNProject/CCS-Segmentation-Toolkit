function basicFunctions = label2ClosedCubicSpline(binImg, neighborSize, nodeNum, opt)

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

try
    methodStr = opt.method;
catch
    methodStr = 'gh';
end

if ~isempty(find(methodStr=='g', 1)) % gillesFlag为gilles算法调用的flag
    gillesFlag = 1;
else
    gillesFlag = 0;
end

if ~isempty(find(methodStr=='h', 1)) % harris算法调用的flag
    harrisFlag = 1;
else
    harrisFlag = 0;
end

if ~isempty(find(methodStr=='l', 1)) % log算法调用的flag
    logFlag = 1;
else
    logFlag = 0;
end

try
    cloThr = opt.cloThr;
catch 
    cloThr = 8;
end

binImg = double(medfilt2(binImg,[neighborSize neighborSize]));
bwFixI = im2bw(binImg);

tmpList = regionprops(bwFixI,'FilledImage');
tmpBw = tmpList.FilledImage;
tmpList = regionprops(bwFixI,'BoundingBox');
tmpBi = tmpList.BoundingBox;
patchX = ceil(tmpBi(1,2));
patchY = ceil(tmpBi(1,1));
patchXSize = ceil(tmpBi(1,4));
patchYSize = ceil(tmpBi(1,3));
[imheight, imwidth] = size(bwFixI);
oneLabelImage = false(imheight, imwidth);
oneLabelImage(patchX:(patchX+patchXSize-1),patchY:(patchY+patchYSize-1)) = tmpBw(:,:);
bwFixI = oneLabelImage;
fixE = edge(bwFixI,'canny');

foundFlag = 0;
for i = 1:size(fixE,1)
    if foundFlag == 1
        break;
    end
    for j = 1:size(fixE,2)
        xCoord = i;
        yCoord = j;        
        if fixE(xCoord, yCoord) == 1
            foundFlag = 1;
            break;
        end
    end
end

fixContour = traceContour(fixE, xCoord, yCoord);

xsize = size(binImg,1);
ysize = size(binImg,2);
for i = 1:size(fixContour,1)
    if binImg(fixContour(i,2),fixContour(i,1))==0
        conPointX = fixContour(i,2);
        conPointY = fixContour(i,1);
        xmin = max((conPointX-1),0);
        ymin = max((conPointY-1),0);
        xmax = min((conPointX+1),xsize);
        ymax = min((conPointY+1),ysize);
        
        outPoint = uint16(zeros(1,2));
        outPoint(1,1) = conPointX;
        outPoint(1,2) = conPointY;
        dist = Inf;
        for j = xmin:xmax
            for k = ymin:ymax
                if binImg(j,k) == 255
                    tmpdist=(conPointX-j)^2 + (conPointY-k)^2;
                    if tmpdist < dist
                        dist = tmpdist;
                        outPoint(1,1) = j;
                        outPoint(1,2) = k;
                    end
                end
            end
        end
        fixContour(i,2) = outPoint(1,1);
        fixContour(i,1) = outPoint(1,2);
    end
end

fixNxtContour = double(zeros(size(fixContour,1),2));
fixNxtContour(1:(size(fixContour,1)-1),:) = fixContour(2:end,:);
fixNxtContour(end,:) = fixContour(1,:);
fixDifContour = double(fixContour) - fixNxtContour;
fixAddDif = abs(fixDifContour(:,1)) + abs(fixDifContour(:,2));
fixAddDif(fixAddDif~=0) = 1;
fixContour = fixContour(fixAddDif==1,:);


com_ind = [];


if logFlag
    log_points = kp_log(uint8(binImg));
    log_points = log_points(:,1:2);
    [~, log_ind] = ContourProjection(fixContour,log_points);
    sorted_log_ind = sort(log_ind);
    com_ind = CombineInd(sorted_log_ind, com_ind);
end

if (isempty(com_ind)&&(gillesFlag==0)&&(harrisFlag==0))
    defFlag = 1;
end

if (gillesFlag||defFlag)
    gil_points = kp_gilles(binImg);
    gil_points = gil_points(:,1:2);
    [~, gil_ind] = ContourProjection(fixContour,gil_points);

    sorted_gil_ind = sort(gil_ind);
    com_ind = CombineInd(sorted_gil_ind, com_ind);
end

if (harrisFlag||defFlag)
    har_points = kp_harris(binImg);
    har_points = har_points(:,1:2);
    [~, har_ind] = ContourProjection(fixContour,har_points);
    sorted_har_ind = sort(har_ind);
    com_ind = CombineInd(sorted_har_ind, com_ind);
end

com_len = size(com_ind,1);
com_points = double(zeros(com_len,2));
for i = 1:com_len
    curInd = com_ind(i,1);
    com_points(i,:) = fixContour(curInd,:);
end

com_ind_neighbors = uint16(zeros(com_len,2));
com_ind_neighbors(2:end,1) = com_ind(1:(com_len-1),1);
com_ind_neighbors(1,1) = com_ind(end,1);
com_ind_neighbors(1:(com_len-1),2) = com_ind(2:end,1);
com_ind_neighbors(end,2) = com_ind(1,1);
filtered_ind = uint16(zeros(com_len,1));
filtered_ind_pre = uint16(zeros(com_len,1));
filtered_ind_nxt = uint16(zeros(com_len,1));
filt_acc = 0;
for i = 1:com_len
    curInd = double(com_ind(i,1));
    preInd = double(com_ind_neighbors(i,1));
    nxtInd = double(com_ind_neighbors(i,2));

    if preInd > curInd
        preDist = abs(curInd+size(fixContour,1)-preInd);
    else
        preDist = abs(curInd-preInd);
    end
    if nxtInd < curInd
        nxtDist = abs(nxtInd+size(fixContour,1)-curInd);
    else
        nxtDist = abs(nxtInd-curInd);
    end
    
    if ((preDist<=cloThr)||(nxtDist<=cloThr))
        filt_acc = filt_acc + 1;
        filtered_ind(filt_acc,1) = curInd;
        filtered_ind_pre(filt_acc,1) = preInd;
        filtered_ind_nxt(filt_acc,1) = nxtInd;
    end
end
filtered_ind = filtered_ind(filtered_ind~=0);
filtered_ind_pre = filtered_ind_pre(filtered_ind_pre~=0);
filtered_ind_nxt = filtered_ind_nxt(filtered_ind_nxt~=0);
filt_len = size(filtered_ind,1);

if filt_len ~= 0
    filtered_label = uint16(zeros(filt_len,1));
    label_acc = 0;
    for i = 1:filt_len
        curLabel = filtered_label(i,1);
        curInd = filtered_ind(i,1);

        if curLabel==0
            label_acc = label_acc + 1;
            filtered_label(i,1) = label_acc;
        end
        if i == filt_len
            nxtInd = filtered_ind(1,1);
        else
            nxtInd = filtered_ind((i+1),1);
        end

        nxtDist = abs(double(curInd)-double(nxtInd));
        if nxtDist<=cloThr
            if i == filt_len
                filtered_label(1,1) = label_acc;
            else
                filtered_label((i+1),1) = label_acc;
            end
        end
    end
    
    fstInd = filtered_ind(1,1);
    lstInd = filtered_ind(filt_len,1);

    tmpDist = abs(double(fstInd)+size(fixContour,1)-double(lstInd));
    if tmpDist<=cloThr
        lstLabel = filtered_label(filt_len,1);
        lastLen = size(filtered_label(filtered_label==lstLabel),1);
        filtered_label(filtered_label==label_acc) = 1;
        filtered_label = [filtered_label((filt_len-lastLen+1):end,1); filtered_label(1:(filt_len-lastLen),1)];
        filtered_ind = [filtered_ind((filt_len-lastLen+1):end,1); filtered_ind(1:(filt_len-lastLen),1)];
        filtered_ind_pre = [filtered_ind_pre((filt_len-lastLen+1):end,1); filtered_ind_pre(1:(filt_len-lastLen),1)];
        filtered_ind_nxt = [filtered_ind_nxt((filt_len-lastLen+1):end,1); filtered_ind_nxt(1:(filt_len-lastLen),1)];
        label_acc = label_acc - 1;
    end

    delete_label = true(com_len,1);
    for i = 1:label_acc
        curGroup = find(filtered_label==i);
        curSize = size(curGroup,1);
        if curSize==2
            curNode1 = curGroup(1,1);
            curNode2 = curGroup(2,1);
            curInd1 = filtered_ind(curNode1,1);
            preInd = filtered_ind_pre(curNode1,1);
            curInd2 = filtered_ind(curNode2,1);
            nxtInd = filtered_ind_nxt(curNode2,1);
            if preInd > curInd1
                preDist = abs(double(curInd1)+size(fixContour,1)-double(preInd));
            else
                preDist = abs(double(curInd1)-double(preInd));
            end
            if nxtInd < curInd2
                nxtDist = abs(double(nxtInd)+size(fixContour,1)-double(curInd2));
            else
                nxtDist = abs(double(curInd2)-double(nxtInd));
            end           
            if (preDist<=nxtDist)
                comInd = com_ind==curInd1;
                delete_label(comInd,1) = false;
            else
                comInd = com_ind==curInd2;
                delete_label(comInd,1) = false;
            end
        elseif curSize==3
            midNode = curGroup(2,1);
            midInd = filtered_ind(midNode,1);
            comInd = com_ind == midInd;
            delete_label(comInd,1) = false;
        elseif curSize>=4
            local_delete_label = true(curSize,1);
            distArray = double(zeros((curSize-1),1));
            for k = 1:(curSize-1)
                if filtered_ind(curGroup(k,1),1)>filtered_ind(curGroup((k+1),1),1)
                    distArray(k) = abs(double(filtered_ind(curGroup(k,1),1))+size(fixContour,1)-double(filtered_ind(curGroup((k+1),1),1)));                    
                else
                    distArray(k) = abs(double(filtered_ind(curGroup(k,1),1))-double(filtered_ind(curGroup((k+1),1),1)));
                end
            end
            
            while (min(distArray)<=cloThr)&&(size(distArray,1)>1)
                deletedList = find(local_delete_label==false);
                [~, minind] = min(distArray);
                tmpDistArray = double(zeros((size(distArray,1)-1),1));
                if minind == 1
                    delInd = 2;
                    tmpDistArray(1,:) = distArray(1,:) + distArray(2,:);
                    tmpDistArray(2:end,:) = distArray(3:end,:);
                elseif minind == size(distArray,1)
                    delInd = size(distArray,1);
                    tmpDistArray(1:(size(tmpDistArray,1)-1),:) = distArray(1:(size(distArray,1)-2),:);
                    tmpDistArray(end,:) = distArray((size(distArray,1)-1),:) + distArray(end,:);
                else
                    preDist = distArray((minind-1),1);
                    nxtDist = distArray((minind+1),1);
                    if preDist < nxtDist
                        delInd = minind;
                        if minind > 2
                            tmpDistArray(1:(minind-2),:) = distArray(1:(minind-2),:);
                        end
                        tmpDistArray((minind-1),:) = distArray((minind-1),:) + distArray(minind,:);
                        tmpDistArray(minind:end,:) = distArray((minind+1):end,:);
                    else
                        delInd = minind+1;
                        tmpDistArray(1:(minind-1),:) = distArray(1:(minind-1),:);
                        tmpDistArray(minind,:) = distArray(minind,:) + distArray((minind+1),:);
                        if minind < (size(distArray,1)-1)
                            tmpDistArray((minind+1):end,:) = distArray((minind+2):end,:);
                        end
                    end
                end
                
                if isempty(deletedList)
                    local_delete_label(delInd,1) = false;
                else
                    segArray = double(zeros((size(deletedList,1)+1),1));
                    for k = 1:size(segArray,1)
                        if k == 1
                            segArray(k) = deletedList(k)-1-1;
                        elseif k == size(segArray,1)
                            segArray(k) = size(local_delete_label,1) - deletedList((k-1),1)-1;
                        else
                            segArray(k) = deletedList(k,1) - deletedList((k-1),1) - 1;
                        end
                    end

                    delInd = delInd - 1; 
                    for k = 1:size(segArray,1)
                        if delInd<=segArray(k,1)
                            if k == 1
                                delInd = delInd + 1;
                            else
                                delInd = delInd + deletedList((k-1),1);
                            end
                            break;
                        else
                            delInd = delInd - segArray(k,1);
                        end
                    end
                    local_delete_label(delInd,1) = false;
                end
                distArray = tmpDistArray;
            end
            
            startInd = curGroup(1,1);
            deletedList = find(local_delete_label==false);
            deletedList = deletedList + startInd-1;
            for k = 1:size(deletedList,1)
                absInd = filtered_ind(deletedList(k,1),1);
                comInd = com_ind==absInd;
                delete_label(comInd,1) = false;          
            end
        end
    end
    com_ind = com_ind(delete_label);
    com_len = size(com_ind,1);
end

curva = LineCurvature2D(fixContour); 
curva(isnan(curva))=0;
curva = abs(curva);

newInd = com_ind;
if nodeNum > com_len
    cmpNum = nodeNum - com_len;
    for i = 1:cmpNum
        newInd_len = size(newInd,1);
        newIndDist = uint16(zeros(newInd_len,1));
        for j = 1:newInd_len
            if j == newInd_len
                newIndDist(j,1) = abs(double(newInd(1,1)) + size(fixContour,1) -double(newInd(end,1)));
            else
                newIndDist(j,1) = abs(double(newInd((j+1),1)) - double(newInd(j,1)));
            end
        end
        
        [maxDist,maxIndex] = max(newIndDist);
        if maxDist <= 2*(cloThr+1)
            break;
        end
        
        tmpInd = uint16(zeros((newInd_len+1),1));
        if maxIndex~=newInd_len
            startVal = newInd(maxIndex,1);
            termiVal = newInd((maxIndex+1),1);
            curCurva = curva((startVal+cloThr+1):(termiVal-cloThr-1));
            [~,mIndVal] = max(curCurva);
            insertedVal = mIndVal + startVal + cloThr;
            tmpInd(1:maxIndex,:) = newInd(1:maxIndex,:);
            tmpInd((maxIndex+1),:) = insertedVal;
            tmpInd((maxIndex+2):end,:) = newInd((maxIndex+1):end,:);
        else
            startVal = newInd(maxIndex,1);
            termiVal = newInd(1,1);
            curCurva = [curva(startVal:end); curva(1:termiVal)];
            tmpLen = size(curCurva,1);
            curCurva = curCurva((cloThr+1):(tmpLen-cloThr));
            [~,mIndVal] = max(curCurva);
            mIndVal = mIndVal + cloThr;
            lstLen = size(fixContour,1)-startVal+1;
            if mIndVal<=lstLen
                insertedVal = mIndVal + startVal - 1;
                tmpInd(1:newInd_len,:) = newInd(:,:);
                tmpInd((newInd_len+1),:) = insertedVal;
            else
                insertedVal = mIndVal + startVal - size(fixContour,1) - 1;
                tmpInd(1,:) = insertedVal;
                tmpInd(2:(newInd_len+1),:) = newInd(:,:);
            end
        end
        newInd = tmpInd;
    end
elseif nodeNum < com_len
    cmpNum = com_len - nodeNum;
    for i = 1:cmpNum
        newInd_len = size(newInd,1);
        newIndDist = uint16(zeros(newInd_len,1));
        for j = 1:newInd_len
            if j == newInd_len
                newIndDist(j,1) = abs(double(newInd(1,1)) + size(fixContour,1) -double(newInd(end,1)));
            else
                newIndDist(j,1) = abs(double(newInd((j+1),1)) - double(newInd(j,1)));
            end
        end
        sidIndDist = uint16(zeros(newInd_len,1));
        for j = 1:newInd_len
            if j == 1
                sidIndDist(j,1) = newIndDist(newInd_len,1) + newIndDist(1,1);
            else
                sidIndDist(j,1) = newIndDist((j-1),1) + newIndDist(j,1);
            end
        end
        [~,minIndex] = min(sidIndDist);
        tmpInd = uint16(zeros((newInd_len-1),1));
        if minIndex ~= 1
            tmpInd(1:(minIndex-1),:) = newInd(1:(minIndex-1),:);
            tmpInd(minIndex:end,:) = newInd((minIndex+1):end,:);
        else
            tmpInd(:,:) = newInd(2:end,:);
        end
        newInd = tmpInd;
    end
end

com_len = size(newInd,1);
com_points = double(zeros(com_len,2));
for i = 1:com_len
    curInd = newInd(i,1);
    com_points(i,:) = fixContour(curInd,:);
end

lopt.plotFlag = false;
if plotFlag == true
    figure;
    imshow(bwFixI);
    len = size(com_points,1);
    for i = 1:len
        hold on;
        xcoord = com_points(i,1);
        ycoord = com_points(i,2);
        plot(xcoord,ycoord,'bo');
    end
    lopt.plotFlag = true;
end

basicFunctions = closedCubicSpline(com_points, lopt);

end

function tracedContour = traceContour(fixE, xCoord, yCoord)
try
    fixContour = bwtraceboundary(fixE, [double(xCoord), double(yCoord)], 'W', 8);
catch
    try
        fixContour = bwtraceboundary(fixE, [double(xCoord), double(yCoord)], 'N', 8);
    catch
        try
            fixContour = bwtraceboundary(fixE, [double(xCoord), double(yCoord)], 'S', 8);
        catch
            try
                fixContour = bwtraceboundary(fixE, [double(xCoord), double(yCoord)], 'E', 8);
            end
        end
    end
end

tracedContour = [fixContour(:,2) fixContour(:,1)];
end

function comIndList = CombineInd(orgIndList, newIndList)
comIndList = [orgIndList; newIndList];
comIndList = unique(comIndList);
comIndList = sort(comIndList);
end

function [proPointList, proIndList] = ContourProjection(contourList, orgPointList)
plen = size(orgPointList,1);
proPointList = double(zeros(plen,2));
proIndList = double(zeros(plen,1));
clen = size(contourList,1);
for i = 1:plen
    curP = orgPointList(i,:);
    curP = double(curP);
    repP = repmat(curP,clen,1);
    repP = repP - double(contourList);
    sqrtP = (repP(:,1).^2 + repP(:,2).^2).^0.5;
    [~,ind] = min(sqrtP);
    proP = contourList(ind,:);
    proPointList(i,:) = proP;
    proIndList(i,1) = ind;
end

proPointList = unique(proPointList,'rows');
proIndList = unique(proIndList,'rows');
end

function k=LineCurvature2D(Vertices,Lines)
% Function is written by D.Kroon University of Twente (August 2011)

if(nargin<2)
    Lines=[(1:(size(Vertices,1)-1))' (2:size(Vertices,1))'];
end

Na=zeros(size(Vertices,1),1); Nb=zeros(size(Vertices,1),1);
Na(Lines(:,1))=Lines(:,2); Nb(Lines(:,2))=Lines(:,1);

checkNa=Na==0; checkNb=Nb==0;
Naa=Na; Nbb=Nb;
Naa(checkNa)=find(checkNa); Nbb(checkNb)=find(checkNb);

Na(checkNa)=Nbb(Nbb(checkNa)); Nb(checkNb)=Naa(Naa(checkNb));

Ta=-sqrt(sum((Vertices-Vertices(Na,:)).^2,2));
Tb=sqrt(sum((Vertices-Vertices(Nb,:)).^2,2)); 

Ta(checkNa)=-Ta(checkNa); Tb(checkNb)=-Tb(checkNb);

x = [Vertices(Na,1) Vertices(:,1) Vertices(Nb,1)];
y = [Vertices(Na,2) Vertices(:,2) Vertices(Nb,2)];
M = [ones(size(Tb)) -Ta Ta.^2 ones(size(Tb)) zeros(size(Tb)) zeros(size(Tb)) ones(size(Tb)) -Tb Tb.^2];
invM=inverse3(M);
a(:,1)=invM(:,1,1).*x(:,1)+invM(:,2,1).*x(:,2)+invM(:,3,1).*x(:,3);
a(:,2)=invM(:,1,2).*x(:,1)+invM(:,2,2).*x(:,2)+invM(:,3,2).*x(:,3);
a(:,3)=invM(:,1,3).*x(:,1)+invM(:,2,3).*x(:,2)+invM(:,3,3).*x(:,3);
b(:,1)=invM(:,1,1).*y(:,1)+invM(:,2,1).*y(:,2)+invM(:,3,1).*y(:,3);
b(:,2)=invM(:,1,2).*y(:,1)+invM(:,2,2).*y(:,2)+invM(:,3,2).*y(:,3);
b(:,3)=invM(:,1,3).*y(:,1)+invM(:,2,3).*y(:,2)+invM(:,3,3).*y(:,3);

k = 2*(a(:,2).*b(:,3)-a(:,3).*b(:,2)) ./ ((a(:,2).^2+b(:,2).^2).^(3/2));
end

function  Minv  = inverse3(M)
adjM(:,1,1)=  M(:,5).*M(:,9)-M(:,8).*M(:,6);
adjM(:,1,2)=  -(M(:,4).*M(:,9)-M(:,7).*M(:,6));
adjM(:,1,3)=  M(:,4).*M(:,8)-M(:,7).*M(:,5);
adjM(:,2,1)=  -(M(:,2).*M(:,9)-M(:,8).*M(:,3));
adjM(:,2,2)=  M(:,1).*M(:,9)-M(:,7).*M(:,3);
adjM(:,2,3)=  -(M(:,1).*M(:,8)-M(:,7).*M(:,2));
adjM(:,3,1)=  M(:,2).*M(:,6)-M(:,5).*M(:,3);
adjM(:,3,2)=  -(M(:,1).*M(:,6)-M(:,4).*M(:,3));
adjM(:,3,3)=  M(:,1).*M(:,5)-M(:,4).*M(:,2);
detM=M(:,1).*M(:,5).*M(:,9)-M(:,1).*M(:,8).*M(:,6)-M(:,4).*M(:,2).*M(:,9)+M(:,4).*M(:,8).*M(:,3)+M(:,7).*M(:,2).*M(:,6)-M(:,7).*M(:,5).*M(:,3);
Minv=bsxfun(@rdivide,adjM,detM);
end

function points = kp_gilles(im,o_radius)
% Extract keypoints using Gilles algorithm
% Author :: Vincent Garcia
% Date   :: 05/12/2007

im = im(:,:,1);
if nargin==1
    radius = 10;
else
    radius = o_radius;
end
mask = fspecial('disk',radius)>0;
loc_ent = entropyfilt(im,mask);
[~,~,tmp] = findLocalMaximum(loc_ent,radius);

[l,c]     = find(tmp>0.95*max(tmp(:)));
points    = [c,l,repmat(radius,[size(l,1),1])];
end

function points = kp_harris(im)
% Extract keypoints using Harris algorithm (with an improvement version)
% Author :: Vincent Garcia
% Date   :: 05/12/2007

im = double(im(:,:,1));
sigma = 1.5;

s_D = 0.7*sigma;
x  = -round(3*s_D):round(3*s_D);
dx = x .* exp(-x.*x/(2*s_D*s_D)) ./ (s_D*s_D*s_D*sqrt(2*pi));
dy = dx';

Ix = conv2(im, dx, 'same');
Iy = conv2(im, dy, 'same');

s_I = sigma;
g = fspecial('gaussian',max(1,fix(6*s_I+1)), s_I);
Ix2 = conv2(Ix.^2, g, 'same');
Iy2 = conv2(Iy.^2, g, 'same');
Ixy = conv2(Ix.*Iy, g, 'same');

cim = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps);
[~,~,max_local] = findLocalMaximum(cim,3*s_I);
t = 0.01*max(max_local(:));
[r,c] = find(max_local>=t);
points = [c,r];
end

function [points] = kp_log(img,o_nb_blobs)
% Extract keypoints using Laplacian of Gaussian (LoG) algorithm
% Author :: Vincent Garcia
% Date   :: 05/12/2007

img = double(img(:,:,1));

if nargin==1
    nb_blobs = 120;
else
    nb_blobs = o_nb_blobs;
end

sigma_begin = 2;
sigma_end   = 15;
sigma_step  = 1;
sigma_array = sigma_begin:sigma_step:sigma_end;
sigma_nb    = numel(sigma_array);

img_height  = size(img,1);
img_width   = size(img,2);

snlo = zeros(img_height,img_width,sigma_nb);
for i=1:sigma_nb
    sigma       = sigma_array(i);
    snlo(:,:,i) = sigma*sigma*imfilter(img,fspecial('log', floor(6*sigma+1), sigma),'replicate');
end

snlo_dil             = imdilate(snlo,ones(3,3,3));
blob_candidate_index = find(snlo==snlo_dil);
blob_candidate_value = snlo(blob_candidate_index);
[~,index]          = sort(blob_candidate_value,'descend');
blob_index           = blob_candidate_index( index(1:min(nb_blobs,numel(index))) );
[lig,col,sca]        = ind2sub([img_height,img_width,sigma_nb],blob_index);
points               = [col,lig,3*reshape(sigma_array(sca),[size(lig,1),1])];
end

function [row,col,max_local] = findLocalMaximum(val,radius)
% Determine the local maximum of a given value
% Author :: Vincent Garcia
% Date   :: 09/02/2007

mask  = fspecial('disk',radius)>0;
nb    = sum(mask(:));
highest          = ordfilt2(val, nb, mask);
second_highest   = ordfilt2(val, nb-1, mask);
index            = highest==val & highest~=second_highest;
max_local        = zeros(size(val));
max_local(index) = val(index);
[row,col]        = find(index==1);
end