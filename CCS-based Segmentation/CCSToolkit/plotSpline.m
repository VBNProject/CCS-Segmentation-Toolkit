function [CCSHan, CCSLeg] = plotSpline(backgroundImage, plotRes, CCSGroup, opt)

try
    cmp = opt.colormap;
    cmpStart = cmp(1,:);
    cmpTermi = cmp(2,:);
    colorFlag = true;
catch
    colorFlag = false;
end

try
    lw = opt.lineWidth;
catch
    lw = 2.0;
end

try
    plotFlag1 = opt.vectorFlag;
catch
    plotFlag1 = false;
end

try
    plotFlag2 = opt.nodeFlag;
catch
    plotFlag2 = false;
end

try
    plotInterval = opt.interval;
catch
    plotInterval = 1;
end

try
    legendFlag = opt.legendFlag;
catch
    legendFlag = true;
end

try
    plotStart = opt.startNum;
catch
    plotStart = 1;
end

try
    noteInt = opt.noteInt;
catch
    noteInt = 10;
end

try
    plotAloneFlag =opt.plotAloneFlag;
catch
    plotAloneFlag = true;
end

try
    saveFlag = opt.saveFlag;
    fileName = opt.fileName;
catch
    saveFlag = false;
end

if plotAloneFlag
    figure;
    imshow(backgroundImage);
end

CCSNum = size(CCSGroup,1);
newNum = floor((CCSNum-1)/plotInterval)+1;

if ((newNum-1)*plotInterval+1)~=CCSNum
    newCCSGroup = cell((newNum+1), 1);
    for i = 1:newNum
        newCCSGroup{i,1} = CCSGroup{((i-1)*plotInterval+1),1};
    end
    newCCSGroup{(newNum+1),1} = CCSGroup{CCSNum,1};
    iterNum = newNum + 1;
else
    newCCSGroup = cell( newNum, 1);
    for i = 1:newNum
        newCCSGroup{i,1} = CCSGroup{((i-1)*plotInterval+1),1};
    end    
    iterNum = newNum;
end

CCSLeg = cell(1,iterNum);
CCSHan = double(zeros(1,iterNum));
for k = 1:iterNum
    basicFunctions = newCCSGroup{k};
    n = size(basicFunctions,1);
    tVec = linspace(0,1,(plotRes+2));
    tVec = tVec';
    B1 = double(zeros(2,1));
    B2 = double(zeros(2,1));
    B3 = double(zeros(2,1));
    B4 = double(zeros(2,1));
    
    if colorFlag
        if iterNum == 1
            color = cmpStart;
        else
            color = cmpStart*double(iterNum-k)/double(iterNum-1) + cmpTermi*double(k-1)/double(iterNum-1);
        end
    else
        color = [rand rand rand];
    end
    
    for i = 1:n
        B1(1,1) = basicFunctions(i,1,1);
        B1(2,1) = basicFunctions(i,1,2);
        B2(1,1) = basicFunctions(i,2,1);
        B2(2,1) = basicFunctions(i,2,2);
        B3(1,1) = basicFunctions(i,3,1);
        B3(2,1) = basicFunctions(i,3,2);
        B4(1,1) = basicFunctions(i,4,1);
        B4(2,1) = basicFunctions(i,4,2);
        
        len = size(tVec,1);
        pointCoord = double(zeros(len,2));
        for j = 1:len
            t = tVec(j);
            pointCoord(j,:) = B1 + t.*B2 + (t*t).*B3 + (t^3).*B4;
        end
        
        hold on;
        if i == 1
            CCSHan(1,k) = plot(pointCoord(:,1),pointCoord(:,2),'Color',color, 'LineWidth', lw);
        else
            plot(pointCoord(:,1),pointCoord(:,2),'Color',color, 'LineWidth', lw);
        end
        
    end

    pointList = extractNodePoints(basicFunctions);
    outFlag = true;
    opt.plotFlag = false;
    normUnitVec = nodeNormal(basicFunctions, outFlag, opt);
    if plotFlag1==true
        hold on;
        xcoord = pointList(:,1);
        ycoord = pointList(:,2);
        xvec = normUnitVec(:,1);
        yvec = normUnitVec(:,2);
        quiver(xcoord, ycoord, xvec, yvec, 0.4, 'Color', color);       
    end
    
    if plotFlag2==true
        hold on;
        xcoord = pointList(:,1);
        ycoord = pointList(:,2);
        plot(xcoord(:,1),ycoord(:,1), 's', 'Color', color);
        for m = 1:noteInt:size(xcoord,1)
            text((xcoord(m,1)),(ycoord(m,1)),num2str(m),'BackgroundColor', [1.0 1.0 1.0], ...
                'EdgeColor', color, 'Color', [0.0 0.0 0.0]);
        end
    end
    
    orgNum = (k-1)*plotInterval+1;
    if orgNum>CCSNum
        orgNum = CCSNum;
    end

    if legendFlag
        try
            CCSLeg{1,k} = ['ID: ' num2str(opt.CCSGroupId) ', SecNum: ' num2str(orgNum+plotStart-1)];
        catch
            CCSLeg{1,k} = ['SecNum: ' num2str(orgNum+plotStart-1)];
        end
    end
end

if plotAloneFlag == true
    if legendFlag
        legend(CCSHan, CCSLeg, 'Location', 'NorthEast');
    end
    
    if saveFlag
        saveas(CCSHan, fileName);
    end
end