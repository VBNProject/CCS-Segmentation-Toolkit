function intersectPoints = raySplineIntersection(basicFunctions, raysParam, opt)

n = size(basicFunctions,1);
len = size(raysParam,1);
intersectPoints = cell(len,1);

try
    plotFlag = opt.plotFlag;
catch
    plotFlag = false;
end

B1 = double(zeros(2,1));
xcoord = double(zeros(n,1));
ycoord = double(zeros(n,1));

for i = 1:n
    B1(1,1) = basicFunctions(i,1,1);
    B1(2,1) = basicFunctions(i,1,2);
    if i == 1
        xcoord(n,1) = B1(1,1);
        ycoord(n,1) = B1(2,1);
    else
        xcoord((i-1),1) = B1(1,1);
        ycoord((i-1),1) = B1(2,1);
    end
end

for i = 1:len
    xorg = raysParam(i,1);
    yorg = raysParam(i,2);
    xvec = raysParam(i,3);
    yvec = raysParam(i,4);
    
    interPoints = [];
    if xvec~=0
        k = yvec/xvec;
        b = yorg - k*xorg;
 
        for j = 1:n
            bonex = basicFunctions(j,1,1);
            boney = basicFunctions(j,1,2);
            btwox = basicFunctions(j,2,1);
            btwoy = basicFunctions(j,2,2);
            bthrx = basicFunctions(j,3,1);
            bthry = basicFunctions(j,3,2);
            bfoux = basicFunctions(j,4,1);
            bfouy = basicFunctions(j,4,2);

            % kone + ktwo*t + kthr*t^2 + kfou*t^3 = 0
            kone = k*bonex-boney+b;
            ktwo = k*btwox-btwoy;
            kthr = k*bthrx-bthry;
            kfou = k*bfoux-bfouy;
            fun = @(t)kone+ktwo*t+kthr*(t^2)+kfou*(t^3);
            [xarray,fxarray,rootnum] = multiRoot(fun,0,1);
            if (isempty(xarray))||(isempty(fxarray))
                [xarray,fxarray,rootnum] = multiRoot(fun,0,1, 1000);
            end
            if (~isempty(xarray))&&(~isempty(fxarray))
                for kl = 1:rootnum
                    tsol = xarray(kl,1);
                    xcoord = bonex+btwox*tsol+bthrx*(tsol^2)+bfoux*(tsol^3);
                    ycoord = boney+btwoy*tsol+bthry*(tsol^2)+bfouy*(tsol^3);
                    interVec = [(xcoord-xorg) (ycoord-yorg)];
                    rayVec = [xvec yvec];
                    tmpProduct = dot(interVec,rayVec);
                    if tmpProduct >= 0
                        interPoints = [interPoints;[xcoord ycoord]];
                    end
                end
            end           
        end
        
    else

        for j = 1:n
            bonex = basicFunctions(j,1,1);
            boney = basicFunctions(j,1,2);
            btwox = basicFunctions(j,2,1);
            btwoy = basicFunctions(j,2,2);
            bthrx = basicFunctions(j,3,1);
            bthry = basicFunctions(j,3,2);
            bfoux = basicFunctions(j,4,1);
            bfouy = basicFunctions(j,4,2);

            fun = @(t)bonex+btwox*t+bthrx*(t^2)+bfoux*(t^3)-xorg;
            [xarray,fxarray,rootnum] = multiRoot(fun,0,1);
            if (isempty(xarray))||(isempty(fxarray))
                [xarray,fxarray,rootnum] = multiRoot(fun,0,1, 1000);
            end
            if (~isempty(xarray))&&(~isempty(fxarray))
                for kl = 1:rootnum
                    tsol = xarray(kl,1);
                    xcoord = bonex+btwox*tsol+bthrx*(tsol^2)+bfoux*(tsol^3);
                    ycoord = boney+btwoy*tsol+bthry*(tsol^2)+bfouy*(tsol^3);
                    interVec = [(xcoord-xorg) (ycoord-yorg)];
                    rayVec = [xvec yvec];
                    tmpProduct = dot(interVec,rayVec);
                    if tmpProduct >= 0
                        interPoints = [interPoints;[xcoord ycoord]];
                    end
                end
            end
       end
        
    end

    interPoints = unique(interPoints, 'rows');
    plen = size(interPoints,1);
    sameArray = true(plen,1);
    if plen ~= 1 
        eps = 1e-3;
        for m = 1:(plen-1)
            curPoint = interPoints(m,:);
            for nn = (m+1):plen
                tarPoint = interPoints(nn,:);
                sameCheck = (abs(tarPoint(1,1)-curPoint(1,1))<eps)&&(abs(tarPoint(1,2)-curPoint(1,2))<eps);
                if sameCheck
                    sameArray(nn,1) = false;
                end
                
            end
        end
    end
    interPoints = interPoints(sameArray,:);

    intersectPoints{i,1} = interPoints;
end

if plotFlag == true
    for i = 1:len
        interPoints = intersectPoints{i,1};
        len = size(interPoints,1);
        for j = 1:len
            hold on;
            xcoord = interPoints(j,1);
            ycoord = interPoints(j,2);
            plot(xcoord,ycoord,'go');
        end
    end
end

end

function [x,fx,iter] = bisect(fun,a,b,eps,varargin)
    if nargin<3
        error('Not enough input parameters');
    end
    
    if (nargin<4)||(isempty(eps))
        eps = 1e-6;
    end

    fa = feval(fun,a,varargin{:});
    fb = feval(fun,b,varargin{:});
    
    k=1;
    if fa*fb>0
        x = NaN;
        fx = NaN;
    elseif fa==0
        x = a;
        fx = fa;
    elseif fb==0
        x = b;
        fx = fb;
    else
        while abs(b-a)>eps
            x = (a+b)/2;
            fx = feval(fun,x,varargin{:});
            if fa*fx>0
                a = x;
                fa = fx;
            elseif fb*fx>0
                b = x;
                fb = fx;
            else
                break;
            end
            k=k+1;
        end
    end
    iter = k;
end

function [x,fx,rootnum] = multiRoot(fun,a,b,n,eps,varargin)
    if (nargin<4)||(isempty(n))
        n = 100;
    end
    if (nargin<5)||(isempty(eps))
        eps = 1e-6;
    end
    
    x = [];
    fx = [];
    rootnum = 0;
    res = 1.0/double(n);
    for i = a:res:b
        ai = i;
        bi = i+res;
        fai = feval(fun,ai,varargin{:});
        fbi = feval(fun,bi,varargin{:});
        subMul = fai*fbi;
        if subMul<=0
            [xi,fxi,~] = bisect(fun,ai,bi,eps);
            x = [x; xi];
            fx = [fx; fxi];
            rootnum = rootnum + 1;
        end
    end
    
end