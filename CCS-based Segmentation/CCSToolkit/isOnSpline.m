function [onSplineFlag,segNum,t] = isOnSpline(basicFunctions, targetPoints, opt)

n = size(basicFunctions,1);
len = size(targetPoints,1);
onSplineFlag = false(len,1);
segNum = NaN(len,1);
t = NaN(len,1);

try
    eps = opt.biEps;
catch
    eps = 1e-6;
end

try
    eqs = opt.eqEps;
catch
    eqs = 0.01;
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
    px = targetPoints(i,1);
    py = targetPoints(i,2);

    [onFlag, segN, tVal] = verifyOnSpline(px, py, xcoord, ycoord, basicFunctions, eps, eqs, true);
    if onFlag == false
        [onFlag, segN, tVal] = verifyOnSpline(px, py, xcoord, ycoord, basicFunctions, eps, eqs, false);
    end
    
    onSplineFlag(i,1) = onFlag;
    segNum(i,1)=segN;
    t(i,1) = tVal;
end

end

function [onFlag, segN, tVal] = verifyOnSpline(px, py, xcoord, ycoord, basicFunctions, eps, eqs, fastSearch)
    n = size(basicFunctions,1);
    foundFlag = false;
    
    for j = 1:n      
        
        if j == 1
            insideFlag = ((px>=xcoord(n))&&(px<(xcoord(1))))||((px>xcoord(1))&&(px<=(xcoord(n))))||((py>=ycoord(n))&&(py<(ycoord(1))))||((py>ycoord(1))&&(py<=(ycoord(n))));
        else
            insideFlag = ((px>=xcoord(j-1))&&(px<(xcoord(j))))||((px>xcoord(j))&&(px<=(xcoord(j-1))))||((py>=ycoord(j-1))&&(py<(ycoord(j))))||((py>ycoord(j))&&(py<=(ycoord(j-1))));
        end
        
        if fastSearch == false
            insideFlag = true;
        end

        if insideFlag
            bonex = basicFunctions(j,1,1);
            boney = basicFunctions(j,1,2);
            btwox = basicFunctions(j,2,1);
            btwoy = basicFunctions(j,2,2);
            bthrx = basicFunctions(j,3,1);
            bthry = basicFunctions(j,3,2);
            bfoux = basicFunctions(j,4,1);
            bfouy = basicFunctions(j,4,2);
            
            xfun = @(t)bonex+btwox*t+bthrx*(t^2)+bfoux*(t^3)-px;
            yfun = @(t)boney+btwoy*t+bthry*(t^2)+bfouy*(t^3)-py;
            
            [xarray,~,~] = multiRoot(xfun,0,1,100,eps);
            [yarray,~,~] = multiRoot(yfun,0,1,100,eps);
            sizexa = size(xarray,1);
            sizeya = size(yarray,1);
            
            tarray = double(zeros(max(sizexa,sizeya),2));
            tacc = 0;
            if (~isempty(xarray))&&(~isempty(yarray))
                for p = 1:sizexa
                    for q = 1:sizeya
                        if abs(xarray(p)-yarray(q))<eqs
                            tacc = tacc + 1;
                            tarray(tacc,1) = (xarray(p)+yarray(q))/2.0;
                            tarray(tacc,2) = abs(xarray(p)-yarray(q));
                        end
                    end
                end
            end
            
            if tacc ~=0
                foundFlag = true;
                onFlag = true;
                segN = j;
                tarray = tarray(1:tacc,:);
                [~, minInd] = min(tarray(:,2));
                tVal = tarray(minInd,1);
            end

        end
        
    end
    
    if foundFlag == false
        onFlag = false;
        segN = NaN;
        tVal = NaN;
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