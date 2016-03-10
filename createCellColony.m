% 创建细胞聚集的模型图片
% 模型数据中，细胞为若干半径不一的小球，而核团为一个特定半径的球形区域
% 在核团区域内部，细胞会以更高的概率出现，而在核团外部，细胞出现的密度则较低

%% cellColony1模型参数
% bufsize = 200; % 图像xyz三个方向的尺寸
% cellGV = [140 255]; % 细胞灰度值范围
% cellRad = [2 3]; % 细胞的半径范围
% ncProb = 0.02; % 核团内细胞出现的概率
% bgProb = 0.00001; % 核团外细胞出现的概率
% ncRad = bufsize*0.3; % 核团区域的半径
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % 核团区域的中心
% dstDir = './cellColony1/cell_';
% maskDstDir = './cellColony1_mask/cell_';

%% cellColony2模型参数
% bufsize = 200; % 图像xyz三个方向的尺寸
% cellGV = [140 255]; % 细胞灰度值范围
% cellRad = [2 3]; % 细胞的半径范围
% ncProb = 0.02; % 核团内细胞出现的概率
% bgProb = 0.001; % 核团外细胞出现的概率
% ncRad = bufsize*0.3; % 核团区域的半径
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % 核团区域的中心
% dstDir = './cellColony2/cell_';
% maskDstDir = './cellColony2_mask/cell_';

%% cellColony3模型参数
% bufsize = 200; % 图像xyz三个方向的尺寸
% cellGV = [140 255]; % 细胞灰度值范围
% cellRad = [2 3]; % 细胞的半径范围
% ncProb = 0.007; % 核团内细胞出现的概率
% bgProb = 0.0005; % 核团外细胞出现的概率
% ncRad = bufsize*0.3; % 核团区域的半径
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % 核团区域的中心
% dstDir = './cellColony3/cell_';
% maskDstDir = './cellColony3_mask/cell_';

%% cellColony4模型参数
bufsize = 200; % 图像xyz三个方向的尺寸
cellGV = [140 255]; % 细胞灰度值范围
cellRad = [2 3]; % 细胞的半径范围
ncProb = 0.007; % 核团内细胞出现的概率
bgProb = 0.001; % 核团外细胞出现的概率
ncRad = bufsize*0.3; % 核团区域的半径
ncCen = [bufsize/2, bufsize/2, bufsize/2]; % 核团区域的中心
dstDir = './cellColony4/cell_';
maskDstDir = './cellColony4_mask/cell_';

%% 构造模型
outBuffer = uint8(zeros(bufsize,bufsize,bufsize));
maskBuffer = uint8(zeros(bufsize,bufsize,bufsize));
for k = 1:bufsize % z
    for j = 1:bufsize % y
        for i = 1:bufsize % x
            r = sqrt((i-ncCen(1))*(i-ncCen(1))+(j-ncCen(2))*(j-ncCen(2))+(k-ncCen(3))*(k-ncCen(3)));
            cellrand = rand();
            % 如果位于核团区域内，且随机值落在了(0 0.4)之间  或者 如果位于核团区域外，且随机值落在了(0 bgProb)之间
            putCheck = ((r<=ncRad)&&(cellrand<ncProb))||((r>ncRad)&&(cellrand<bgProb));
            if putCheck == true   
                % 获取一个在cellGV区间内随机分布的灰度值
                cgv = uint8(cellGV(1) + rand()*(cellGV(2)-cellGV(1)+1)); 
                % 获取一个随机分布的半径值
                crd = int16(cellRad(1) + rand()*(cellRad(2)-cellRad(1)+1));
                % 在该半径内填充前面得到的灰度
                for m = -crd:1:crd % z
                    for n = -crd:1:crd % y
                        for p = -crd:1:crd % x
                            % cr = sqrt(double((abs(m)-crd)*(abs(m)-crd)+(abs(n)-crd)*(abs(n)-crd)+(abs(p)-crd)*(abs(p)-crd)));
                            cr = sqrt(double((abs(m))*(abs(m))+(abs(n))*(abs(n))+(abs(p))*(abs(p))));
                            if cr<=crd
                                currentx = i + p;
                                currenty = j + n;
                                currentz = k + m;
                                inCheck = (currentx>0)&&(currentx<=bufsize)&&(currenty>0)&&(currenty<=bufsize)&&(currentz>0)&&(currentz<=bufsize);
                                if inCheck == true
                                    outBuffer(currentx, currenty, currentz) = max(outBuffer(currentx, currenty, currentz),cgv);
                                    maskBuffer(currentx, currenty, currentz) = 255;
                                end                                
                            end
                            
                        end
                    end
                end
                
            end
        end
    end
    disp(k);
end

for k = 1:bufsize
    imwrite(outBuffer(:,:,k),[dstDir num2str(k,'%03d') '.tif']);
    imwrite(maskBuffer(:,:,k),[maskDstDir num2str(k,'%03d') '.tif']);
end
