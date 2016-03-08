function [global_pos, arr_out, local_out_img] = NormSearch2(croppedImage, croppedImage2, invM, a)
% 用途：沿矢量方向寻找最优解
% 输入参数
% croppedImage：原始图像（初始轮廓所在的图像，pic1）中，某个CCS节点周围沿节点矢量方向截取的局部区域
% croppedImage2：原始图像的下一张图像（即要演化出新轮廓的图像，pic2）中，某个CCS节点周围沿节点矢量方向截取的局部区域
% invM：从pic2的截取区域croppedImage2的平面坐标系，逆映射回pic2的平面坐标系的线性变换矩阵。
% a：控制内外力平衡的参数，取值范围[0 1.0]，越接近1，外力作用越大
% 输出参数
% global_pos：节点演化到的新位置在pic2的坐标系下的坐标值，double类型
% arr_out：根据用户输入的选项，输出滑动窗口沿矢量方向滑动到不同位置时，得到的欧氏距离、夹角余弦、协方差、pearson相关系数和相关系数
% 分别为arr_out.arr_euc, arr_out.arr_cos, arr_out.arr_cov, arr_out.arr_pea, arr_out.arr_cor
% 每个arr_xxx的大小均为1*arrSize

%% 处理基本参数
height = size(croppedImage,1);
width = size(croppedImage,2);
arrSize = size(croppedImage2,1) - size(croppedImage,1) + 1;
orgProfile = double(reshape(croppedImage, height*width, 1));

%% 在截取区域内使用滑动窗口沿矢量方向寻找最优匹配位置
arr_euc = double(zeros(arrSize,1));
arr_cos = double(zeros(arrSize,1));
arr_cov = double(zeros(arrSize,1));
arr_pea = double(zeros(arrSize,1));
arr_cor = double(zeros(arrSize,1));

for i = 1:arrSize
    tmpImage = croppedImage2(i:(i+height-1),:);
    cmpProfile = double(reshape(tmpImage, height*width, 1));
    % 欧氏距离
    arr_euc(i,1) = 1-pdist([orgProfile';cmpProfile'],'euclidean');
    % 夹角余弦
    arr_cos(i,1) = 1-pdist([orgProfile';cmpProfile'],'cosine');
    % 协方差
    arr_cov(i,1) = mean((double(orgProfile) - mean(double(orgProfile))).*(double(cmpProfile) - mean(double(cmpProfile))));
    % pearson相关系数
    arr_pea(i,1) = arr_cov(i,1)/(std(double(orgProfile))*std(double(cmpProfile)));
    % 相关系数（实际上就是互相关系数）
    arr_cor(i,1) = 1 - pdist([orgProfile';cmpProfile'], 'correlation');
    % figure; imshow(croppedImage);
    % figure; plot(orgProfile,'r'); hold on; plot(cmpProfile,'b');
end

min_arr_euc = min(arr_euc);
max_arr_euc = max(arr_euc);
min_arr_cov = min(arr_cov);
max_arr_cov = max(arr_cov);

for i = 1:arrSize
    arr_euc(i) = (arr_euc(i)-min_arr_euc)/(max_arr_euc-min_arr_euc);
    arr_cov(i) = (arr_cov(i)-min_arr_cov)/(max_arr_cov-min_arr_cov);
end

[~,ind_euc] = max(arr_euc);
[~,ind_cos] = max(arr_cos);
[~,ind_cov] = max(arr_cov);
[~,ind_pea] = max(arr_pea);
[~,ind_cor] = max(arr_cor);

% 绘制滑动窗口的相关性曲线
%     figure; plot(arr_euc,'y');
%     hold on;plot(arr_cos,'r');
%     hold on;plot(arr_cov,'k');
%     hold on;plot(arr_pea,'g');
%     hold on;plot(arr_cor,'b');

% deopt.Image = true; deopt.Profile = true; deopt.DifPro = true; pos = 16;
% DebugPlotProfile(croppedImage, croppedImage2, pos, 10, deopt);

in_pos = median([ind_euc ind_cos ind_cov ind_pea ind_cor]); % 获取5个相关性指标计算的最佳位置的中位数
in_pos = in_pos+height/2-1; % 从arr_size的index转换为图像的纵坐标

%% 计算外力将节点拉伸到的位置，该位置仅于croppedImage2有关
meanDifArray = double(zeros(arrSize,1));
for i = 1:arrSize
    if mod(height,2) == 0 % 如果为偶数
        halfH = height/2;
        inImg = croppedImage2(1:(i+halfH-1), :);
        outImg = croppedImage2((i+halfH):end, :);
    else % 如果为奇数
        halfH = (height-1)/2;
        inImg = croppedImage2(1:(i+halfH-1), :);
        outImg = croppedImage2((i+halfH+1):end, :);
    end        
    meanDifArray(i,1) = abs(mean(mean(inImg))-mean(mean(outImg)));
end
[~, maxMeanInd] = max(meanDifArray);
if mod(height,2) == 0
    out_pos = halfH + maxMeanInd - 0.5;
else
    out_pos = halfH + maxMeanInd;
end

%% 计算节点演化到的新位置
% 控制内外力平衡的参数
% a：外力权重，a越大，外力起的作用越大，取值[0 1]
try
    local_a = a;
catch
    local_a = 0.5;
end
bal_pos = local_a*out_pos + (1-local_a)*in_pos;

%% 计算节点演化到的新位置
% ind_avg = median([ind_euc ind_cos ind_cov ind_pea ind_cor]); % 获取5个相关性指标计算的最佳位置的中位数
% local_pos = [double(width/2), double(ind_avg+height/2-1), 1.0]';
local_pos = [double(width/2), double(bal_pos), 1.0]';
local_out_img = croppedImage2;
local_out_img(uint16(bal_pos), :) = 255;
global_pos = invM*local_pos;

arr_out.arr_euc = arr_euc';
arr_out.arr_cos = arr_cos';
arr_out.arr_cov = arr_cov';
arr_out.arr_pea = arr_pea';
arr_out.arr_cor = arr_cor';