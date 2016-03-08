function [global_pos, arr_out, local_out_img] = NormSearch2(croppedImage, croppedImage2, invM, a)
% ��;����ʸ������Ѱ�����Ž�
% �������
% croppedImage��ԭʼͼ�񣨳�ʼ�������ڵ�ͼ��pic1���У�ĳ��CCS�ڵ���Χ�ؽڵ�ʸ�������ȡ�ľֲ�����
% croppedImage2��ԭʼͼ�����һ��ͼ�񣨼�Ҫ�ݻ�����������ͼ��pic2���У�ĳ��CCS�ڵ���Χ�ؽڵ�ʸ�������ȡ�ľֲ�����
% invM����pic2�Ľ�ȡ����croppedImage2��ƽ������ϵ����ӳ���pic2��ƽ������ϵ�����Ա任����
% a������������ƽ��Ĳ�����ȡֵ��Χ[0 1.0]��Խ�ӽ�1����������Խ��
% �������
% global_pos���ڵ��ݻ�������λ����pic2������ϵ�µ�����ֵ��double����
% arr_out�������û������ѡ��������������ʸ�����򻬶�����ͬλ��ʱ���õ���ŷ�Ͼ��롢�н����ҡ�Э���pearson���ϵ�������ϵ��
% �ֱ�Ϊarr_out.arr_euc, arr_out.arr_cos, arr_out.arr_cov, arr_out.arr_pea, arr_out.arr_cor
% ÿ��arr_xxx�Ĵ�С��Ϊ1*arrSize

%% �����������
height = size(croppedImage,1);
width = size(croppedImage,2);
arrSize = size(croppedImage2,1) - size(croppedImage,1) + 1;
orgProfile = double(reshape(croppedImage, height*width, 1));

%% �ڽ�ȡ������ʹ�û���������ʸ������Ѱ������ƥ��λ��
arr_euc = double(zeros(arrSize,1));
arr_cos = double(zeros(arrSize,1));
arr_cov = double(zeros(arrSize,1));
arr_pea = double(zeros(arrSize,1));
arr_cor = double(zeros(arrSize,1));

for i = 1:arrSize
    tmpImage = croppedImage2(i:(i+height-1),:);
    cmpProfile = double(reshape(tmpImage, height*width, 1));
    % ŷ�Ͼ���
    arr_euc(i,1) = 1-pdist([orgProfile';cmpProfile'],'euclidean');
    % �н�����
    arr_cos(i,1) = 1-pdist([orgProfile';cmpProfile'],'cosine');
    % Э����
    arr_cov(i,1) = mean((double(orgProfile) - mean(double(orgProfile))).*(double(cmpProfile) - mean(double(cmpProfile))));
    % pearson���ϵ��
    arr_pea(i,1) = arr_cov(i,1)/(std(double(orgProfile))*std(double(cmpProfile)));
    % ���ϵ����ʵ���Ͼ��ǻ����ϵ����
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

% ���ƻ������ڵ����������
%     figure; plot(arr_euc,'y');
%     hold on;plot(arr_cos,'r');
%     hold on;plot(arr_cov,'k');
%     hold on;plot(arr_pea,'g');
%     hold on;plot(arr_cor,'b');

% deopt.Image = true; deopt.Profile = true; deopt.DifPro = true; pos = 16;
% DebugPlotProfile(croppedImage, croppedImage2, pos, 10, deopt);

in_pos = median([ind_euc ind_cos ind_cov ind_pea ind_cor]); % ��ȡ5�������ָ���������λ�õ���λ��
in_pos = in_pos+height/2-1; % ��arr_size��indexת��Ϊͼ���������

%% �����������ڵ����쵽��λ�ã���λ�ý���croppedImage2�й�
meanDifArray = double(zeros(arrSize,1));
for i = 1:arrSize
    if mod(height,2) == 0 % ���Ϊż��
        halfH = height/2;
        inImg = croppedImage2(1:(i+halfH-1), :);
        outImg = croppedImage2((i+halfH):end, :);
    else % ���Ϊ����
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

%% ����ڵ��ݻ�������λ��
% ����������ƽ��Ĳ���
% a������Ȩ�أ�aԽ�������������Խ��ȡֵ[0 1]
try
    local_a = a;
catch
    local_a = 0.5;
end
bal_pos = local_a*out_pos + (1-local_a)*in_pos;

%% ����ڵ��ݻ�������λ��
% ind_avg = median([ind_euc ind_cos ind_cov ind_pea ind_cor]); % ��ȡ5�������ָ���������λ�õ���λ��
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