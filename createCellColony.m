% ����ϸ���ۼ���ģ��ͼƬ
% ģ�������У�ϸ��Ϊ���ɰ뾶��һ��С�򣬶�����Ϊһ���ض��뾶����������
% �ں��������ڲ���ϸ�����Ը��ߵĸ��ʳ��֣����ں����ⲿ��ϸ�����ֵ��ܶ���ϵ�

%% cellColony1ģ�Ͳ���
% bufsize = 200; % ͼ��xyz��������ĳߴ�
% cellGV = [140 255]; % ϸ���Ҷ�ֵ��Χ
% cellRad = [2 3]; % ϸ���İ뾶��Χ
% ncProb = 0.02; % ������ϸ�����ֵĸ���
% bgProb = 0.00001; % ������ϸ�����ֵĸ���
% ncRad = bufsize*0.3; % ��������İ뾶
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % �������������
% dstDir = './cellColony1/cell_';
% maskDstDir = './cellColony1_mask/cell_';

%% cellColony2ģ�Ͳ���
% bufsize = 200; % ͼ��xyz��������ĳߴ�
% cellGV = [140 255]; % ϸ���Ҷ�ֵ��Χ
% cellRad = [2 3]; % ϸ���İ뾶��Χ
% ncProb = 0.02; % ������ϸ�����ֵĸ���
% bgProb = 0.001; % ������ϸ�����ֵĸ���
% ncRad = bufsize*0.3; % ��������İ뾶
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % �������������
% dstDir = './cellColony2/cell_';
% maskDstDir = './cellColony2_mask/cell_';

%% cellColony3ģ�Ͳ���
% bufsize = 200; % ͼ��xyz��������ĳߴ�
% cellGV = [140 255]; % ϸ���Ҷ�ֵ��Χ
% cellRad = [2 3]; % ϸ���İ뾶��Χ
% ncProb = 0.007; % ������ϸ�����ֵĸ���
% bgProb = 0.0005; % ������ϸ�����ֵĸ���
% ncRad = bufsize*0.3; % ��������İ뾶
% ncCen = [bufsize/2, bufsize/2, bufsize/2]; % �������������
% dstDir = './cellColony3/cell_';
% maskDstDir = './cellColony3_mask/cell_';

%% cellColony4ģ�Ͳ���
bufsize = 200; % ͼ��xyz��������ĳߴ�
cellGV = [140 255]; % ϸ���Ҷ�ֵ��Χ
cellRad = [2 3]; % ϸ���İ뾶��Χ
ncProb = 0.007; % ������ϸ�����ֵĸ���
bgProb = 0.001; % ������ϸ�����ֵĸ���
ncRad = bufsize*0.3; % ��������İ뾶
ncCen = [bufsize/2, bufsize/2, bufsize/2]; % �������������
dstDir = './cellColony4/cell_';
maskDstDir = './cellColony4_mask/cell_';

%% ����ģ��
outBuffer = uint8(zeros(bufsize,bufsize,bufsize));
maskBuffer = uint8(zeros(bufsize,bufsize,bufsize));
for k = 1:bufsize % z
    for j = 1:bufsize % y
        for i = 1:bufsize % x
            r = sqrt((i-ncCen(1))*(i-ncCen(1))+(j-ncCen(2))*(j-ncCen(2))+(k-ncCen(3))*(k-ncCen(3)));
            cellrand = rand();
            % ���λ�ں��������ڣ������ֵ������(0 0.4)֮��  ���� ���λ�ں��������⣬�����ֵ������(0 bgProb)֮��
            putCheck = ((r<=ncRad)&&(cellrand<ncProb))||((r>ncRad)&&(cellrand<bgProb));
            if putCheck == true   
                % ��ȡһ����cellGV����������ֲ��ĻҶ�ֵ
                cgv = uint8(cellGV(1) + rand()*(cellGV(2)-cellGV(1)+1)); 
                % ��ȡһ������ֲ��İ뾶ֵ
                crd = int16(cellRad(1) + rand()*(cellRad(2)-cellRad(1)+1));
                % �ڸð뾶�����ǰ��õ��ĻҶ�
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
