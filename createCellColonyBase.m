bufsize = 200;
outBuffer = uint8(zeros(bufsize,bufsize,bufsize));
ncRad = bufsize*0.3;

for k = 1:bufsize
    for j = 1:bufsize
        for i = 1:bufsize
            r = sqrt((i-bufsize/2)*(i-bufsize/2)+(j-bufsize/2)*(j-bufsize/2)+(k-bufsize/2)*(k-bufsize/2));
            if r<=ncRad
                outBuffer(i,j,k) = 255;
            else
                outBuffer(i,j,k) = 0;
            end
        end
    end
    imwrite(outBuffer(:,:,k),['./cellColony1_base/cell_' num2str(k,'%03d') '.tif']);
end