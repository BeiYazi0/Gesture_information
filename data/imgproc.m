path1='./data_test';
path2='./data_test2';
gesture=['\left     ';'\down     ';'\up       ';'\right    ';...
    '\leftright';'\updown   ';'\true     ';'\false    '];
for i=1:8
    source=[path1,strip(gesture(i,:))];
    dst=[path2,strip(gesture(i,:))];
    for j=1:10
        id=['/',num2str(j)];
        improc(source,dst,id);
    end
end

function improc(source,dst,id)
img1=imread([source,id,'/range.png']);
img2=imread([source,id,'/doppler.png']);
img3=imread([source,id,'/angle.png']);
img=uint8(zeros(size(img1)));
img(:,:,1)=sum(img1,3)/3;
img(:,:,2)=sum(img2,3)/3;
img(:,:,3)=sum(img3,3)/3;
imwrite(img,[dst,id,'.png']);
end