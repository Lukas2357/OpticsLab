% Read spe files from the Princeton EMCCD. 
[file,path,indx] = uigetfile('*.spe');
if isequal(file,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(path, file),... 
         ' and filter index: ', num2str(indx)])
end
a=[path,file];
I=readSPE(a);
I=double(I);
I=I-min(min(I));
figure
imagesc(I)
colorbar
colormap('jet')
shading interp
clc

pbaspect([1 1 1])
