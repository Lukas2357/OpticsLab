load('lajolla.mat')



[file,path,indx] = uigetfile('*.tif'); % Here instead of tif you can put any image format you like eg .png, .jpg
if isequal(file,0)
   disp('User selected Cancel')
else
   disp(['User selected ', fullfile(path, file),... 
         ' and filter index: ', num2str(indx)])
end
a=[path,file];
[image3,map]=(imread(a));
D=(imread(a));
%D=D-min(min(D));
s=5;
txt=[num2str(s),'µm'];
imagesc(D)
colorbar
colormap((jet)) %flipud for inverting the colormap
u=max(max(D));
l=min(min(D));

 caxis([l,u])
 
 pbaspect([1 1 1])
