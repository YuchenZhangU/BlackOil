clear all 
close all
NX = 60;
NY = 60;

load field.out
phi =reshape( field(:,1),NX,NY);
kx = reshape( field(:,2),NX,NY);
ky = reshape( field(:,3),NX,NY);
kz = reshape( field(:,4),NX,NY);
figure;
imagesc(phi),colormap('jet');
axis square
figure;
imagesc(kx),colormap('jet');
axis square
figure;
imagesc(ky),colormap('jet');
figure;
imagesc(kz),colormap('jet');
axis square
clear field

load results.out
p = reshape( results(:,1),NX,NY);
sw = reshape( results(:,2),NX,NY);
sg = reshape( results(:,3),NX,NY);
so = reshape( results(:,4),NX,NY);
Rs = reshape( results(:,5),NX,NY);
clear results
figure;
imagesc(p),colormap('jet');
axis square
colorbar
figure;
imagesc(sw),colormap('jet');
axis square
colorbar
figure;
imagesc(sg),colormap('jet');
axis square
colorbar
figure;
imagesc(so),colormap('jet');
axis square
colorbar
