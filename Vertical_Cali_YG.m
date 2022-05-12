format long; clc; clear
%%vertical calibration Yuming Gu
n = 10;%define the number of height position.
H = 1000; %define the range of height in um.
phi = cell(n,1);
uphi = cell(n,1);

for i = 1:n
    path=['/Users/gyuming/Senior_Design/Photo1/Calibration/',num2str(i),''];
    [phi_i,II]= p_shift(path);
    phi{i,1} = phi_i;
end
uphi_r = itoh_unwrap_c(cell2mat(phi(n,1)));
uphi_r=fitting(uphi_r);
[w,h] = size(uphi_r);
parfor i = 1:n
    U = unwrap_ref(cell2mat(phi(i,1)),uphi_r);
    uphi{i,1} = U;  

end
%%
Kk = zeros(w,h,2);
y = zeros(n,1);
yy = zeros(n,1);
%calibration size
for i = 1
    for j = 1440
        for z = 1:n
            M = cell2mat(uphi(z,1));
            y(z,1) = M(i,j);
        end
        x = linspace(0,9,n)';
        K = polyfit(x,y,1);
        Kk(i,j,1) = K(1);
        Kk(i,j,2) = K(2);
    end 
end
scatter(y,x)


%% the phase shift algorithm
function [M, II]= p_shift(path)
 %3-step phase shifting
%read mage:
fp = fullfile(path, '*.jpg'); % The image format is .jpg
img = dir(fp);
Y = 0;
X = 0;
n = length(img);
for k = 1:n
    fpName = img(k).name;
    imgName = fullfile(path, fpName);
    II = imread(imgName);
     %II = imresize(II,[1000 1000]);
    %II = II(454:2297, 950:2822);
    if size(II,3)==3
        II = rgb2gray(II);
    end 
    II=double(II);
%     II = 2.6522e-5*II.^3 -0.0111*II.^2 + 1.9212*II +1.1274;
%     [w,h] = size(II);
%     Cxy = zeros(w,h);
%     for m = 1:w
%         parfor n = 1:h
%             Cxy(m,n) = f(II(m,n));% need to obtain the Cxy (ratio of I_proj/I_cam)
%         end
%     end
%     I{1,k} = II./Cxy;
    sigma = 2*(k-1)*pi/n;
    Y = Y + sin(sigma) * II;
    X = X + cos(sigma) * II;
end
M = (atan2(-Y,X));
end 
%% 
%Unwrap the image using the Itoh algorithm: the second method is
%performed
%by sequentially unwrapping the all columns, one at a time.
function [Iunwrapped] = itoh_unwrap_c(Iwrapped)
[m n] = size(Iwrapped);
Iunwrapped = Iwrapped;
%%unwrapped columns
for i=1:n
 Iunwrapped(end:-1:1,i) = unwrap(Iunwrapped(end:-1:1,i));
end
% Then sequentially unwrap all the rows
for i=1:m
 Iunwrapped(i,:) = unwrap(Iunwrapped(i,:));
end
end
%%
% referece-guided phaes unwrapping
%%
% referece-guided phaes unwrapping
function phi_u=unwrap_ref(phi,phi_ref)
    phi = phi +2*pi*round((phi_ref-phi)/(2*pi));
    [w,h] = size(phi);
    dphi = zeros(w, h);
    %detrend columns and rows
    for i = 1:w;
        dphi(i, :) = detrend(phi(i, :), 3);
    end
    for j = 1:h
        dphi(:, j) = detrend(phi(:, j), 3);
    end
    phi2 = phi - dphi;
    dphi = F_filter(dphi);
    uphi = dphi + phi2;
    phi_u = uphi;
end
%%
function P = fitting(M)
[w,h] = size(M);
x = 1:h;
x2 = (1:w)';
parfor i = 1:w
A = polyfit(x,M(i,:),2);
y = A(1).*x.^2+A(2).*x.^1+A(3);%.*x+A(4);
M(i,:) = y;
end
parfor j = 1:h
A = polyfit(x2,M(:,j),2);
y = A(1).*x2.^2+A(2).*x2.^1+A(3);%.*x2+A(4);
M(:,j) = y;
end
P = M;
end

%% Fourier filtering
function uphi=F_filter(phi)
[w,h]=size(phi);
% define filted spatial frequency (number of fringe cycle in each direction)
fyy = 32;
fxx = 5; 
% define filter radius
rx = 20; 
ry = 5;
% define centern point
cY=w/2+1;
cX=h/2+1;
% fast 2D fourier transform
FT=fftshift(fft2(phi));
% 7 order filtering
for k = 1:3
    fx = k*fxx; fy = k*fyy;
FT(floor(cY-ry):ceil(cY+ry),floor(cX-fx-rx):ceil(cX-fx+rx))=0;
FT(floor(cY+fy-ry):ceil(cY+fy+ry),floor(cX-rx):ceil(cX+rx))=0;
FT(floor(cY-ry):ceil(cY+ry),floor(cX+fx-rx):ceil(cX+fx+rx))=0;
FT(floor(cY+fy-ry):ceil(cY+fy+ry),floor(cX-rx):ceil(cX+rx))=0;
% D1 = k*fxx-rx;
% mask1= fspecial('disk', D1) == 0;
% mask1 = imresize(padarray(mask1, [floor((x/2)-D1) floor((y/2)-D1)], 1, 'both'), [x y]);
% 
% D2 = k*fxx+rx;
% mask2= fspecial('disk', D2) == 0;
% mask2 = imresize(padarray(mask2, [floor((x/2)-D2) floor((y/2)-D2)], 1, 'both'), [x y]);
% mask = ~(mask1.* ~mask2);
% FT = FT.*mask;
end
% delete high frequency noise
FT(1:ry,floor(cX-rx):ceil(cX+rx))=0;
FT(end-ry:end,floor(cX-rx):ceil(cX+rx))=0;
FT(floor(cY-ry):ceil(cY+ry), 1:rx)=0;
FT(floor(cY-ry):ceil(cY+ry), end-rx:end)=0;
% inverse fourier transform
uphi = real(ifft2(ifftshift(FT)));
end
