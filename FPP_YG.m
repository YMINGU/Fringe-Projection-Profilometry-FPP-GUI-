clear;clc;format long;
z = cell(3, 1);
delta = cell(3,1);
%[f_5_rc,gof_right] = fit(coeff_img_int_r',proj_inten','smoothingspline');
for i = 1:3%1:5
%path of target images
path = ['/Users/gyuming/Senior_Design/Photo1/',''];
%path of reference images
%path2 = ['/Users/haolinzhang/Documents/FPP_2021/12_07_2021/ref',num2str(i),''];
path2 = ['/Users/gyuming/Senior_Design/Photo1/Ref/',''];
%%
% filePattern = fullfile(path, '*.Bmp');
% imgFiles = dir(filePattern);
% N = length(imgFiles);
% pattern_num = N/3;
% phi_t_all = cell(pattern_num,1);

%%
phi_t = p_shift(path); %need to obtain f from sensorMode
phi_r = p_shift(path2);
uphi_rr = unfold(phi_r);
uphi_r=fitting(uphi_rr);
uphi_t = unwrap_ref(phi_t,uphi_r);
delta{i,1} =F_filter(uphi_t-uphi_r);
z{i,1} = cell2mat(delta(i,1)).*0.8151; % obtain height matrix (need Kk from the vertical calibration)
height_c_rxy = cell2mat(z(i,1));
aheight=abs(height_c_rxy);
% surf(y,x,z,'EdgeColor','interp'); colorbar;
% xlabel('x-direction (mm)')
% ylabel('y-direction (mm)')
end

mesh(aheight)

%% the phase shift algorithm
function M = p_shift(path)
 %3-step phase shifting
%read mage:
%%
fp = fullfile(path, '*.jpg'); % The image format is .jpg
img = dir(fp);
Y = 0;
X = 0;
Y1 = 0
X1 = 0;
n = 3
pattern = length(img)/3;
%if pattern == 1
for k = 1:3
    fpName = img(k).name;
    imgName = fullfile(path, fpName);
    II = imread(imgName);
     %II = imresize(II,[1000 1000]);
    %II = II(454:2297, 950:2822);
    if size(II,3)==3
        II = rgb2gray(II);
    end 
    II=double(II);
    [w,h] = size(II);
%     for i = 1:w
%         parfor j =1:h
%             II(i,j) = f(II(i,j));
%         end 
%     end 
%     imshow(uint8(II));
    sigma = 2*(k-1)*pi/n;
    Y = Y + sin(sigma) * II;
    X = X + cos(sigma) * II;
end
M = (atan2(-Y,X));
end 

% if pattern == 2
% for k = 1:3
%     fpName = img(k).name;
%     imgName = fullfile(path, fpName);
%     II = imread(imgName);
%      %II = imresize(II,[1000 1000]);
%     %II = II(454:2297, 950:2822);
%     if size(II,3)==3
%         II = rgb2gray(II);
%     end 
%     II=double(II);
% %     [w,h] = size(II);
% %     for i = 1:w
% %         parfor j =1:h
% %             II(i,j) = f(II(i,j));
% %         end 
% %     end 
% %     imshow(uint8(II));
%     sigma = 2*(k-1)*pi/n;
%     Y = Y + sin(sigma) * II;
%     X = X + cos(sigma) * II;
% end
% M1 = (atan2(-Y,X));
% for k = 4:6
%     fpName = img(k).name;
%     imgName = fullfile(path, fpName);
%     II = imread(imgName);
%      %II = imresize(II,[1000 1000]);
%     %II = II(454:2297, 950:2822);
%     if size(II,3)==3
%         II = rgb2gray(II);
%     end 
%     II=double(II);
% %     [w,h] = size(II);
% %     for i = 1:w
% %         parfor j =1:h
% %             II(i,j) = f(II(i,j));
% %         end 
% %     end 
% %     imshow(uint8(II));
%     sigma = 2*(k-3-1)*pi/n;
%     Y1 = Y1 + sin(sigma) * II;
%     X1 = X1 + cos(sigma) * II;
% end 
% M2 = (atan2(-Y1,X1));
% M = (M1+M2)/2;
% end 
%end 

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
%% unfold phase jump of reference plane
function uphi_ref = unfold(phi)
[w h] = size(phi);
uphi_ref = phi;
% unwrapped columns
    for i=1:h
        uphi_ref(:,i) = unwrap(uphi_ref(:,i));
    end
% unwrap rows
    for j=1:w
        uphi_ref(j,:) = unwrap(uphi_ref(j,:));
    end
end
%%
% referece-guided phaes unwrapping
function phi_u=unwrap_ref(phi,phi_ref)
%fourier
    phi = phi +2*pi*round((phi_ref-phi)/(2*pi));
    [w,h] = size(phi);
    phi2 = fitting(phi);
    dphi = phi-phi2;
    dphi = F_filter(dphi);
    uphi = dphi + phi2;
    phi_u = uphi;


%reference_guided
%     phi2=phi+pi;
%     phi2 = phi2 +2*pi*round((phi_ref-phi2)/2/pi);
%     phi_u=phi2;

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
%fxx = 57; 
fxx= 5;
% define filter radius
%rx = 20;
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
