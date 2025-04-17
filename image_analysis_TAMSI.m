clc
clear all
close all

% Read in timeseries of images
    ND2file = 'Pathname';
    data = bfopen(ND2file);
    nT = length(data{1}); % count how many images (time points)
    [nR,nC] = size(data{1}{1}); % Grab image size (rows, columns)
    dt = 5/nT; % change dt if timelapse changes from 5 sec

% Look at first images 
figure(1);
for t = 1:50    
    IM01 = imcomplement( imadjust(data{1}{t}) );
    imshow(IM01, 'InitialMagnification', 'fit'); % Resize image to see it better 
    
    drawnow
    pause(0.05)
end

%-----------------------------------------

for t = 1:nT       
        IM01 = double( imcomplement( imadjust(data{1}{t}) ) );
        IM01 = IM01./max(IM01(:)); % normalize image
        
        IM02{t} = IM01;
end

% Plot changes to pixels over time -----------------------------------------------
        
[init_pixelchange,s] = pixelchange(IM02,dt);

init_pixelchange_sum = sum(init_pixelchange,3);

figure(2)
x = 1:size(init_pixelchange_sum,2); 
y = 1:size(init_pixelchange_sum,1); 
heatmap = imagesc(x,y,init_pixelchange_sum, [0 80]); 
h = colorbar;
axis equal
    
% Set polygon ROI
    roi = drawpolygon; % Draw polygon
    ROIxy = roi.Position;
    BW = poly2mask(ROIxy(:,1),ROIxy(:,2),nR,nC); % Convert polygon to binary mask

for t = 1:nT       
        IM011 = double( imcomplement( imadjust(data{1}{t}) ) );
        IM011 = IM011./max(IM011(:)); % normalize image
        IM011(~BW) = 0; % Apply regional mask to FV boundary
        IM022{t} = IM011;
end

[final_pixelchange,s] = pixelchange(IM022,dt);
    
final_pixelchange_sum = sum(final_pixelchange,3);
final_pixelchange_sumindiv = squeeze(sum(final_pixelchange, [1 2]));

nonzeros = nnz(final_pixelchange_sum); % determine how many pixels are within FV boundary


y1 = final_pixelchange_sumindiv./nonzeros;
x1(:,1) = s(1,1:end)';

% fit power law function to change in pixel intensity^2/elapsed time
    fitfun = fittype('C*( 1 - 1/( ((x1/t0).^0.5) +1 ))', 'independent', 'x1', 'dependent', 'y1');
    fit_line = fit(x1,y1,fitfun, 'Start', [0.05 0.00005]);
    coeffvals(1,:) = coeffvalues(fit_line);

% plot change in pixel intensity^2/elapsed time with power-law fit
figure(3)
plot(x1,y1,'Linewidth',8,'Color','#a559aa')
hold on
fit_semi = plot(fit_line, ':');
set(fit_semi,'Linewidth',3,'Color', '#000000')

TAMSI = [x1, y1]; % column 1 = elapsed time, column 2 = avg pixel change at elapsed time
Dist = abs(0.18-x1);
[M,I] = min(Dist); 

pixelchange2_per_second = TAMSI(I,2)/(I*dt);  % find value closest to 0.18 s elapsed time 
pixelchange_per_second_arb = pixelchange2_per_second*100; % compute velocity and arbitrarily multiply by factor of 100 
