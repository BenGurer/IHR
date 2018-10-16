%% notes
% need change color map - can i do that after? nope - do on overlays
% use function already created after saving figure and change colormap and then save again

% two loops
% loop one
% select hg/ center of cricle
% mrPRint
% anti-alias
% select hg
% check
% if yes save - if no do again
% load next hemi
% load next subject

% loop two
% load coords
% loop through overlays
% change overlay colormap
% mrPrint
% set figure size
% anti-alias
% mask by cirlce
% anti-alias
% crop
% anti-alias
% save as png 300 dpi - black as transparent
function [thisView,coords,radius] = print_parammap(thisView,coords,radius,saveName,filetype,spotlighted)
%      usage: print_paramsmap_spotlight(thisView,radius)
%         by: Ben Gurer
%       date: 11/09/2018
%    purpose: print mrViews current parameter map. Save 'spotlight'
%               of parameter map - the centre defined by user.

%     inputs: view, centre of spotlgiht coords, radius of spotlight

% A4 paper sizes
width = 11.69;     % Width in inches
height = 8.27;    % Height in inches

done = 0;

% spotlight size
% radius = 750;
% coords = [];
while done == 0
    % print mrFigure
    mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
    % set size to A4
    set(gcf,'PaperUnits', 'inches','PaperType', 'a4', 'Position', [0 0 width*100, height*100]);
    
    set(gcf, 'Position', [0 0 width*100, height*100]);
    
    disp('running anti-aliasing...')
    myaa('publish');  % Render an anti-aliased version
    set(gcf,'PaperUnits', 'inches','PaperType', 'a4', 'Position', [0 0 width*100, height*100]);
    disp('Done.')
    
    I = getimage(gcf);
    imshow(I);
    set(gcf,'name','Spotlighted image', 'PaperUnits', 'inches','PaperType', 'a4', 'Position', [0 0 width*100, height*100]);
    
    %     figure('name','Spotlighted image', 'Position', [150 1100 width*100, height*100],...
    %         'PaperType', 'a4','PaperUnits', 'inches');
    if spotlighted == 0
        disp('select centre of spotlight')
        coords = ginput(1);
    end
    imageSize = size(I);
    ci = [coords(2), coords(1), radius];     % center and radius of circle ([c_row, c_col, r])
    [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
    mask = double(uint8((xx.^2 + yy.^2)<ci(3)^2));
    spotlightImage = zeros(size(I));
%     spotlightImage(:,:,1) = I(:,:,1).*mask;
%     spotlightImage(:,:,2) = I(:,:,2).*mask;
%     spotlightImage(:,:,3) = I(:,:,3).*mask;
 
    background = mask;
    background(mask==0) = 1;    
    background(mask==1) = 0;
        
    spotlightImage(:,:,1) = I(:,:,1) + background;
    spotlightImage(:,:,2) = I(:,:,2) + background;
    spotlightImage(:,:,3) = I(:,:,3) + background;
    
    imshow(spotlightImage); 
    
    fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
    close(fh)
    
    fh = findobj( 'Type', 'Figure', 'Name', 'Spotlighted image' );
    close(fh)    
    
    %% crop image
    % if statement for if defined
    iy = coords(2);
    ix = coords(1);
    
    minx = ix-radius;
    miny = iy-radius;
    maxx = ix+radius-minx;
    maxy = iy+radius-miny;
    
    croppedImage = spotlightImage;
    croppedImage(:,[1:minx],:)=[];
    croppedImage([1:miny],:,:)=[];
    croppedImage(:,[maxx:end],:)=[];
    croppedImage([maxy:end],:,:)=[];
    croppedSize = size(croppedImage);
    
    imshow(croppedImage)
    set(gcf,'name','Cropped image');
    if spotlighted == 0
        %% Ask if user is happy with spotlighting
        temp = questdlg('Are you happy with the spotlight?', ...
            'Are you happy with the spotlight?', 'Yes', 'No', 'Yes');
        done = strcmp(temp, 'Yes');
    else
        done = 1;
    end
end

disp('running anti-aliasing...')
myaa('publish');  % Render an anti-aliased version
disp('Done.')

set(gcf,'name','Exported image'); % change figure name
disp('saving...')
export_fig(saveName, ['-',filetype], '-native', '-a1') % save figure

% close figure
fh = findobj( 'Type', 'Figure', 'Name', 'Exported image' );
close(fh)

disp('Done.')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% dev code
% 
% I = imread('test.png');
% 
% % select centre of HG
% p=ginput(2)
% p1max=max(p(:,1));p2max=max(p(:,2));p1min=min(p(:,1));p2min=min(p(:,2));
% 
% 
% 
% myaa('publish');  % Render an anti-aliased version
% 
% % set figure size to a4
% a4_x = 210; % 8.27
% a4_y = 297; % 11.69
% width = 11.69;     % Width in inches
% height = 8.27;    % Height in inches
% set(gcf,'PaperUnits', 'centimeters');
% set(gcf,'PaperType', 'a4');
% papersize = get(gcf, 'PaperSize');
% pos = get(gcf, 'Position');
% get(gcf,'PaperUnits')
% set(gcf, 'Position', [150 1100 width*100, height*100]); %<- Set size
% 
% % Here we preserve the size of the image when we save it.
% % set(gcf,'InvertHardcopy','on');
% 
% left = (papersize(1)- width)/2;
% bottom = (papersize(2)- height)/2;
% myfiguresize = [left, bottom, width, height];
% set(gcf,'PaperPosition', myfiguresize);
% 
% myaa('publish');  % Render an anti-aliased version
% print('test','-dpng','-r450');
% 
% 
% % GO
% 
% % print mrFigure
% mrPrint(thisView,'useDefault=1','roiSmooth=0','roiLabels=0')
% % set size to A4
% fh = findobj( 'Type', 'Figure', 'Name', 'Print figure' );
% 
% myaa('publish');  % Render an anti-aliased version
% 
% 
% I = getimage(gcf);
% 
% radius = 2.5 ;
% dpi = 450;
% 
% figure('Position', [150 1100 width*100, height*100],...
%     'PaperType', 'a4','PaperUnits', 'inches');
% imshow(I);
% p=ginput(1)
% 
% imageSize = size(I);
% ci = [p(2), p(1), 750];     % center and radius of circle ([c_row, c_col, r])
% [xx,yy] = ndgrid((1:imageSize(1))-ci(1),(1:imageSize(2))-ci(2));
% mask = double(uint8((xx.^2 + yy.^2)<ci(3)^2));
% croppedImage = zeros(size(I));
% croppedImage(:,:,1) = I(:,:,1).*mask;
% croppedImage(:,:,2) = I(:,:,2).*mask;
% croppedImage(:,:,3) = I(:,:,3).*mask;
% imshow(croppedImage);
% 
% iy = p(2);
% ix = p(1);
% minx = ix-r;
% miny = iy-r;
% maxx = ix+r-minx;
% maxy = iy+r-miny;
% r = 750;
% Asize = size(croppedImage);
% A = croppedImage;
% A(:,[1:minx],:)=[];
% A([1:miny],:,:)=[];
% A(:,[maxx:end],:)=[];
% A([maxy:end],:,:)=[];
% figure('Position', [150 1100 Asize(1), Asize(2)],...
%     'PaperType', 'a4','PaperUnits', 'inches', 'Color', 'k');
% imshow(A)
% 
% myaa('publish');  % Render an anti-aliased version
% export_fig test.png -native -a1
% 
% tightfig(gcf)
% 
% myaa('publish');  % Render an anti-aliased version
% % print('test','-dpng','-r450');
% fname = 'test';
% format = 'png';
% resolution = 450;
% resolution = round(resolution * 100 / 2.54);
% 
% imwrite(A, fname, format, 'ResolutionUnit', 'meter', ...
%     'XResolution', resolution, 'YResolution', resolution,'Transparency',[0,0,0],'Background',[0,0,0]);
% 
% 
