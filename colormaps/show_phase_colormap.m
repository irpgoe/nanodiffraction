
% specify a grid and radial mask
e = nanodiffraction('Ny',512,'Nz',512,'pby',256,'pbz',256);
m1 = e.radial_mask('r1',-1,'r2',256,'grid',e.R);
m2 = e.radial_mask('r1',100,'r2',256,'grid',e.R);


%% cmap
fig = figure(1);
fig.Position = [100 100 1000 1000];
r=3;c=3;
ax1=subplot(r,c,1);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax2=subplot(r,c,2);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax3=subplot(r,c,3);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax4=subplot(r,c,4);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax5=subplot(r,c,5);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax6=subplot(r,c,6);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax7=subplot(r,c,7);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax8=subplot(r,c,8);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax9=subplot(r,c,9);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

colormap(ax1,cmap('C1'));
colormap(ax2,cmap('C2'));
colormap(ax3,cmap('C3'));
colormap(ax4,cmap('C4'));
colormap(ax5,cmap('C5'));
colormap(ax6,cmap('C6'));
colormap(ax7,cmap('C7'));
colormap(ax8,cmap('C8'));
colormap(ax9,cmap('C9'));

print(gcf,'peterkovesi.png','-dpng');

%% phasemap

fig = figure(2);
fig.Position = [100 100 330 330];
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);
colormap(phasemap);

print(gcf,'phasemap.png','-dpng');

%% colorblind

fig = figure(3);
fig.Position = [100 100 1000 330];
Z = peaks(512);
r=1;c=3;
ax1=subplot(r,c,1);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax2=subplot(r,c,2);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax3=subplot(r,c,3);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

colormap(ax1,ametrine);
colormap(ax2,isolum);
colormap(ax3,morgenstemning);

print(gcf,'colorblind.png','-dpng');

%% pmkmp

fig = figure(4);
fig.Position = [100 100 1000 1000];
Z = peaks(512);
r=3;c=3;
ax1=subplot(r,c,1);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax2=subplot(r,c,2);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax3=subplot(r,c,3);
h=imagesc(e.phi);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m2);

ax4=subplot(r,c,4);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax5=subplot(r,c,5);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax6=subplot(r,c,6);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax7=subplot(r,c,7);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax8=subplot(r,c,8);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

ax9=subplot(r,c,9);
h=imagesc(Z);axis image;box off;cb=colorbar;set(cb,'YTick',[])
set(gca,'Visible','off')
set(h,'AlphaData',m1);

colormap(ax1,pmkmp(256,'IsoL'));
colormap(ax2,pmkmp(256,'IsoAZ'));
colormap(ax3,pmkmp(256,'IsoAZ180'));
colormap(ax4,pmkmp(256,'LinearL'));
colormap(ax5,pmkmp(256,'LinLhot'));
colormap(ax6,pmkmp(256,'CubicYF'));
colormap(ax7,pmkmp(256,'CubicL'));
colormap(ax8,pmkmp(256,'Swtth'));
colormap(ax9,pmkmp(256,'Edge'));

print(gcf,'pmkmp.png','-dpng');


%% preferred linear maps
bgrLevel = 20;
variance = 0.1;
meanValue = 0;

x = linspace(-5,5,128);
y = linspace(-5,5,128);
[X,Y] = meshgrid(x,y);
Z = Y.*sin(X) - X.*cos(Y) + bgrLevel;
DF = Z + sqrt(variance)*randn(size(Z)) + meanValue;


d = display();
fig = figure(5);
fig.Position = [100 100 1000 800];
ax1=subplot(2,2,1);
    imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    title('Matlab Parula colormap');
    
ax2=subplot(2,2,2);
    imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    title('Colormap LinearL by Matteo Niccoli');
    
ax3=subplot(2,2,3);
    imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    title('Colormap CubicYF by Matteo Niccoli');
    
ax4=subplot(2,2,4);
    imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);   
    title('Colormap Morgenstemning by Matthias Geissbuehler');

colormap(ax1,parula);    
colormap(ax2,pmkmp(256,'LinearL'));
colormap(ax3,pmkmp(256,'CubicYF'));    
colormap(ax4,morgenstemning);

print(gcf,'preferred_linear_colormaps.png','-dpng');


%% preferred linear with alpha

fig = figure(6);
fig.Position = [100 100 1000 400];
subplot(1,2,1);
    imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);   
    title('No masking');
    
subplot(1,2,2);
    h=imagesc(DF);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);   
    title('Masking of background using alpha level');
    set(h,'AlphaData',DF>16);
colormap(morgenstemning);

print(gcf,'masking_alpha_level.png','-dpng');

%% preferred cyclic maps
x = linspace(-1.5,2.5,128);           % coordinates
y = linspace(-2,2,128);              % -2:2 is the central part
[X,Y] = meshgrid(x,y);
Z = (Y.^2 + (X - 2).^2).^(-1/2) - ... % force field (electric field between
    (Y.^2 + (X - 1/2).^2).^(-1/2)/2;  % two charges
[FX,FY] = gradient(Z);                % gradient

O = atan2d(-FY,FX);                   % orientation (note that in a figure 
                                      % the y-direction goes 'down'
O(O<0) = O(O<0) + 180;                % project onto range [0 180)

d = display();
fig = figure(7);
fig.Position = [100 100 1000 800];
ax1=subplot(2,2,1);
    imagesc(O);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    caxis([0 180]);
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));
    title('Matlab HSV colormap');
    
ax2=subplot(2,2,2);
    imagesc(O);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    caxis([0 180]);
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));
    title('Colormap IsoAZ by Matteo Niccoli');
    
ax3=subplot(2,2,3);
    imagesc(O);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    caxis([0 180]);
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));
    title('Cyclic colormap C2 by Peter Kovesi');
    
ax4=subplot(2,2,4);
    imagesc(O);axis image;cb=colorbar;set(cb,'YTick',[])
    set(gca,'XTick',[],'YTick',[]);
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    caxis([0 180]);
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));        
    title('Cyclic colormap C8 by Peter Kovesi');

colormap(ax1,hsv);    
colormap(ax2,pmkmp(256,'IsoAZ'));
colormap(ax3,cmap('C2'));    
colormap(ax4,cmap('C8'));

print(gcf,'preferred_cyclic_colormaps.png','-dpng');