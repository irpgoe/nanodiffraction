clear all; close all; clc
addpath(genpath('toolbox'));

%% data storage class
fp = files();

% example pattern
ff = fp.read(1);

%% experiment class
ex = nanodiffraction('energy',13.8E3,...
                    'detDistance',5.1,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1600,...
                    'pbz',1600); 

% display class
ds = display();
link(ex,ds);

ds.diffraction(ex.qr .* ~(ff == max(ff(:))));
ds.add_circle(1.047,'circleColor','red');

print(gcf,'ex3.png','-dpng');

%%
figure(1);
subplot(1,3,1);ds.stxm(result1.stxm.df,[]);                     % show stxm scan
subplot(1,3,2);ds.diffraction(frame);                           % show a diffraction pattern
               ds.add_circle(5);                                % show a circle at 5 nm^-1
subplot(1,3,3);ds.saxs(result2.b1d,'limits',[0 5 0.1 100]);     % show 1d saxs curve with custom limits


