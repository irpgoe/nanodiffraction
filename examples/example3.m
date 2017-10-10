clear all; close all; clc
addpath(genpath('toolbox'));

% data storage class
fp = files( 'beamline','id13',...
            'prepath','/home/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_LS2522/DATA/AUTO-TRANSFER/eiger1',...
            'newfile','n12_PD102_sm06_roi7',...
            'detector','eiger',...
            'eigerNr',69);

% experiment class
ex = nanodiffraction('energy',13E3,...
                    'detDistance',0.13852,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1194.49,...
                    'pbz',1214.048); 

% display class
ds = display();

% link classes
link(fp,ex,ex,ds);

% more experimental details
ex.set_scan_info('SNy',101,...
                 'SNz',101);
             
ex.set_roi_and_binning('detectorRoi','off','detRoiY',[1000 1400],'detRoiZ',[1000 1400],'binning','off','biny',2,'binz',2)

% analysis: scan
result1 = ex.analyze_scan('method','stxm','yCrop',[1 5],'zCrop',[1 5]);

% analysis: single frame
frame = fp.read(1);
result2 = ex.analyze_scan('method','b1d','bins',500,'yCrop',[1 1],'zCrop',[1 1]);

% display
figure(1);
subplot(1,3,1);ds.stxm(result1.stxm.df,[]);                     % show stxm scan
subplot(1,3,2);ds.diffraction(frame);                           % show a diffraction pattern
               ds.add_circle(5);                                % show a circle at 5 nm^-1
subplot(1,3,3);ds.saxs(result2.b1d,'limits',[0 5 0.1 100]);     % show 1d saxs curve with custom limits


