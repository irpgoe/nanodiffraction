% Analysis of run sc4304 > ID13 > ESRF > 2016 
clear all; close all; clc

% are you on santa or on dancer?
addpath(genpath('toolbox'));
prepath = which_computer('santa');

% make path definitions, add custom directories
fp = files( 'beamline','id13',...
            'prepath',[prepath 'Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1'],...
            'newfile','herz2_roi2',...
            'detector','eiger',...
            'eigerNr',201);   
        
% experiments
e = nanodiffraction('energy',14.6407E3,...
                    'detDistance',1.928650,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1500.499,...
                    'pbz',1372.509);          
         
% setup display
ds = display();

% link classes
link(fp,e,e,ds);
                
% pca init: Note that pca init should be done before binning is activated
% or pixels should be given in new values!
e.pca_init('r1',86,'r2',230);

% set detector mask
e.set_roi_and_binning('binning','on','biny',4,'binz',4,'detectorRoi','on','detRoiY',round([(e.p.pby-200) (e.p.pby+200)]),'detRoiZ',round([(e.p.pbz-200) (e.p.pbz+200)]));
e.set_mask('detector',fp.read(1)==4294967295);
e.set_mask('pca',e.p.detMaskPCA);

% scan parameters
e.set_scan_info('SNy',101,...
                'SNz',101,...
                'stepy',0.5e-6,...
                'stepz',0.5e-6);

%% analyze scan
res = e.analyze_scan('method','stxm+pca');

%% show result
close all
f1 = figure(1);
ds.stxm(res.stxm.df,struct('scale',0.02,'sampl',10,'unit','mm'));
caxis([5000 10000]);
set(imhandles(f1),'AlphaData',res.stxm.df>2e3);

f2 = figure(2);
ds.stxm(res.pca.angle,struct('scale',0.02,'sampl',10,'unit','mm'));
hold on
ds.pca(res.pca.angle);
set(imhandles(f2),'AlphaData',res.stxm.df>2e3);

f3 = figure(3);
ds.stxm(res.pca.angle,struct('scale',0.02,'sampl',10,'unit','mm'));
hold on
ds.pca(res.pca.angle,struct('quiver','on'));
set(imhandles(f3),'AlphaData',res.stxm.df>2e3);

%% clicktool
e.clicktool(ds,[])
