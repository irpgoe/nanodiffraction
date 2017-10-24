% Example analysis script based on the data set published in Nicolas et
% al. (2017). Scanning X-ray Diffraction on Cardiac Tissue: Automatized 
% Data Analysis and Processing, J. Synchrotron Rad. 24.

% add the toolbox to the path
addpath(genpath('nanodiffraction')); % edit

% file module
f = files( 'beamline','id13',...
           'prepath','/homegroups/AG_Salditt/Projects_cellular_diffraction_and_actin/Analysis/test',... % edit
           'newfile','herz2_roi2',...
           'detector','eiger',...
           'scan',201);   
        
% experiment module
e = nanodiffraction('energy',14.6407E3,...
                    'detDistance',1.928650,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1500.499,...
                    'pbz',1372.509);          
         
% display module
d = display();

% link modules (first f with e and then e with d)
link(f,e,e,d);

% masks
e.set_mask(f.read(1)==(2^32-1));         
e.set_selection(e.radial_mask('r1',86,'r2',230));

% roi and binning
e.set_roi_and_binning('binning','on','biny',4,'binz',4,'detectorRoi','on','detRoiY',round([(e.pby_orig-200) (e.pby_orig+200)]),'detRoiZ',round([(e.pbz_orig-200) (e.pbz_orig+200)]));

% scan parameters
e.set_scan_info('SNy',101,...
                'SNz',101,...
                'stepy',0.5e-6,...
                'stepz',0.5e-6);

%% analyze scan / methods can be combined using the '+' notation
% as methods: stxm | pca | crystal | symmetry | average | heal | sum can be
% used
res = e.analyze_scan('method','stxm+pca');

%% show result
close all
f1 = figure(1);
d.stxm(res.stxm.df,struct('scale',0.02,'sampl',10,'unit','mm'));
caxis([5000 10000]);
set(imhandles(f1),'AlphaData',res.stxm.df>2e3);

f2 = figure(2);
d.stxm(res.pca.angle,struct('scale',0.02,'sampl',10,'unit','mm'));
hold on
d.pca(res.pca.angle);
set(imhandles(f2),'AlphaData',res.stxm.df>2e3);

f3 = figure(3);
d.stxm(res.pca.angle,struct('scale',0.02,'sampl',10,'unit','mm'));
hold on
d.pca(res.pca.angle,struct('quiver','on'));
set(imhandles(f3),'AlphaData',res.stxm.df>2e3);

%% click on a map to get the corresponding diffraction pattern
e.clicktool(d,[])
