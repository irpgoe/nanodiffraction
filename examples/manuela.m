%% %%
% % % 1) Start the analysis
% % %                
clear all; close all; clc;

% add toolbox
addpath(genpath('toolbox'));

% are you on santa or on dancer?
prepath = which_computer('santa');

% detector modules
f = files( 'beamline','id13',...
            'prepath','/home/Labdata-Archive/AG_Koester/Messzeiten_Rohdaten/2016/ESRF_ID13_SC4406_Dec2016/sc4406/DATA/AUTO-TRANSFER/eiger1',...
            'newfile','device21_colloids',...
            'detector','eiger',...
            'eigerNr',1451);

e = nanodiffraction('energy',13.9E3,...
                    'detDistance',0.95,...
                    'pixelsize',75E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1354,...
                    'pbz',1261);

e.set_scan_info('SNy',41,'SNz',101);

% set detector mask 
ff = f.read(1); 
ff(4157307) = 2e9; ff(948,1291) = 2e9;
ff(1,:) = 2e9; ff(end,:) = 2e9; ff(:,1) = 2e9; ff(:,end) = 2e9;
ff([514 552 1065 1103 1616 1654],:) = 2e9;
ff(:,[1030 1041]) = 2e9;
e.set_mask('detector',ff>1e9);

% display routines
d = display();     % (s)how (s)axs
link(f,e,e,d);

% background
f.set_eiger_scan('device21_buffer',1448);
selection_buffer = e.analyze_scan('method','average','yCrop',[18 20],'zCrop',[29 37]);

% stxm analysis
f.set_eiger_scan('device21_colloids',1451);
test = e.analyze_scan('method','stxm');
bloeder_bs = e.analyze_scan('method','average','yCrop',[18 27],'zCrop',[90 100]);
selection = e.analyze_scan('method','average','yCrop',[18 20],'zCrop',[29 37]);
d.diffraction(selection.avg)

% azimuthal average
az = e.azimuthal_mask('phi1',120,'phi2', 160);
selection_avg = b1d(selection.avg-selection_buffer.avg,e.masks | az,[],e.qr);
d.saxs(selection_avg)

% visualize
d.stxm(test.stxm.df,struct('sampl',10,'scale',1,'unit','um'))
e.clicktool(d,[])

% composite
e.set_roi_and_binning('detectorRoi','on','detRoiY',[1250 1450],'detRoiZ',[1110 1400]);
e.set_mask('detector',e.process_mask(e.masks));
e.set_mask('correction',e.process_mask(e.corr));
[comp, p] = e.calculate_composite(struct('ySkip',[1 4],'zSkip',[1 4]));
d.composite(comp,p);caxis([0 3]);

% save progress
save_progress('manuela','results/manuela','analysis of device21_colloids');