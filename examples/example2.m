clear all; close all; clc;

% add toolbox
addpath(genpath('toolbox'));

% make path definitions, add custom directories
fp = files( 'beamline','id13',...
            'prepath', 'homegroups/Labdata-Archive/AG_Salditt/Messzeiten_Rohdaten/2016/extern/ESRF_ID13_SC4304/id13/inhouse/DATA/AUTO-TRANSFER/eiger1',...
            'newfile','detdistcalib',...
            'detector','eiger',...
            'eigerNr',210);

% experiments
e = nanodiffraction('energy',13E3,...
                    'detDistance',5.12,...
                    'pixelsize',172E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1500.499,...
                    'pbz',1372.509); 

% connect to data
link(fp,e)
                
% set detector mask
e.set_mask('detector',fp.read(1)>4e9);

% scan parameters
e.set_scan_info('SNy',31,...
                'SNz',11);

% analyze scan
detdistcalc = e.analyze_scan('ZCrop',[6 8],'method','average+stxm');

% show average pattern
ds = display();
ds.diffraction(detdistcalc.avg);


% write file to tif, then read tif with fabio
% fp.export_tif(detdistcalc.avg,'test_for_pyfai')

%%  healing
e = nanodiffraction('energy',13E3,...
                    'detDistance',5.12,...
                    'pixelsize',172E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',900,...
                    'pbz',900); 

% set pb
e.pby = 25;
e.pbz = 25;

e.set_roi_and_binning('binning','on','biny',36,'binz',36);
C = e.phi;
% indexlist
il = transpose(1:numel(e.phi));
il = reshape(il,size(e.phi,1),size(e.phi,2));
% imagesc(il)

bad = zeros(size(C));
bad(10:20,10:20) = 1;
[symm,good,bad] = e.make_mask_symmetric(bad);
figure(1);
subplot(1,3,1);imagesc(bad);axis image;
subplot(1,3,2);imagesc(good);axis image;
subplot(1,3,3);imagesc(symm);axis image;

% new indexlist
symm = logical(symm);
good = logical(good);
bad = logical(bad);
il_heal = il;
il_bad = il(bad);
il_good = rot90(il(good),2);
il_heal(bad) = il_good;

% show result
figure(2);
subplot(2,2,1);imagesc(il);axis image;
subplot(2,2,2);imagesc(reshape(il_bad,11,11));axis image;
subplot(2,2,3);imagesc(il_heal);axis image;
subplot(2,2,4);imagesc(reshape(il_good,11,11));axis image;

figure(3);
subplot(1,2,1);imagesc(C);axis image;
subplot(1,2,2);imagesc(C(il_heal));axis image;

%% a practical example
e = nanodiffraction('energy',13E3,...
                    'detDistance',5.12,...
                    'pixelsize',172E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1500.499,...
                    'pbz',1372.509); 
e.set_mask('detector',fp.read(1)>4e9);
% e.set_roi_and_binning('detectorRoi','on','detRoiY',[1100 1900],'detRoiZ',[972 1772]);
% e.set_mask('detector',e.process_mask(e.masks));
% data = e.process(detdistcalc.avg);
data = detdistcalc.avg;

% indexlist
il = transpose(1:e.Ny*e.Nz);
il = reshape(il,e.Nz,e.Ny);

bad = e.masks;
[symm,good,bad,good_zero,bad_zero] = e.make_mask_symmetric(bad);
f1 = figure(1);
subplot(2,3,1);imagesc(bad);axis image;title('bad')
subplot(2,3,2);imagesc(good);axis image;title('good')
subplot(2,3,3);imagesc(symm);axis image;title('symmetric')
subplot(2,3,4);imagesc(bad_zero);axis image;title('bad zero')
subplot(2,3,5);imagesc(good_zero);axis image;title('good zero')

% new indexlist
symm = logical(symm);
good = logical(good);
bad = logical(bad);
bad_zero = logical(bad_zero);
good_zero = logical(good_zero);
il_heal = il;
il_heal(bad_zero) = flipud(il(good_zero));

% show result
f2 = figure(2);
subplot(1,2,1);imagesc(il);axis image;title('initial indexlist');
subplot(1,2,2);imagesc(il_heal);axis image;title('repaired indexlist');

f3 = figure(3);
subplot(1,3,1);imagesc(data);axis image;caxis([0 1]);title('broken frame');
subplot(1,3,2);imagesc(data(il_heal));axis image;caxis([0 1]);title('repaired frame');
subplot(1,3,3);imagesc(e.heal(data,e.masks));axis image;caxis([0 1]);title('repaired frame - shorthand');


%% now heal real data
% experiments
e = nanodiffraction('energy',13E3,...
                    'detDistance',5.12,...
                    'pixelsize',172E-6,...
                    'Ny',2070,...
                    'Nz',2167,...
                    'pby',1500.499,...
                    'pbz',1372.509); 

% connect to data
link(fp,e)
                
% set detector mask
e.set_mask('detector',fp.read(1)>4e9);

% scan parameters
e.set_scan_info('SNy',31,...
                'SNz',11);

% % analyze scan
% detdistcalc = e.analyze_scan('ZCrop',[6 8],'method','average+stxm');
% 
% % show average pattern
% ds = display();
% ds.diffraction(detdistcalc.avg);

% analyze scan
detdistcalc = e.analyze_scan('ZCrop',[6 8],'method','average+stxm+heal');

% show average pattern
ds = display();
ds.diffraction(detdistcalc.avg);