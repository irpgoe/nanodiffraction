% Currently, read_from_dada is deactivated in production.

clear all; close all; clc;

addpath(genpath('toolbox'))

% make path definitions, add custom directories
fp = files();
        
frame1 = fp.read_from_dada(); % default image
frame2 = fp.read_from_dada('experiment','ls2522',...
                            'detector','eiger1',...
                            'filenumber',11,...
                            'framenumber',2,...
                            'instrument','ESRF',...
                            'roi','1090,1111,200,200',...
                            'binning','1,1');

figure(1); imla(frame1); axis image; c=colorbar; ylabel(c,'log_{10}(intensity/counts)');
figure(2); imla(frame2); axis image; c=colorbar; ylabel(c,'log_{10}(intensity/counts)');