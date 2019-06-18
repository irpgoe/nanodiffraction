function init = initialize()

if exist('nanodiffraction.m','file') ~= 2
    error('Please add the nanodiffraction toolbox to your MATLAB search path');
end

path_to_script = which('initialize.m');
path_to_script = fileparts(path_to_script);
path_to_dataset = fullfile(path_to_script,'data',filesep);

fprintf(1,'Data will be simulated and stored in the following directory:\n');
fprintf(1,'%s\n',path_to_dataset);

sny = 128;
snz = 128;

% initialize orientation field
x = linspace(-1.5,2.5,sny);           
y = linspace(-5,-1,snz);              
[X,Y] = meshgrid(x,y);
ZO = (Y.^2 + (X - 2).^2).^(-1/2) - ... % force field (electric field between
    (Y.^2 + (X - 1/2).^2).^(-1/2)/2;  % two charges)
[FX,FY] = gradient(ZO);                % gradient
O = atan2d(-FY,FX);                   % orientation (note that in a figure 
                                      % the y-direction goes 'down'
O(O<0) = O(O<0) + 180;                % project onto range [0 180)

% initialize darkfield
bgrLevel = 20;
variance = 0.1;
meanValue = 0;

x = linspace(-5,5,128);
y = linspace(-5,5,128);
[X,Y] = meshgrid(x,y);
Z = Y.*sin(X) - X.*cos(Y) + bgrLevel;   
DF = Z + sqrt(variance)*randn(size(Z)) + meanValue;

% calculate beamstop, air scattering and load detector mask
e = nanodiffraction();   
ny = e.Ny;
nz = e.Nz;

AIRBGR = e.simulate_airscattering(struct('I0',200,'l_post',0.05,'l_pre',20,'l_sd',e.detDistance*1000,'l_sampl',0.01));
BS = e.radial_mask('r1',-1,'r2',e.simulate_bs(0.2,150),'grid',e.qr);
MASK = zeros(nz,ny);
MASK([196:213 408:424],:) = 1;  % modular gaps
MASK(302:304,249:end) = 1;      % beamstop holder
MASK = BS | MASK;


init = struct();
init.sny = sny;
init.snz = sny;
init.mask = MASK;
init.bs = BS;
init.air = AIRBGR;
init.df = DF;
init.orientation = O;

save([path_to_script '/defaultValues.mat'],'init');


% create dataset
for d = 1:snz
    % data file
    mname = sprintf('%s%s%0.6d.h5',path_to_dataset,'unit_test_00001_data_',d);
    h5create(mname,'/entry/data/data',[ny nz sny],'Datatype','uint16','ChunkSize',[32 32 1],'Deflate',3);

    tic;
    for j = 1:sny
        I = AIRBGR + e.simulate_actomyosin(O(d,j)); % air background and structure
        I = I * DF(d,j);                            % sample thickness
        I = imnoise(I*1e-12, 'poisson')*1e12;       % noise
        I = I .* ~MASK;                             % mask (with beamstop)
        I = transpose(I);                           % eiger transposes data
        
        start = [1 1 j];
        count = [ny nz 1];
        h5write(mname,'/entry/data/data',uint16(I),start,count);
    end
    toc;
end

% master file
mname = sprintf('%s%s.h5',path_to_dataset,'unit_test_00001_master');
fid = H5F.create(mname);
g1 = H5G.create(fid,'entry','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
g2 = H5G.create(g1,'data','H5P_DEFAULT','H5P_DEFAULT','H5P_DEFAULT');
H5G.close(g1);
H5G.close(g2);
for d = 1:snz
    fname = sprintf('%s%s%0.6d.h5',path_to_dataset,'unit_test_00001_data_',d);
    H5L.create_external(fname,'/entry/data/data',fid,sprintf('%s%0.6d','/entry/data/data_',d),'H5P_DEFAULT','H5P_DEFAULT');
end
H5F.close(fid);
h5disp([path_to_dataset 'unit_test_00001_master.h5']);


              
d = display();
figure(1);
ax1 = subplot(1,3,1); 
    d.stxm(O,struct('sampl',20,'scale',1,'unit','um'));hold on;
    d.pca(O,struct('dir','v1'))
    hold on;
    contour(e.helper.clipLH(Z,-5,5),40,'LineWidth',2,'color','white'); 
    caxis([0 180]);    
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));
    title('simulated orientation field');
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    c = colorbar; ylabel(c,'orientation (degrees)');
    
subplot(1,3,2);
    d.stxm(ZO,struct('sampl',20,'scale',1,'unit','um')); caxis([0 0.2]); 
    hold on;
    contour(e.helper.clipLH(ZO,-5,5),40,'LineWidth',2,'color','white'); 
    d.add_quiver(O,struct('sampling',4,'scale',0.3,'dir','v1'));
    title('simulated force field')
    colormap(pmkmp);
    xlim([1 numel(x)]);
    ylim([1 numel(y)]);
    c = colorbar; ylabel(c,'intensity (a.u.)');
                
subplot(1,3,3);
    d.stxm(DF,struct('sampl',20,'scale',1,'unit','um'));
    title('simulated darkfield')
    c = colorbar; ylabel(c,'intensity (a.u.)');

    
colormap(ax1,cmap('C2'));


disp('To test the analysis, please execute the following lines:')
disp('f = files(''eigerConfig'',struct(''size'',[487 619]));');
disp('e = nanodiffraction();');
disp('d = display();');
disp('e.attach(f);');
disp('d.attach(e);');
disp('load(''nanodiffraction/initialize/defaultPilatusMask.mat'')')
disp('e.set_mask(mask);');
disp('e.set_selection(e.radial_mask(15, 40));');
disp('e.set_scan_info(128, 128);');
disp('e.set_roi_and_binning(struct(''roi'',''on'',''roiY'',e.around_y(150),''roiZ'',e.around_z(150)));');
disp('[df,angle,w] = split_struct(e.analyze_scan(''stxm+pca'',struct(''ySkip'',4,''zSkip'',4)),{''df'',''angle'',''w''});');
disp('figure(2);d.stxm(df);');
disp('figure(3);d.pca(angle);');

end

