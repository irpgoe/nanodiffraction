classdef display<handle
    % DISPLAY  Class for displaying standard STXM and XRF maps.
    %
    %   d = display()
    %
    % Note, that the display class should be linked to an instance of the
    % nanodiffraction class either through
    %   d.exp = (handle to nanodiffraction class);
    % or by using the link() function (see help link for usage information)
    %
    % The following arguments are accepted:
    %   axis:: [image]
    %       Default axis of images plotted
    %
    %   cmap:: [{gray,inverse}]
    %       Default colormap
    %
    %   figure:: [1]
    %       Figure number
    %
    %   colorbar:: [on]
    %       Adds a colorbar by default
    %
    %   caxis:: [auto]
    %       Automatic color scaling by default. Another useful option is to
    %       use the helper function autoc which scales to the 5 - 95
    %       percentile.
    %
    %   xlabel:: []
    %       Standard x-label.
    %
    %   ylabel:: []
    %       Standard y-label.
    %
    %  It contains the following functionality:
    %   cluster(): shows the result of a cluster analysis
    %   fluo(): displays fluorescence maps
    %   stxm(): displays darkfield images
    %   pca(): shows a darkfield with arrows superimposed to indicate the
    %   director of the anisotropic field
    %   composite(): shows a composite image
    %   diffraction(): shows a diffraction pattern
    %   add_circle(): adds circles at certain q_r to the diffraction
    %   pattern
    %   create_movie(): creates an avi movie
    %   saxs(): shows one or more 1d saxs curves
        
    % Copyright 2017 Institute for X-ray Physics (University of Göttingen)

    % Permission is hereby granted, free of charge, to any person obtaining 
    % a copy of this software and associated documentation files (the "Software"), 
    % to deal in the Software without restriction, including without limitation 
    % the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    % and/or sell copies of the Software, and to permit persons to whom the 
    % Software is furnished to do so, subject to the following conditions:

    % The above copyright notice and this permission notice shall be included 
    % in all copies or substantial portions of the Software.

    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
    % EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
    % MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
    % IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
    % DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    % TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
    % OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    properties
        p = struct; % holds parameters for plotting
        exp = nanodiffraction; % link to experiment for acces to axes and so forth
        lims = @(xlow,xhigh,ylow,yhigh) set(gca,'XLim',[xlow xhigh],'YLim',[ylow yhigh]); % wrapper for limits
        autoc = @(im) caxis([prctile(reshape(im,1,[]),5) prctile(reshape(im,1,[]),95)]); % autocontrast
        autonorm = @(im) (im - mean(reshape(im,1,[])))/(std(reshape(im,1,[]))); % autonormalization
    end
    
    methods
        % constructor
        function obj = display(varargin)
            
        	p = inputParser;
            % set defaults and add optional arguments
            addOptional(p,'axis','image');
            addOptional(p,'cmap',{'gray','inverse'});
            addOptional(p,'figure',1);
            addOptional(p,'q',1);
            addOptional(p,'title','');
            addOptional(p,'colorbar','on');
            addOptional(p,'caxis','auto');
            addOptional(p,'xlabel','');
            addOptional(p,'ylabel','');
            addOptional(p,'zlabel','');
            
            % parse arguments
            parse(p,varargin{:});
            % copy p.Results to obj.p
            names = fieldnames(p.Results);
            for ii = 1:numel(names)
                obj.p.(names{ii}) = p.Results.(names{ii});
            end
        end
        
        function add_line(obj,qr_values)
            % ADD_LINE  plots a vertical line at a specific qr value or
            % values.
            %   
            %   ADD_LINE(QR_VALUES) 
            %
            %   The following options are supported:
            %
            %     QR_VALUES:: []
            %       A single qr value or a 1d vector of values.
            %
            
            y_limits = ylim(gca);
            y = linspace(y_limits(1),y_limits(2),2);
            
            hold on;
            for ii = 1:numel(qr_values)
                x = linspace(qr_values(ii),qr_values(ii),2);
                line(x,y,'color','r');
            end
            hold off;
        end
        
        function cluster(obj,cluster_data,varargin)
            % CLUSTER  plots cluster map with relative contributions.
            %   
            %   CLUSTER(CLUSTER_DATA, MEANS) 
            %
            %   The following options are supported:
            %
            %     CLUSTER_DATA:: []
            %       Cluster map.
            %
            %     MEANS:: []
            %       Means of each cluster.
            %

            if nargin == 3
                means = varargin{1};
            end
            
            % number of clusters
            n_clusters = max(cluster_data(:));
            
            % show cluster as map
            subplot(1,2,1);
            imagesc(clusterres.clusters);
            colormap(pmkmp(n_clusters));
            axis image; 
            c=colorbar;
            ylabel(c,'cluster label');

            % show result as bar plot
            ax2 = subplot(1,2,2);
            colormap(ax2,gray(n_clusters));
            b = bar(means,'stacked');
            ylim([0 1]);
            xlim([0 n_clusters+1]);
            title('contribution of coefficients to clusters');xlabel('cluster label');ylabel('percentage'); 
            
            % calculate legend
            l_content = cell(1,n_clusters);
            for ii = 1:n_clusters
                l_content(ii) = sprintf('c_%d',ii-1);
            end
            l = legend(l_content);
            set(l,'Location','EastOutside');
            
        end
        
        
        % shows stxm scan results
        function stxm(obj, data, varargin)
            % STXM  plots parameter mappings with scale.
            %   
            %   STXM(DATA, OPTS) 
            %
            %   The following options are supported:
            %
            %     data:: []
            %       Parameter map.
            %
            %     opts:: []
            %       Struct that should contain a sampling value and a
            %       scale. It can in addition contain a custom unit. The
            %       default is microns.
            %
            %   Example:
            %   (displayobj).stxm(data,struct('sampl',<number>,...
            %                                 'scale',<number>,...
            %                                 'unit',<um|mm>));

            defaults = struct('sampl',1,'scale',1,'unit','um');
            fields = fieldnames(defaults);
            if nargin > 2
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            % short-hand notation
            sampl = opts.('sampl');
            scale = opts.('scale');
            unit = opts.('unit');
            
            % isotropic or anisotropic scaling
            if numel(scale) > 1
                try
                    scx = scale(1);
                    scy = scale(2);
                catch e
                    error(e);
                end
            else 
                scx = scale;
                scy = scale;
            end
            
            % show data
            imagesc(data);
                axis image;
                c = colorbar;
                % contrast and color
                obj.autoc(data);
                colormap(pmkmp(128));
                % label
                if strcmp(unit,'mm')
                    xlabel('y [mm]');
                    ylabel('z [mm]');
                else
                    % default: um
                    xlabel('y [\mum]');
                    ylabel('z [\mum]');
                end
                ylabel(c,'intensity [counts]');
                % axis tic labels
                xticks([1:sampl:size(data,2)]);
                xticklabels(([1:sampl:size(data,2)] - 1)*scx);
                yticks([1:sampl:size(data,1)]);
                yticklabels(([1:sampl:size(data,1)] - 1)*scy);
                
        end
        
        function add_quiver(obj,orientation,varargin)
            % ADD_QUIVER  Adds quiver lines onto an image, based on the
            % orientation.
            %
            %   QUIVER(ORIENTATION, VARARGIN)
            %   
            %   The following arguments are accepted:
            %       orientation:: []
            %           2d-array that contains angles in degrees confined
            %           to the interval [-90 90].
            %
            %       varargin:: []
            %           First argument: sampling ratio.
            %           Second argument: scale.
            %
            %   Example:
            %       (display_object).add_quiver(angles,struct('sampling',1,'scale',0.1));
            
            defaults = struct('sampling',1,'scale',0.1,'selection',[]);
            fields = fieldnames(defaults);
            if nargin > 2
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            % shorthand
            sampling = opts.sampling;
            scale = opts.scale;
            
            [x,y] = meshgrid(1:size(orientation,2),1:size(orientation,1));
            u = 1 .* sind(orientation); 
            v = 1 .* cosd(orientation);
            
            % sampling
            sampl = zeros(size(x));
            sampl(1:sampling:end,1:sampling:end) = 1;
            sampl = find(sampl);
            x = x(sampl);
            y = y(sampl);
            u = u(sampl);
            v = v(sampl);
            
            % selection
            if ~isempty(opts.selection)
                x = x(find(opts.selection));
                y = y(find(opts.selection));
                u = u(find(opts.selection));
                v = v(find(opts.selection));
            end
            
            hold on
            quiver(x,y,-v,u, scale, 'k','LineWidth',2,'ShowArrowHead','off');
            quiver(x,y,v,-u, scale, 'k','LineWidth',2,'ShowArrowHead','off');
            hold off
        end
        
        
        
        % create composite figure from scan
        function obj = composite(obj,comp,p)
            % COMPOSITE  shows a collection of diffraction patterns as a
            % composite with default styling.
            %   
            %   COMPOSITE(COMP, P) 
            %
            %   The following options are supported:
            %
            %     composite:: []
            %       A two-dimensional composite image.
            %
            %     p:: []
            %       Parameters that were used for composite calculation.
            %       p is obtained as a second argument from
            %       e.calculate_composite().
            
            figure(1); 
            imla(comp); 
            colormap(flipud(hot));
            c=colorbar;
            ylabel(c,'log_{10}(intensity/cps)');
            xlabel('frame along x')
            ylabel('frame along y')
            axis image;
            ax = gca;
            ax.XTick = [round(p.Ny/2):p.Ny:p.Ny*p.imsY];
            ax.YTick = [round(p.Nz/2):p.Nz:p.Nz*p.imsZ];
            ax.XTickLabel = 1:p.imsY;
            ax.YTickLabel = 1:p.imsZ;
            
        end
        
%         % create movie from data series
%         function obj = movie_from_figure(obj,filename,figurePosition,frameRate)
%             
%             % open figure
%             figure(f)
%             f=figure(1);
%             f.Position = [0 500 1000 350];
%             
%             % create video object and set frame rate
%             v = VideoWriter('newfile2.avi','Uncompressed AVI');
%             v.FrameRate = 15;
%             open(v);
%             
%             % get frame
%             frame = getframe(f);
%             
%             % write frame to file
%             writeVideo(v,frame)
%             
%             % close video object
%             close(v);
%         end
        
        
        function figureHandle = diffraction(obj,data,varargin)
            % DIFFRACTION  Displays diffraction patterns in a pre-defined
            % fashion.
            %
            %   DIFFRACTION(DATA, VARARGIN)
            %   
            %   The following arguments are accepted:
            %       data:: []
            %           Frame to be drawn.
            %
            %       process:: [on]
            %           Calls nanodiffraction.process when reading the
            %           data
            %
            %       qRange:: [-max(max(obj.exp.qr)) max(max(obj.exp.qr))]
            %           q-Range
            %
            %       qSteps:: [round((2*max(max(obj.exp.qr)))) / 10]
            %           Delta q between ticks
            %
            
            if nargin == 1
                warning('No arguments, choose either process, qRange, qSteps');
            end
            
            % qRange and Steps in rec. nanometers
            % parse input
            p = inputParser;
            addOptional(p,'process','on');
            addOptional(p,'qRange',[-max(max(obj.exp.qr)) max(max(obj.exp.qr))]);
            addOptional(p,'qSteps',round((2*max(max(obj.exp.qr)))) / 10 );
            parse(p,varargin{:});

            % q tick values
            qTicks = (floor(p.Results.qRange(1)/p.Results.qSteps)*p.Results.qSteps):p.Results.qSteps:(floor(p.Results.qRange(2)/p.Results.qSteps)*p.Results.qSteps);  
            
            % N as a function of Q
            N = obj.exp.helper.n_of_q(qTicks);
            
            % shift N according to primary beam
            Ny = round(N) + obj.exp.pby;
            Nz = round(N) + obj.exp.pbz;

            % calculate x/y q-axis
            qx = (obj.exp.helper.q_of_n((0:size(data,2)-1) - obj.exp.pby));
            qy = (obj.exp.helper.q_of_n((0:size(data,1)-1) - obj.exp.pbz));

            % show file
            if strcmp(p.Results.process,'on')
                obj.imlap(data);
            elseif strcmp(p.Results.process,'off')
                imla(data);
            else
                error('Argument has to be either "on" or "off"');
            end
            
            % tweak display
            colormap gray;
            colormap(flipud(colormap))
            c=colorbar;
            ylabel(c,'log_{10}(intensity/cps)');
            xlabel('q_r [nm^{-1}]')
            ylabel('q_r [nm^{-1}]')
            axis image;
            ax = gca;
            ax.XTick = Ny;
            ax.YTick = Nz;
            ax.XTickLabel = qTicks;
            ax.YTickLabel = -qTicks;
            
            % save data size that is currently holded
            obj.p.cols = size(data,2);
            obj.p.rows = size(data,1);
        end
        
        
        function add_circle(obj, qCircles, varargin)
            % ADD_CIRCLE  Superimposes circles on diffraction pattern.
            %
            %   ADD_CIRCLE(QCIRCLES, VARARGIN)
            %
            % The following arguments are accepted:
            %   qCircles:: []
            %       Qr value or vector of Qr values. Defines the radius of
            %       the circles to be plotted. Values are in reciprocal
            %       nanometers.
            %
            %   circleColor:: [black]
            %       Color of the circles to be plotted. Any default Matlab
            %       color can be chosen.
            
            
            % parse input
            p = inputParser;
            addOptional(p,'circleColor','black');
            parse(p,varargin{:});   
            
             % calculate q circles
            phi = 0:360;

            % assure that display is linked to experiment
            if isempty(obj.exp)
                warning('did you set display.exp to your experiment?')
            end
            
            for ii = 1 : numel(qCircles)

                % calculate circles
                R = obj.exp.helper.n_of_q(qCircles(ii)); % radius in pixels
                x = R*cosd(phi) + obj.exp.pby;
                y = R*sind(phi) + obj.exp.pbz;

                % filter all values that are outside the boundaries
                y(x>obj.p.cols) = []; x(x>obj.p.cols) = []; 
                x(y>obj.p.rows) = []; y(y>obj.p.rows) = []; 
                y(x<1) = []; x(x<1) = [];
                x(y<1) = []; y(y<1) = [];
                
                % show q circles
                hold on 
                line(x,y,'LineStyle','--','color',p.Results.circleColor,'LineWidth',1);
                hold off
            end
            figure(gcf);
        end
        
        function imlap(obj,frame)
            % IMLAP  Reading and processing of a frame.
            %
            %   IMLAP(FRAME)
            %
            % The following arguments are accepted:
            %   FRAME:: []
            %       Frame number within a scan.
            
            imla(obj.exp.process(frame));axis image;
        end
        
        function saxs(obj,saxsresult,varargin)
            % SAXS  Creates a loglog or semilog plot of the one-dimensional
            % structure factor I(qr) with default labeling of x- and
            % y-axis.
            %
            %   SAXS(SF, VARARGIN)
            %
            % The following arguments are accepted:
            %   SF:: []
            %       Structure factor, given as a struct that contains at
            %       least the fields "dat_1d" and "qr".
            %
            %   LIMITS:: []
            %       Limits [xlow xhigh ylow yhigh]. By default, this is set
            %       to "auto".
            %
            %   PIXEL:: [off]
            %       If the structure factor should be plotted as a function
            %       of detector pixel, PIXEL should be set to "on". In this
            %       case, the SF should contain a field "pixel".
            %
            %   MODE:: [semilogy]
            %       Either semilogy or loglog are accepted.
            
            
            % parse input
            p = inputParser;
            addOptional(p,'limits',[]);
            addOptional(p,'pixel','off');
            addOptional(p,'mode','semilogy');
            parse(p,varargin{:});  
            
            switch p.Results.mode
                case 'semilogy'
                    plot_fn = @semilogy;
                case 'loglog'
                    plot_fn = @loglog;
            end
            legendstr = {};
            for ii = 1:numel(saxsresult)
                if iscell(saxsresult)
                tmp = saxsresult{ii};
                elseif ~iscell(saxsresult) && numel(saxsresult) == 1
                    tmp = saxsresult;
                else
                    error('something is wrong');
                end
                if strcmp(p.Results.pixel,'off')
                    plot_fn(tmp.qr,tmp.dat_1d,'o-','MarkerSize',2);hold all;grid on;
                    xlabel('q_r [nm^{-1}]');
                elseif strcmp(p.Results.pixel,'on')
                    plot_fn(tmp.pixel,tmp.dat_1d,'o-','MarkerSize',2);hold all;grid on;
                    xlabel('pixel');
                else
                end
                ylabel('log_{10}(intensity/counts)');
                if isempty(p.Results.limits)
                    axis auto;
                else
                    axis(p.Results.limits);
                end
                legendstr{ii} = num2str(ii);
            end
            hold off;
            legend(legendstr);
        end
        
        function cmap = single_hue_colormap(obj,hue,sat,cmap) 
            % SINGLE_HUE_COLORMAP  Creates a colormap based on a single hue
            % value.
            %
            %   SINGLE_HUE_COLORMAP(HUE, SAT, COLORMAP)
            %
            % The following arguments are accepted:
            %   HUE:: []
            %       Hue given in degrees.
            %
            %   SAT:: []
            %       Saturation is given a value between 0 and 1.
            %
            %   COLORMAP:: []
            %       desc.
            
%             h = ones(n,1).*hue/360;
%             s = ones(n,1).*sat;
%             v = linspace(0,1,n)';
%             s = linspace(1,0,n)';
%             cmap = [h s v];
            
            test = rgb2hsv(rgb2gray(cmap));
            test(:,1) = hue/360; test(:,2) = sat;
            
            % shorten
            test(1:32,:) = test(1:2:end,:);
            test(33:end,:) = test(1:2:end,:);
            test(33:end,2) = linspace(sat,0,32);
            test(33:end,3) = 1;

            cmap = hsv2rgb(test);
        end
        
        function newmap = azimuthal_colormap(obj,varargin)
            % AZIMUTHAL_COLORMAP  Returns a colormap suitable for
            % orientation visualization. Note, that by default, hsv is used
            % as colormap, however, the more uniform colormap by Peter
            % Kovesi is recommended.
            %   
            %   AZIMUTHAL_COLORMAP(COLORMAP_NAME)
            %   
            %   The input argument COLORMAP_NAME can take the following 
            %   values:
            %       'pmkmp'     :: based on the matlab package 
            %                      'Perceptually improved colormaps' by 
            %                      Matteo Niccoli and the adaptation by 
            %                      Mike Bostock
            %                      (https://bl.ocks.org/mbostock/310c99e53880faec2434)
            %                      (https://github.com/d3/d3-scale)
            %       'isoAz'     :: uses pmkmp
            %       'custom'    :: manually defined colormap
            %       'custom2'   :: manually defined colormap, with black
            %                      and white interchanged
            %       'aike'      :: Based on http://cgm.technion.ac.il/people/Viki/figure1&5_cluster/computeColor.m
            %                      and http://people.seas.harvard.edu/~dqsun/#home
            %       'peterkovesi' :: Based on "Peter Kovesi. Good Colour
            %                        Maps: How to Design Them. 
            %                        arXiv:1509.03700 [cs.GR] 2015"
            
            
            if nargin == 2
                base_map = varargin{1};
            else
                base_map = 'hsv';
            end
            switch base_map
                case 'hsv'
                    newmap = hsv;
                case 'pmkmp'
                    H = rgb2hsv(pmkmp(180));
                    icut = find((H(1,1) - H(:,1)) > 0.5,1);
                    H = H(1:icut,:);
                    Hold = hsv2rgb(H);
                    H(:,1) = mod(H(:,1) + 0.5,1);
                    H(:,2) = flipud(H(:,2));
                    H(:,3) = flipud(H(:,3));
                    H = hsv2rgb(H);
                    newmap = [Hold;H];
                case 'isoAz'
                    newmap = pmkmp(180,'IsoAZ').^0.7;
                case 'custom'
                    newmap = customangles;
                case 'custom2'
                    newmap = customangles2;
                case 'aike'
                    newmap = aikescolor;
                case 'peterkovesi'
                    newmap = cmap('C2');
                otherwise
                    error('Wrong input, please refer to help display.azimuthal_colormap for more information');
            end
        end
        
        function pca(obj,angle,varargin)
            % PCA  Creates a visual representation of the map of
            % orientations. By default, angles are visualized using color,
            % however, lines can be superimposed. Color can also be
            % weighted by saturation.
            %
            %   PCA(ANGLE, VARARGIN)
            %
            % The following arguments are accepted:
            %   ANGLE:: []
            %       Map of angles from PCA analysis.
            %
            %   VARARGIN:: []
            %       Struct, that can contain several plotting options.
            %       Option 1: 'dir' ('v1' | 'v2')
            %           By default, the input angle represents the
            %           orientation of the diffraction
            %           streak/modulation/peak. However, by default, the
            %           angle is rotated by 90 degrees to reflect the
            %           orientation of the fibre axis.
            %           Using the 'dir' option, one can choose 'v2' or 'v1'
            %           as options that indicate which eigenvector
            %           direction should be used.
            %       Option 2: 'quiver' ('on' | 'off')
            %           'quiver': 'on' will superimpose lines on each data
            %           point to indicate fibre or streak orientation.
            %       Option 3: 'alpha' (alphaMap)
            %           'alpha' can be given a map of alpha values that can
            %           be used for mapping color saturation to another
            %           parameter map. If an alpha map is given, it will be
            %           automatically used.
            
            switch nargin 
                case 1
                    error('Please provide a nxm matrix with entries within the range [-pi,pi)');
                case 2
                    opts = struct();
                    % everything ok. Quiver option not chosen.
                case 3
                    opts = varargin{1};
            end
            
            % default direction is along the fiber axis
            if isfield(opts,'dir')
                if strcmp(opts.dir,'v2')
                    angle = angle - 90;
                    normDir = true;
                else
                    normDir = false;
                end
            else 
                % default
                angle = angle - 90;
                normDir = true;
            end
            
            if isfield(opts,'quiver')
                if strcmp(opts.quiver,'on')
                    % add quiver to plot
                    doQuiver = true;
                else 
                    doQuiver = false;
                end
            else
                % default
                doQuiver = false; 
            end
            
            if isfield(opts,'alpha')
                doAlpha = true;
            else 
                doAlpha = false;
            end

            
            % show orientation
            if isfield(opts,'axHandle')
                % update figure
                set(get(opts.axHandle,'Children') ,'CData',angle);
                h = gca;
            else
                % create new figure
                h = imagesc(angle);
            end
            
            c = colorbar; 
            axis image;
            colormap(obj.azimuthal_colormap('peterkovesi'))
            
            % alpha mask
            if doAlpha
                set(h,'AlphaData',flipud(rot90(opts.alpha)));
            end
            
            % real or reciprocal space direction?
            if normDir
                caxis([-90 90]);
                ylabel(c,'fibre orientation [degrees]');
            else
                caxis([0 180]);
                ylabel(c,'reflex orientation [degrees]');
            end

            % add quiver lines
            if doQuiver
                obj.add_quiver(angle,struct('sampling',3,'scale',0.3));
            end
        end
        
        function image_overlay(obj,im1,im2,transp)
            % IMAGE_OVERLAY  Overlays to images. The second image is placed
            % on top of the first image, and the transparency channel of
            % the second image can be tuned for overlay.
            %
            %   IMAGE_OVERLAY(IMAGE1, IMAGE2, TRANSP)
            %   
            %   The following arguments are accepted:
            %       image1:: []
            %           Base image.
            %
            %       image2:: []
            %           Overlay image.
            %
            %       transp:: []
            %           Transparency of image2.
            %            
            %   to create rgb images use e.g.:
            %       rgbimage1 = ind2rgb(round(mat2gray(data)*63 + 1),testColormap);
            %   or for a single color overlay
            %       rgbimage2 = ind2rgb(ones(256)*16,testColormap);
            %   a colormap can be defined by e.g.
            %       testColormap = ds.single_hue_colormap(265,0.4,gray);

            % autoscale
            imagesc(im1);

            % overlay second image
            hold on
            h = imagesc(im2);axis image;
            hold off

            % tune transparency
            set(h,'AlphaData',transp);
        end
        
        
        function show_location(obj,scanpoint,varargin)
            % SHOW_LOCATION  highlights the location of a given scanpoint
            % in a 2d parameter map.
            %   
            %   SHOW_LOCATION(SCANPOINT, OPTS) 
            %
            %   The following options are supported:
            %
            %     scanpoint:: []
            %       A single index, corresponding to the location in the
            %       scan. Note, that the index starts with 1 in the
            %       top-left corner of the scan.
            %
            %     opts:: []
            %       Work in progress. Currently, only a struct with the
            %       direction as a parameter can be given.
            %       
            %       Usage: struct('dir',<'col'|'row'>))
            %
            %   Example: 
            %   imagesc(rand(60,21));
            %   ds.show_location(150,struct('dir','row'));
            %   Note that ds.show_location(150) would also work.
            
            
            defaults = struct('dir','row');
            fields = fieldnames(defaults);
            if nargin > 2
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            % get current figure
            axhandle = get(gca,'Children');
            if numel(axhandle) > 1
                axhandle = axhandle(end);
            end
            [sizeY,sizeX] = size(axhandle.CData);
            
            switch opts.dir
                case 'col'
                [y,x] = ind2sub([sizeY, sizeX],scanpoint);
                case 'row'
                [x,y] = ind2sub([sizeX, sizeY],scanpoint);
                otherwise
                error('Please specify the direction using <col> or <row>');
            end
            
            figure(gcf);
            hold on;
            line([x-0.5 x+0.5],[y+0.5 y+0.5],'color','r','LineWidth',2);
            line([x-0.5 x+0.5],[y-0.5 y-0.5],'color','r','LineWidth',2);
            line([x-0.5 x-0.5],[y-0.5 y+0.5],'color','r','LineWidth',2);
            line([x+0.5 x+0.5],[y-0.5 y+0.5],'color','r','LineWidth',2);
            hold off;
            
        end
    end
end