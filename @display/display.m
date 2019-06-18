classdef display<handle
    % DISPLAY  Class for displaying standard STXM and XRF maps.
    %
    %   ``d = display()``
    %
    % Note, that the display class should be linked to an instance of the
    % nanodiffraction class by using the ``display.attach()`` function (see
    % ``help attach`` for usage information).
    %
    % The following arguments are accepted:
    %   axis:: ['image']
    %       Default axis of images plotted
    %
    %       Options: Default Matlab options.
    %
    %   cmap:: [{gray,inverse}]
    %       Default colormap
    %
    %   colorbar:: ['on']
    %       Adds a colorbar by default
    %
    %       Options: 'on'|'off'
    %
    %   caxis:: ['auto']
    %       Automatic color scaling by default. Another useful option is to
    %       use the helper function autoc which scales to the 5 - 95
    %       percentile.
    %
    %   xlabel:: ['']
    %       Standard x-label.
    %
    %   ylabel:: ['']
    %       Standard y-label.
    %
    %   zlabel:: ['']
    %       Standard z-label.
    %
    %  It contains the following functionality:
    %   add_circle():
    %       adds circles at certain q_r to the diffraction pattern
    %   add_line():
    %       adds a vertical line to a SAXS curve
    %   add_quiver():
    %       superimposes a PCA result with quiver lines
    %       indicating the orientation of the scattering of the structure
    %       orientation.    
    %   autoc():
    %       Auto-contrast based on the 5% percentile.
    %   autonorm(): 
    %       Normalizes a distribution (z-transform).
    %   azimuthal_colormap(): 
    %       Several colormaps to plot phase angles are
    %       available.
    %   composite(): 
    %       displays a composite image.
    %   diffraction(): 
    %       shows a diffraction pattern    
    %   imlap(): 
    %       processes an image and then displays it on a log. scale.
    %   image_overlay(): 
    %       Overlays an image with a transparent second image.
    %   saxs(): 
    %       shows one or more 1d saxs curves    
    %   scalebar(): 
    %       Adds a scalebar (in pixel units) to an image.
    %   show_location(): 
    %       Highlights a scan point in a map based on its
    %       linear index.
    %   stxm(): 
    %       displays darkfield images    
    %   pca(): 
    %       shows a darkfield with arrows superimposed to indicate the
    %       director of the anisotropic field
    %    
    % Copyright 2017 Institute for X-ray Physics (University of Göttingen)
    %
    % Permission is hereby granted, free of charge, to any person obtaining 
    % a copy of this software and associated documentation files (the "Software"), 
    % to deal in the Software without restriction, including without limitation 
    % the rights to use, copy, modify, merge, publish, distribute, sublicense, 
    % and/or sell copies of the Software, and to permit persons to whom the 
    % Software is furnished to do so, subject to the following conditions:
    %
    % The above copyright notice and this permission notice shall be included 
    % in all copies or substantial portions of the Software.
    %
    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
    % EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF 
    % MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
    % IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
    % DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
    % TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE 
    % OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
    properties (Access = private)
        exp = [];               % handle to experiment class
    end
    
    properties (Access = public)
        p = struct;             % holds parameters for plotting
        
        axis = 'image';
        cmap = {'gray','inverse'};
        colorbar = 'on';
        caxis = [];
        xlabel = '';
        ylabel = '';
        zlabel = '';
        scale = 'log';
        
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
            addOptional(p,'colorbar','on');
            addOptional(p,'caxis',[]);
            addOptional(p,'xlabel','');
            addOptional(p,'ylabel','');
            addOptional(p,'zlabel','');
            addOptional(p,'scale','log');
            
            % parse arguments
            parse(p,varargin{:});
            
            % copy structure fields to properties of class
            names = fieldnames(p.Results);
            for ii = 1:numel(names)
                obj.(names{ii}) = p.Results.(names{ii});
            end
        end
        
        function new = copy(obj)
            % COPY  Makes a copy of a handle object. This keeps track of 
            % all the given parameters so that multiple instances of a 
            % display module can be used.
            %
            %   ``new_instance = copy(display_object)``
            %
            % The following arguments are supported:
            %   
            %   display_object: [] (required)
            %       This function requires a display object that should be
            %       copied.
            %
            % Example:
            %   See the following example for help::
            %            
            %      d_generic = display();
            %      d_saxsdata = copy(d_generic);
            %      d_waxsdata = copy(d_generic);
            %
            % Output arguments:
            %   new_instance: 
            %       A copy of the initial display_object.
            %
            
            
            % Instantiate new object of the same class.
            new = feval(class(obj));
            
            % Copy all non-hidden properties.
            prop = properties(obj);
            for i = 1:length(prop)
                new.(prop{i}) = obj.(prop{i});
            end
        end
        
        
        function attach(obj,nanodiffraction_handle)
            % ATTACH  Attaches a nanodiffraction-module to this class.
            %
            %   ``attach(nanodiffraction_handle)``
            %
            % The following arguments are supported:
            %   files_handle: [] (required)
            %       Handle to a files-module.
            %
            % Example: 
            %   See the following example for help::
            %
            %       e = nanodiffraction(...);
            %       d = display(...);
            %       d.attach(e);
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
            
            obj.exp = nanodiffraction_handle;
        end       
        
        
        function add_line(obj,qr_values)
            % ADD_LINE  plots a vertical line at a specific location given
            % in units of the reciprocal wavevector (inverse nanometers).
            % This function is used to plot vertical lines in SAXS plots.
            % The vertical lines run through the entire y-range.
            %   
            %   ``add_line(qr_values)``
            %
            % The following arguments are supported:
            %
            %     qr_values: [] (required)
            %       A one-dimensional vector of values in units of inverse
            %       nanometers. Each entry specifies the location of a
            %       vertical line that is superimposed onto the current
            %       plot.
            %
            % Example:
            %   See the following example for help::
            %            
            %       d = display();
            %       d.saxs(some_saxs_data);
            %       d.add_line([0.2 0.6 1.2]);
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            y_limits = ylim(gca);
            y = linspace(y_limits(1),y_limits(2),2);
            
            hold on;
            for ii = 1:numel(qr_values)
                x = linspace(qr_values(ii),qr_values(ii),2);
                line(x,y,'color','r');
            end
            hold off;
        end
        
        
        function add_uncertainty(obj,x_pos,y_pos,sigma,varargin)
            % ADD_LINE  plots a vertical line at a specific location given
            % in units of the reciprocal wavevector (inverse nanometers).
            % This function is used to plot vertical lines in SAXS plots.
            % The vertical lines run through the entire y-range.
            %   
            %   ``add_line(qr_values)``
            %
            % The following arguments are supported:
            %
            %     x_pos: [] (required)
            %       A one-dimensional vector of values in units of inverse
            %       nanometers. Each entry specifies the location of a
            %       vertical line that is superimposed onto the current
            %       plot.
            %
            % Example:
            %   See the following example for help::
            %            
            %       d = display();
            %       d.saxs(some_saxs_data);
            %       d.add_line([0.2 0.6 1.2]);
            %
            % Output arguments:
            %   This function does not return any arguments.

            defaults = struct('dx',0.01,'lineWidth',2,'xScale','lin');
            fields = fieldnames(defaults);
            if nargin > 4
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
            dx = opts.dx;
            
            hold on;
            for ii = 1:numel(x_pos)
                line([x_pos(ii) x_pos(ii)]  ,[-sigma(ii) +sigma(ii)] + y_pos(ii),'color','k','LineWidth',opts.lineWidth); % vertical line
                if strcmp(opts.xScale,'lin')
                    line([-dx +dx] + x_pos(ii),[+sigma(ii) +sigma(ii)] + y_pos(ii),'color','k','LineWidth',opts.lineWidth); % top line
                    line([-dx +dx] + x_pos(ii),[-sigma(ii) -sigma(ii)] + y_pos(ii),'color','k','LineWidth',opts.lineWidth); % bottom line
                elseif strcmp(opts.xScale,'log')
                    ldx = 10^(log10(x_pos(ii))-log10(x_pos(ii)-dx));    % lower dx
                    udx = 10^(log10(x_pos(ii)+dx) - log10(x_pos(ii)));  % upper dx
                    line([x_pos(ii)-ldx x_pos(ii)+udx],[+sigma(ii) +sigma(ii)] + y_pos(ii),'color','k','LineWidth',opts.lineWidth); % top line
                    line(10.^log10([x_pos(ii)-dx x_pos(ii)+dx]),[-sigma(ii) -sigma(ii)] + y_pos(ii),'color','k','LineWidth',opts.lineWidth); % bottom line
                else
                    warning('no horizontal bars');
                end
            end
            hold off;
            
            
        end
        
        
        function f = figure_size(obj,varargin)
            
            pos = [0 0 10 8];
            if nargin == 2
                pos = varargin{1};
            end
            
            f = gcf;
            f.PaperUnits = 'centimeters';
            f.PaperPosition = pos;
            f.Units = 'centimeters';
            f.Position = pos;
        end
        
        function ax = remove_axis(obj,varargin)
            % REMOVE_AXIS  removes the axis and ticks from the active
            % figure.
            %   
            %   ``remove_axis()``
            %
            % The function does not require any arguments.
            %
            % Example:
            %   See the following example for help::
            %               
            %       d = display();
            %       d.stxm(some_map);
            %       d.remove_axis();
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            opt = 'xy';
            if nargin > 1
                opt = varargin{1};
            end

            
            switch opt
                case 'xy'
                    ax = gca;
                    ax.XTick = [];
                    ax.XTickLabel = [];
                    ax.YTick = [];
                    ax.YTickLabel = [];
                    ax.XLabel = [];
                    ax.YLabel = [];
                case 'all'
                    ax = gca;
                    ax.XTick = [];
                    ax.XTickLabel = [];
                    ax.YTick = [];
                    ax.YTickLabel = [];
                    ax.XLabel = [];
                    ax.YLabel = [];
                    colorbar(ax,'delete');
            end
        end
        
        
        function add_scalebar(obj, sl_pixel, varargin)
            % ADD_SCALEBAR  adds a scale bar to the current figure.
            %   
            %   ``add_scalebar(sl_pixel, opts)`` 
            %
            % The following options are accepted:
            %     sl_pixel: [] (required)
            %       Scale bar length in pixel units.
            %
            %     opts: [see default values below] (optional)
            %       Structure that can contain the following fields. Note,
            %       that all fields are optional. Default values that are
            %       otherwise used are as usual given in angle brackets:
            %
            %       sb_height: [0.2]
            %           height of the scalebar in figure units.
            %
            %       margin_right: [0.05]
            %           right margin of the scalebar in figure units.
            %
            %       margin_bottom: [0.05]
            %           bottom margin of the scalebar in figure units.
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       d.stxm(some_map);
            %       d.add_scalebar(100); % scale bar length = 100 scan points * 2 µm steps
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            defaults = struct('sb_height',0.2,'margin_right',0.05,'margin_bottom',0.05,'color','k');
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
            
            axhandle = get(gca,'Children');
            if numel(axhandle) > 1
                axhandle = axhandle(end);
            end
            [sizeY,sizeX] = size(axhandle.CData);
            
            % scalebar length and height in pixels
            sb_len = round(sl_pixel);
            sb_height = round(opts.sb_height*sb_len);
            
            % scalebar margin right and bottom
            sb_off_x = round(sizeX*opts.margin_right);
            sb_off_y = round(sizeY*opts.margin_bottom);
            
            % superimpose an image consisting of -infs to ensure the sb is 
            % black and plot scale bar as transparency layer
            hold on; 
            plot([sizeX-sb_off_x-sb_len; sizeX-sb_off_x], ...
                 [sizeY-sb_off_y; sizeY-sb_off_y], ...
                 '-', 'color', opts.color, 'LineWidth', sb_height);
            hold off;
        end

        
        function stxm(obj, data, varargin)
            % STXM  plots parameter maps with with pre-defined settings.
            %   
            %   ``stxm(data, opts)``
            %
            % The following options are supported:
            %
            %     data: [] (required)
            %       Parameter map.
            %
            %     opts: [see default values below] (optional)
            %       Structure that can contain the following fields. Note,
            %       that all fields are optional. Default values that are
            %       otherwise used are as usual given in angle brackets:
            %
            %           sampl: [100] 
            %               A tick is added to every N-th pixel/scan point.
            %
            %           scale: [1]
            %               Length of one pixel/scan point in units 'unit'.
            %
            %           unit: [mm]
            %               Either 'um' or 'mm'.
            %
            %           alpha: []
            %               An alpha/transparency map, values of the map should
            %               range within (0,1).
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       d.stxm(data,struct('sampl',100,...
            %                          'scale',10,...
            %                          'unit','um'));
            %
            % Output arguments:
            %   This function does not return any arguments.

            defaults = struct('sampl',10,'scale',1,'unit','pixel','alpha',[]);
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
            
            if ~isempty(opts.alpha)
                doAlpha = true;
            else 
                doAlpha = false;
            end
            
            % show data
            h=imagesc(data);
                axis image;
                
                % colorbar
                c = colorbar;
                c.LineWidth = 1.5;
                c.TickDirection = 'out';
                ylabel(c,'intensity (counts)','Interpreter','Latex');
                
                % contrast and color
                obj.autoc(data);
                colormap(pmkmp(128));
                set(gca,'TickDir','out');
                set(gca,'LineWidth',1.5,'TickLength',[0.02 0.02]);
            
                % label
                if strcmp(unit,'mm')
                    xlabel('y (mm)','Interpreter','Latex');
                    ylabel('z (mm)','Interpreter','Latex');
                elseif strcmp(unit,'um')
                    xlabel('y ($\mathrm{\mu m}$)','Interpreter','Latex');
                    ylabel('z ($\mathrm{\mu m}$)','Interpreter','Latex');
                else
                    xlabel('y (pixel)','Interpreter','Latex');
                    ylabel('z (pixel)','Interpreter','Latex');
                end
                
                % axis tic labels
                xticks([1:sampl:size(data,2)]);
                xticklabels(([1:sampl:size(data,2)] - 1)*scx);
                yticks([1:sampl:size(data,1)]);
                yticklabels(([1:sampl:size(data,1)] - 1)*scy);
                
            % alpha mask
            if doAlpha
                set(h,'AlphaData',opts.alpha);
            end
                
        end
        
        
        function [limits] = round_limits(obj,varargin)
            % ROUND_LIMITS  rounds the caxis limits of the current figure to 
            % a given decimal digit level:
            %
            % | e.g.: round(1342,-3) => 1000
            % | e.g.: round(0.042,2) => 0.04
            %
            %   ``round_limits(level)``
            %
            % The following arguments are supported:
            %   level: [1] (optional)
            %       The rounding level. The same level as used in the
            %       Matlab-round function.
            %
            % Example:
            %   See the following example for help::
            %
            %      d = display();
            %      d.stxm(some_map);
            %      d.round_limits(-3);
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            if nargin > 1
                level = varargin{1};
            else
                level = 1;
            end
            
            % get limits
            tmp = gca;
            limits = tmp.CLim;
            
            % adjust limits
            lim_tmp = round(limits,level);
            
            i = 1;
            while lim_tmp(1) == lim_tmp(2) && i < 10
                lim_tmp = round(limits,level+i);
                i  = i + 1;
            end
            
            limits = lim_tmp;
            caxis([limits(1) limits(2)]);
        end
        
        
        function add_title(obj, varargin)
            % ADD_TITLE  Adds a title to the image that includes the caxis
            % range and the length of a scalebar.
            %
            %   ``add_title(opts)``
            %   
            % The following arguments are supported:
            %   opts: [see default values below] (optional)
            %       Structure that can contain the following fields.
            %       Note, that all fields are optional. Default values
            %       that are otherwise used are as used given in angle
            %       brackets:
            %
            %       limits: [auto]
            %           Per default choses the limits of the current axis.
            %           Otherwise, a two-element vector with a range can be
            %           used, e.g. [1 2].
            %
            %       sb_len: []
            %           Length of the scale bar in µm.
            %
            %       subtitle: [false]
            %           The information can be plotted as a regular title
            %           (false, the default)
            %           to a figure or subplot or can be inserted as an
            %           xlabel (true).
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       d.add_title(struct('sb_len',5,'subtitle',true));
            %
            % Output arguments:
            %       This function does not return any arguments.
            
            defaults = struct('limits','auto','sb_len',[],'subtitle',false);
            fields = fieldnames(defaults);
            if nargin > 1
                opts = varargin{1};
                for f = 1:numel(fields)
                    if ~isfield(opts,fields{f})
                        opts.(fields{f}) = defaults.(fields{f});
                    end
                end
            else
                opts = defaults;
            end
            
            title_fun = @title;
            if opts.subtitle
                title_fun = @xlabel;
            end
            
            if strcmp(opts.limits,'auto')
                ca = gca;
                opts.limits = ca.CLim;
            end
            
            if strcmp(opts.limits,'none')
                limit_str = '';
            else
                limit_str = sprintf('%g -> %g', opts.limits(1), opts.limits(2));
            end
                
            if ~isempty(opts.sb_len)
                sb_str = sprintf('%g µm', opts.sb_len);
            else
                sb_str = '';
            end
            title_fun([limit_str ' : ' sb_str])
            
        end

        
        function add_quiver(obj,orientation,varargin)
            % ADD_QUIVER  Adds quiver lines onto an image, based on the
            % orientation.
            %
            %   ``add_quiver(orientation, opts)``
            %   
            % The following arguments are supported:
            %   orientation: [] (required)
            %       2d-array that contains angles in degrees confined
            %       to the interval [-90 90].
            %
            %   opts: [see default values below] (optional)
            %       Structure that can contain the following fields.
            %       Note, that all fields are optional. Default values
            %       that are otherwise used are as used given in angle
            %       brackets:
            %
            %       sampling: [1]
            %           plot quiver lines every n-th pixel/scan point.
            %
            %       scale: [0.1]
            %           length of quiver lines.
            %
            %       selection: []
            %           display quiver lines only on selected
            %           scanpoints. selection is a 2d logical array.
            %
            %       dir: ['v2']
            %           Can either be set to 'v2' or 'v1' where 'v2'
            %           denotes the fiber direction and 'v1' the scattering
            %           direction, respectively.
            %
            % Example:
            %   See the following example for help::
            %            
            %       d = display();
            %       d.add_quiver(orientation,struct('sampling',1,'scale',0.1));
            %
            % Output arguments:
            %       This function does not return any arguments.
            
            
            defaults = struct('sampling',1,'scale',0.1,'selection',[],'dir','v2');
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
        
            % default direction is along the fiber axis
            if strcmp(opts.dir,'v2')
                orientation = orientation - 90;
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
                opts.selection = opts.selection(sampl);
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
        
        
        
        function obj = composite(obj,comp,p)
            % COMPOSITE  shows a collection of diffraction patterns as a
            % composite with default styling.
            %   
            %   ``composite(comp, p)``
            %
            % The following options are supported:
            %   composite: [] (required)
            %       A two-dimensional composite image.
            %
            %   p: [see default values below] (required)
            %       Parameters that were used for composite calculation.
            %       p is obtained as a second argument from
            %       e.calculate_composite(). It can also be defined
            %       manually. p is a structure that can contain the
            %       following fields:
            %       
            %       Ny: [] (required)
            %           Number of pixels along the horizontal dimension of a 
            %           single frame.
            %
            %       Nz: [] (required)
            %           Number of pixels along the vertical dimension of a 
            %           single frame.
            %
            %       imsY: [] (required)
            %           Number of images displayed along the horizontal
            %           dimension of the composite image.
            %
            %       imsZ: [] (required)
            %           Number of images displayed along the vertical
            %           dimension of the composite image.
            %
            %       yCrop: [1 obj.scan.SNy] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       zCrop: [1 obj.scan.SNz] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       ySkip: [1] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       zSkip: [1] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       correction: ['off'] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       emptySub: ['off'] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       empty: [] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       heal: ['off'] (optional)
            %           Switch that toggled healing on or off for the 
            %           calculation of the composite image. This value is 
            %           optional, but useful for future reference.
            %
            %       healmask: [] (optional)
            %           Mask that has been used for healing of the diffraction
            %           patterns. This value is optional, but useful for future
            %           reference.
            %
            %       binning: [] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       detRoiY: [] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            %       detRoiZ: [] (optional)
            %           Description. This value is optional, but useful for 
            %           future reference.
            %
            % Example:
            %   See the following example for help::
            %
            %       e = nanodiffraction(<some settings here>);
            %       e.set_scan_info(<some settings here>);
            %       [comp,p] = e.calculate_composite(...
            %                       struct('ySkip',4,'zSkip',4));
            %       d = display();
            %       d.composite(comp,p);
            %
            % Output arguments:
            %   This function does not return any arguments.
            
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
            % DIFFRACTION  Displays a diffraction patterns with pre-defined
            % settings. Simplifies the correct addition of labels to the x-
            % and y-axis of the figure. By default, diffraction patterns
            % are shown in log-scale and display in grayscale (a linear
            % scale that does not highlight any particular q-range).
            %
            %   ``diffraction(data, opts)``
            %   
            % The following arguments are supported:
            %   data: [] (required)
            %       A two-dimensional diffraction pattern. By default,
            %       this diffraction pattern is automatically processed
            %       to be in accordance with the currently active
            %       settings for a detector ROI and a binning factor.
            %       Processing can however be switched off using the
            %       second optional argument 'opts'.
            %
            %   opts: [see default values below] (optional)
            %       Structure that can contain the following fields.
            %       Note, that all fields are optional. Default values
            %       that are otherwise used are as usual given in angle
            %       brackets:
            %
            %       process: ['on'] 
            %           Either 'on', 'off' or 'force'. If activated
            %           ('on'), then data of the size of a raw detector
            %           frame is expected and the data will be
            %           automatically binned and cropped based on the
            %           detector roi and binning settings
            %
            %       qRange: [-max(max(e.qr)) max(max(e.qr))]
            %           Lower and upper limit for the detector q-range to 
            %           be plotted on the axis of the figure. Typically,
            %           the full q-range of the detector is used. The
            %           q-range of the currently linked nanodiffraction
            %           object (shorthand: e) is used for this purpose.
            %
            %       qSteps: [round((2*max(max(obj.exp.qr)))) / 10]
            %           Distance between ticks in units of the radial 
            %           wavevector transfer. The q-range of the currently 
            %           linked nanodiffraction object (shorthand: e) is 
            %           used for this purpose.
            %
            %       alpha: []
            %           A two-dimensional alpha map.
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       e = nanodiffraction();
            %       frame = e.read(1); 
            %       d.diffraction(frame,struct('process','off','qSteps',5);
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            if nargin == 1
                warning('Please provide a data set.');
            end
            
            % qRange and Steps in rec. nanometers
            defaults = struct('process','on',...
                              'qRange',[-max(max(obj.exp.qr)) max(max(obj.exp.qr))],...
                              'qSteps',round((2*max(max(obj.exp.qr)))) / 10,...
                              'alpha',[]);
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

            % q tick values
            qTicks = (floor(opts.qRange(1)/opts.qSteps)*opts.qSteps):opts.qSteps:(floor(opts.qRange(2)/opts.qSteps)*opts.qSteps);  
            
            % N as a function of Q
            N = obj.exp.n_of_q(qTicks);
            
            % shift N according to primary beam
            Ny = round(N) + obj.exp.pby;
            Nz = round(N) + obj.exp.pbz;

            % calculate x/y q-axis
            qx = (obj.exp.q_of_n((0:size(data,2)-1) - obj.exp.pby));
            qy = (obj.exp.q_of_n((0:size(data,1)-1) - obj.exp.pbz));

            % check whether data needs to be processed
            [dz,dy] = size(data);
            [rz,ry] = deal(obj.exp.roiZ(2) - obj.exp.roiZ(1) + 1,...
                obj.exp.roiY(2) - obj.exp.roiY(1) + 1);
            if isequal([dz,dy],[rz,ry]) && ~strcmp(opts.process,'force')
                warning('Data appears to be already processed, processing will be turned off. To enforce processing, use "force" as an argument.')
                opts.process = 'off';
            end
            
            switch obj.scale
                case 'log'
                    im = @(data) imla(data);
                case 'lin'
                    im = @(data) imagesc(data);
                otherwise
                    im = @(data) imla(data);
            end 
            
            % show file
            if strcmp(opts.process,'on') || strcmp(opts.process,'force')
                data = obj.exp.process(data);
            elseif strcmp(opts.process,'off')
            else
                error('Argument has to be either "on", "off" or "force"');
            end
            im(data);
%             drawnow;
            h = gca; 
            
            % tweak display
            axis image;
            colormap gray;
            colormap(flipud(colormap))
            c=colorbar;
            ylabel(c,'$\log_{10}$(intensity/cps)','Interpreter','Latex');
            xlabel('$\mathrm{q_r}$ ($\mathrm{nm^{-1}}$)','Interpreter','Latex')
            ylabel('$\mathrm{q_r}$ ($\mathrm{nm^{-1}}$)','Interpreter','Latex')
            axis image;
            if isempty(obj.caxis)
                caxis('auto');
            else
                caxis(obj.caxis);
            end
            ax = gca;
            ax.XTick = Ny;
            ax.YTick = Nz;
            ax.XTickLabel = qTicks;
            ax.YTickLabel = -qTicks;
            
            % alpha mask
            if ~isempty(opts.alpha)
                set(h.Children,'AlphaData',opts.alpha);
            end
            
            % save data size that is currently holded
%             obj.p.cols = size(data,2);
%             obj.p.rows = size(data,1);
        end
        
        
        function add_circle(obj, qCircles, varargin)
            % ADD_CIRCLE  Superimposes circles on a diffraction pattern.
            % Multiple circle radii can be given to the function. The radii
            % are expected in units of rec. nanometers.
            %
            %   ``add_circle(qCircles, opts)``
            %
            % The following arguments are accepted:
            %   qCircles: [] (required)
            %       One-dimensional array of circle radii in units of the
            %       reciprocal wavevector (inverse nanometers).
            %
            %   opts: [see default values below] 
            %       Structure that can contain the following fields. Note
            %       that all fields are optional. Default values that are
            %       otherwise used are given as usual in angle brackets
            %       below:
            %
            %       circleColor:: [black]
            %           Color of the circles to be plotted. Any default Matlab
            %           color can be chosen.
            %
            %       wedge:: [0 360] 
            %           One-dimensional array of circle radii in units of the
            %           reciprocal wavevector (inverse nanometers).
            %       
            %
            % Example:
            %   See the following example for help::
            %   
            %       d = display();
            %       e = nanodiffraction(); 
            %       frame = e.read(1);
            %       d.diffraction(frame);
            %       d.add_circle([0.2 0.6 1.2],'white');
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            if nargin == 1
                warning('Please provide circle positions.');
            end
            
            % qRange and Steps in rec. nanometers
            defaults = struct('circleColor','black',...
                              'wedge',[0 360],...
                              'nrWedges',2,...
                              'circleStyle','--');
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
            
             % calculate q circles
            phi = opts.wedge(1):1:opts.wedge(2);
            
            switch opts.nrWedges
                case 1
                case 2
                    opts.nrWedges = 1;
                    opts.wedge = opts.wedge+180;
                    obj.add_circle(qCircles,opts);
                case 4
                    opts.nrWedges = 1;
                    wedge1 = opts.wedge+180;
                    wedge2 = opts.wedge+90;
                    wedge3 = opts.wedge-90;
                    opts.wedge = wedge1;
                    obj.add_circle(qCircles,opts);
                    opts.wedge = wedge2;
                    obj.add_circle(qCircles,opts);
                    opts.wedge = wedge3;
                    obj.add_circle(qCircles,opts);
                otherwise
                    error('nrWedges should be either 1,2 or 4');
            end

            % assure that display is linked to experiment
            if isempty(obj.exp)
                warning('Did you set display.exp to your experiment?')
            end
            
            % size of image
            axhandle = get(gca,'Children');
            if numel(axhandle) > 1
                axhandle = axhandle(end);
            end
            [rows,cols] = size(axhandle.CData);
            
            for ii = 1 : numel(qCircles)

                % calculate circles
                R = obj.exp.n_of_q(qCircles(ii)); % radius in pixels
                x = R*cosd(phi) + obj.exp.pby;
                y = R*sind(phi) + obj.exp.pbz;

                % filter all values that are outside the boundaries
                y(x>cols) = []; x(x>cols) = []; 
                x(y>rows) = []; y(y>rows) = []; 
                y(x<1) = []; x(x<1) = [];
                x(y<1) = []; y(y<1) = [];
                
                % show q circles
                hold on 
                line(x,y,'LineStyle',opts.circleStyle,'color',opts.circleColor,'LineWidth',1);
                hold off
            end
            figure(gcf);
        end
        
        function imlap(obj,frame)
            % IMLAP  Displays a detector image in logscale. The detector
            % image is hearby automatically processed, based on the currect
            % setting of the detector ROI and binning factors. 
            %
            %   ``imlap(frame)``
            %
            % The following arguments are supported:
            %   frame: [] (required)
            %       Detector image. The image should have identical
            %       dimensions as the detector. 
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       e = nanodiffraction();
            %       frame = e.read(1);
            %       d.imlap(frame);
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            imla(obj.exp.process(frame));axis image;
        end
        
        function saxs(obj,saxsresult,varargin)
            % SAXS  Creates a loglog or semilog plot of the one-dimensional
            % structure factor I(qr) with default labeling of x- and
            % y-axis.
            %
            %   ``saxs(sf, opts)``
            %
            % The following arguments are supported:
            %   sf: [] (required)
            %       Structure factor, given as a structure that contains at
            %       least the fields "dat_1d" and "qr". Both fields should
            %       contain a one-dimensional array that should have the
            %       same lengths. A cell array of structures is also
            %       supported if multiple one-dimensional saxs curves
            %       should be plotted.
            %
            %   opts: [see default values below] (optional)
            %       Structure that can contain the following fields. Note,
            %       that all fields are optional. Default values that are
            %       otherwise used are as usual given in angle brackets:
            %
            %       limits: []
            %           Limits [xlow xhigh ylow yhigh]. By default, this is set
            %           to "auto".
            %
            %       pixel: ['off']
            %           If the structure factor should be plotted as a function
            %           of detector pixel, 'pixel' should be set to "on". In 
            %           this case, the SF should contain a field "pixel".
            %
            %       mode: ['semilogy']
            %           Either 'semilogy' or 'loglog' are accepted.
            %
            % Example:
            %   See the following example for help::
            %            
            %       d = display();
            %       d.saxs({my_1d_averaged_data,...
            %          struct('dat_1d',some_1d_data,'qr',some_qr_axis)},...
            %          struct('mode','loglog'));
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            % parse input
            defaults = struct('limits',[],...
                              'pixel','off',...
                              'mode','semilogy',...
                              'markerSize',2,...
                              'colormap',@pmkmp);
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
            
            
            switch opts.mode
                case 'lin'
                    plot_fn = @plot;
                case 'semilogy'
                    plot_fn = @semilogy;
                case 'loglog'
                    plot_fn = @loglog;
            end
            legendstr = {};
            
            % pmkmp requires at least 2 colors
            if numel(saxsresult) == 1
                cm = opts.colormap(2);
                cm = cm(1,:);
            else
                cm = opts.colormap(numel(saxsresult));
            end
            
            for ii = 1:numel(saxsresult)
                if iscell(saxsresult)
                tmp = saxsresult{ii};
                elseif ~iscell(saxsresult) && numel(saxsresult) == 1
                    tmp = saxsresult;
                else
                    error('something is wrong');
                end
                
                % check for correct fieldnames
                if isfield(tmp,'qr')
                    x = tmp.qr;
                elseif isfield(tmp,'x')
                    x = tmp.x;
                else
                    error('x-axis not found in structure');
                end
                    
                if isfield(tmp,'dat_1d')
                    y = tmp.dat_1d;
                elseif isfield(tmp,'y')
                    y = tmp.y;
                else
                    error('y-axis not found in structure');
                end
                
                % pixel option
                if strcmp(opts.pixel,'off') 
                    plot_fn(x,y,'o','color',cm(ii,:),'MarkerSize',opts.markerSize,'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor',cm(ii,:));hold all;grid on;
                    xlabel('q_r (nm^{-1})');
                elseif strcmp(opts.pixel,'on')
                    plot_fn(tmp.pixel,y,'o','color',cm(ii,:),'MarkerSize',opts.markerSize,'MarkerFaceColor',cm(ii,:),'MarkerEdgeColor',cm(ii,:));hold all;grid on;
                    xlabel('pixel');
                else
                end
                ylabel('intensity (counts)');
                if isempty(opts.limits)
                    axis auto;
                else
                    axis(opts.limits);
                end
                legendstr{ii} = num2str(ii);
            end
            hold off;
            legend(legendstr);
        end
        
        
        function newmap = azimuthal_colormap(obj,varargin)
            % AZIMUTHAL_COLORMAP  Returns a colormap suitable for
            % orientation visualization. Note, that by default, hsv is used
            % as colormap, however, the more uniform colormap by Peter
            % Kovesi is recommended.
            %   
            %   ``azimuthal_colormap(colormap_name)``
            %   
            % The following arguments are supported:
            %   colormap_name: ['hsv'] (required)
            %       The following colormap names are currently implemented:
            %
            %       :'hsv': default HSV colormap
            %       :'pmkmp': based on the matlab package 
            %                 'Perceptually improved colormaps' by 
            %                 Matteo Niccoli and the adaptation by 
            %                 Mike Bostock
            %                 (https://bl.ocks.org/mbostock/310c99e53880faec2434)
            %                 (https://github.com/d3/d3-scale)
            %       :'isoAz': uses pmkmp
            %       :'peterkovesi': Based on "Peter Kovesi. Good Colour
            %                       Maps: How to Design Them. 
            %                       arXiv:1509.03700 [cs.GR] 2015"
            %
            % Example:
            %   See the following example for help::
            %
            %       d = display();
            %       d.stxm(some_stxm_map);
            %       d.remove_axis();
            %       cmap = d.azimuthal_colormap('pmkmp');
            %       colormap(gca,cmap);
            %
            % Output arguments:
            %   new_map: 
            %       Colormap based on the choice made using the
            %       colormap name.
            %
            
            
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
            %   ``pca(angle, opts)``
            %
            % The following arguments are supported:
            %   angle: [] (required)
            %       Map of angles that is e.g. output from a PCA analysis.
            %
            %   opts: [see default values below] (optional)
            %       Structure that can contain the following fields. Note,
            %       that all fields are optional. Default values that are
            %       otherwise used are as usual givin in angle brackets.
            %
            %       dir: ['v2']
            %           By default, the input angle represents the
            %           orientation of the diffraction
            %           streak/modulation/peak. However, by default, the
            %           angle is rotated by 90 degrees to reflect the
            %           orientation of the fibre axis.
            %           Using the 'dir' option, one can choose 'v2' or 'v1'
            %           as options that indicate which eigenvector
            %           direction should be used.
            %
            %       quiver: ['off']
            %           If set to 'on' quiver lines will be superimposed on
            %           each data point to indicate fibre or streak 
            %           orientation.
            %           Note, that default quiver arguments are used. For a
            %           more custom quiver display, please use
            %           display.add_quiver.
            %
            %       alpha: []
            %           alpha can be given a map of alpha values that can
            %           be used for mapping color saturation to another
            %           parameter map. If an alpha map is given, it will be
            %           automatically used.
            %
            % Example:
            %   See the following example for help::
            %
            %      d.display();
            %       d.pca(some_orientation_map,struct('quiver','on');
            %
            % Output arguments:
            %   This function does not return any arguments.
            %
            
            if nargin == 1
                warning('Please provide a data set.');
            end
            
            defaults = struct('sampl',10,...
                              'scale',1,...
                              'unit','pixel',...
                              'dir','v2',...
                              'quiver','off',...
                              'alpha',[],...
                              'axHandle',[]);
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
            
            % default direction is along the fiber axis
            if strcmp(opts.dir,'v2')
                angle = angle - 90;
                normDir = true;
            else
                normDir = false;
            end

            % add quiver to plot
            if strcmp(opts.quiver,'on')
                doQuiver = true;
            else 
                doQuiver = false;
            end
            
            % add overlay
            if ~isempty(opts.alpha)
                doAlpha = true;
            else 
                doAlpha = false;
            end

            % show orientation
            if ~isempty(opts.axHandle)
                % update figure
                set(get(opts.axHandle,'Children') ,'CData',angle);
                h = gca;
            else
                % create new figure
                h = imagesc(angle);
                axis image;
            end
            
            % colorbar            
            c = colorbar; 
            c.LineWidth = 1.5;
            c.TickDirection = 'out';
             
            % real or reciprocal space direction?
            if normDir
                caxis([-90 90]);
                ylabel(c,'fibre orientation (degrees)','Interpreter','Latex');
            else
                caxis([0 180]);
                ylabel(c,'reflex orientation (degrees)','Interpreter','Latex');
            end
            
            % contrast and color
            colormap(obj.azimuthal_colormap('peterkovesi'))
            set(gca,'TickDir','out');
            set(gca,'LineWidth',1.5,'TickLength',[0.02 0.02]);
            
            % label
            if strcmp(unit,'mm')
                xlabel('y (mm)','Interpreter','Latex');
                ylabel('z (mm)','Interpreter','Latex');
            elseif strcmp(unit,'um')
                xlabel('y ($\mathrm{\mu m}$)','Interpreter','Latex');
                ylabel('z ($\mathrm{\mu m}$)','Interpreter','Latex');
            else
                xlabel('y (pixel)','Interpreter','Latex');
                ylabel('z (pixel)','Interpreter','Latex');
            end

            % axis tic labels
            xticks([1:sampl:size(angle,2)]);
            xticklabels(([1:sampl:size(angle,2)] - 1)*scx);
            yticks([1:sampl:size(angle,1)]);
            yticklabels(([1:sampl:size(angle,1)] - 1)*scy);
            
            % alpha mask
            if doAlpha
                set(h,'AlphaData',opts.alpha);
            end
            
            % add quiver lines
            if doQuiver
                obj.add_quiver(angle,struct('sampling',1,'scale',0.3));
            end
        end
        
        function image_overlay(obj,im1,im2,transp)
            % IMAGE_OVERLAY  Overlays two images. The second image is 
            % placed on top of the first image, and the transparency 
            % channel of the second image can be tuned for overlay.
            % Note, that the images should be given in rgb color. To create 
            % rgb images use e.g.:
            %
            %       rgbimage1 = ind2rgb(round(mat2gray(data)*63 + 1),gray);
            %
            % or for a single color overlay
            %
            %       rgbimage2 = ind2rgb(ones(256)*16,gray);
            %
            % a colormap can be defined by e.g.::
            %
            %       d = display();
            %       testColormap = d.single_hue_colormap(265,0.4,gray);
            %
            %   ``image_overlay(image1, image2, transparency)``
            %   
            % The following arguments are supported:
            %       image1: [] (required)
            %           Base image.
            %
            %       image2: [] (required)
            %           Overlay image.
            %
            %       transparency: []
            %           Transparency of image2. Values should range within
            %           the interval (0,1).
            %            
            % Example:
            %   Example is missing.
            %
            % Output arguments:
            %   This function does not return any arguments.


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
            % in a two-dimensional parameter map.
            %   
            %   ``show_location(scanpoint, opts)`` 
            %
            % The following options are supported:
            %     scanpoint: [] (required)
            %       A single index, corresponding to the location in the
            %       scan. Note, that the index starts with 1 in the
            %       top-left corner of the scan.
            %
            %     opts: [] (optional)
            %       Structure that can contain the following field. Note 
            %       that currently, only a single field can be set, which
            %       is likely to be made more flexible in the future.
            %
            %       dir: ['row']
            %           Specifies the direction along the single index should
            %           be incremented. Either 'col' or 'row' are allowed.
            %
            % Example: 
            %   See the following example for help::
            %   
            %       d = display();
            %       imagesc(rand(60,21));
            %       d.show_location(150,struct('dir','row'));
            %
            % Output arguments:
            %   This function does not return any arguments.
            
            
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