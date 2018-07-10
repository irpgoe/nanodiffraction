function [colored,hue,sat,val,ph] = hsv_orientation(varargin)
% Two basic arguments are accepted: phase and magnitude. Note, that the
% phase is expected to be scaled between 0 and 360 degrees
% The magnitude will be scaled between 0 and 1 automatically
% basic scaling functions
%
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

    scale = @(amp) (amp - min(min(amp)))/(max(max(amp)) - min(min(amp)));
    scaleH = @(amp,high) scale(amp).*high;
    scaleLH = @(amp,low,high) scaleH(amp,high-low) + low;
    

    if nargin == 0
        warning('No arguments given, using example data');
        
        % test data
        ph = scaleLH(peaks,0,360);
        ph  = ph(1:48,1:48)/360;
        amp = diff(peaks);
        amp = amp(1:48,1:48);
        
        hue = scale(ph);
        sat = ones(size(hue));
        val = scale(amp);
        
    elseif nargin == 1
        warning('Interpreting data as phase');
        
        ph = varargin{1}/360;
        amp = ones(size(ph));
        
        hue = scaleLH(ph,min(min(ph)),max(max(ph)));
        sat = ones(size(hue));
        val = amp;
        
    elseif nargin == 2 
        ph = varargin{1}/360;
        amp = varargin{2};
        
        hue = scaleLH(ph,min(min(ph)),max(max(ph)));
        sat = ones(size(hue));
        val = scale(amp);
        
    elseif nargin > 2
        
        % amplitude and phase
        ph = varargin{1};
        % check second argument carefully
        noAmpGiven = 0;
        testvar = varargin{2};
        type = whos('testvar');
        if strcmp(type.class,'char')
            noAmpGiven = 1;
            amp = ones(size(ph));
        else
            amp = varargin{2};
        end
        
        pin = inputParser;
        % set defaults and add optional arguments
        % clipping
        addOptional(pin,'clipAmp',[]);
        addOptional(pin,'clipAmpLow',0);
        addOptional(pin,'clipAmpHigh',0);
        addOptional(pin,'clipPhase',[]);
        addOptional(pin,'clipPhaseLow',0);
        addOptional(pin,'clipPhaseHigh',0);
        % colormap
        addOptional(pin,'colormap','');
        % limits
        addOptional(pin,'val',[0 1]);              
        addOptional(pin,'hue',[0 360]);
        % parse arguments
        if noAmpGiven
            parse(pin,varargin{2:end}); 
        else
            parse(pin,varargin{3:end}); 
        end
        
        
        % clip amplitude
        if pin.Results.clipAmpLow 
            amp(amp < pin.Results.clipAmpLow) = pin.Results.clipAmpLow;
        elseif pin.Results.clipAmpHigh
            amp(amp > pin.Results.clipAmpHigh) = pin.Results.clipAmpHigh;
        elseif ~isempty(pin.Results.clipAmp) 
            amp(amp < pin.Results.clipAmp(1)) = pin.Results.clipAmp(1);
            amp(amp > pin.Results.clipAmp(2)) = pin.Results.clipAmp(2);
        end

        % clip phase
        if pin.Results.clipPhaseLow 
            ph(ph < pin.Results.clipPhaseLow) = pin.Results.clipPhaseLow;
        elseif pin.Results.clipPhaseHigh
            ph(ph > pin.Results.clipPhaseHigh) = pin.Results.clipPhaseHigh;
        elseif ~isempty(pin.Results.clipPhase) 
            ph(ph < pin.Results.clipPhase(1)) = pin.Results.clipPhase(1);
            ph(ph > pin.Results.clipPhase(2)) = pin.Results.clipPhase(2);
        end
        
        % scaling hue
        if ~strcmp(pin.Results.colormap,'')
            cmap = colormap(pin.Results.colormap);
            tmp = round(scaleLH(ph,1,64)); % value is now between 1 and 64
            
            hue = zeros(size(ph));
            [rows,cols] = size(hue);
            for rr = 1:rows
                for cc = 1:cols
                    hue(rr,cc) = rgb_to_hsv(cmap(tmp(rr,cc),:));
                    % hue is now between 0  and 360
                end
            end
            hue = hue./360;
        else
            hue = scaleLH(ph,pin.Results.hue(1)/360,pin.Results.hue(2)/360);
            % hue is now between 0  and 360
        end
        
        % scaling 
        if noAmpGiven == 0
            sat = scaleLH(amp,pin.Results.val(1),pin.Results.val(2));
        elseif noAmpGiven
            sat = amp;
        end
        
        % scaling 
        if noAmpGiven == 0
            val = scaleLH(amp,pin.Results.val(1),pin.Results.val(2));
        elseif noAmpGiven
            val = amp;
        end
        
    end

    colored = zeros(size(hue,1),size(hue,2),3);
    colored(:,:,1) = hue;
    colored(:,:,2) = sat;
    colored(:,:,3) = val;

    figure,set(gcf,'Position',[0 0 1400 500]);
    ax1 = subplot(1,5,1);imagesc(val);c=colorbar;axis image;
    colormap(ax1,'gray');c.Location = 'SouthOutside';title('amplitude');
    ax2 = subplot(1,5,2);imagesc(ph);c=colorbar;axis image;
    colormap(ax2,'parulaangles');c.Location = 'SouthOutside';title('phase angle');
    ax3 = subplot(1,5,3);imagesc(hsv2rgb(colored));c=colorbar;axis image;
    colormap(ax3,'parulaangles');c.Location = 'SouthOutside';title('amplitude and phase');
    ax4 = subplot(1,5,4);histogram(ph);axis square;
    colormap(ax4,'parulaangles');c.Location = 'SouthOutside';title('distribution of angles');
    ax5 = subplot(1,5,5);histogram(val);axis square;
    colormap(ax5,'parulaangles');c.Location = 'SouthOutside';title('distribution of amplitude');
end