function [hIm,hCB]=imageseis(Data,Time,X)
% [hIm,hCB]=ImageSeis(Data,Time,X) displays seismic data Data, using a
% default blue-white-red interactive colorbar. Data should be a matrix of
% time by trace number.
%
% Optional Parameters:
% - Time should be a vector of time samples used for the y-axis of the plot.
% - X should be a vector of space samples used for the y-axis of the plot.
%
% Optional Output:
% - hIm is a handle to the image
% - hCB is a handle to the interactive colorbar
%
% Using the Interactive Colorbar
% - Dragging the black bars at the top and bottom of the colorbar
%   symmetricaly changes the clipping level.
% - When right clicking on the colorbar one has the following options:
%   - Gray:              Changes the colormap to grayscale
%   - Inverted Gray:     Changes the colormap to inverted grayscale
%   - Seismic:           Changes the colormap to seismic blue-white-red
%                        This is the default colormap
%   - Colormap Editor:   Opens the matlab colormap editor
%   - Set Limits:        Sets the current colorbar clipping level to the
%                        colorbar limits. This allows for more fine tuned
%                        clipping around the center of the colorbar.
%   - Reset Limits:      Resets colorbar limits to their original values.
%   - Symetric Clipping: Toggles if clipping levels should be symmetric
%                        around the center of the colorbar.
%   - Connect Colorbar:  Connect interactive colorbar with other
%                        interactive colorbars so that they all clip the
%                        same.
%
% To connect multiple interactive colorbars connect their handles using:
% - set(hCB,'ConnectedColorBarObjects',hCB_other);
%   - hCB is the handle to the interactive colorbar which you want to
%     connect to hCB_other.
%   - hCB_other is the handle to the interactive colorbar to which you want
%     to connect hCB.
% To deactivate a connection:
% - set(hCB,'UpdateConnectedColorBars',false);
%   - hCB is the handle to the colorbar which you do not want to be updated
%     anymore and which you do not want to update other colorbars with.

% Max Holicki - designed for AESB2140
% October 2014
% m.e.holicki@tudelft.nl

% This software makes use of:
% InteractiveColorbar.m (Modified)
% hey dude
% http://www.mathworks.nl/matlabcentral/fileexchange/45363-interactivecolorbar (13.10.2014)
%
% seismic.m (Modified, incorporated into InteractiveColorbar.m)
% Community Surface Dynamics Modeling System
% https://csdms.colorado.edu/svn/sedflux/sedflux-mfiles/sedflux/seismic.m (13.10.2014)
%% Check Matrix Size
if ~ismatrix(Data)||isscalar(Data)||~isnumeric(Data);
    fprintf(2,'ImageSeis: Data is not a 2D numeric matrix\n');
    return;
end;
if exist('Time','var')&&(~isvector(Time)||numel(Time)~=size(Data,1)||~isnumeric(Time));
    fprintf(2,'ImageSeis: Time is not a number of time samples numeric vector\n');
    return;
end;
if exist('X','var')&&(~isvector(X)||numel(X)~=size(Data,2)||~isnumeric(X));
    fprintf(2,'ImageSeis: X is not a number of space samples numeric vector\n');
    return;
end;
%% Plot
if exist('Time','var');
    if ~exist('X','var');
        hIm=imagesc(1:size(Data,2),Time,Data,max(max(abs(Data)))*[-1,1]);
        xlabel('Trace #');
    else
        hIm=imagesc(X,Time,Data,max(max(abs(Data)))*[-1,1]);
        xlabel('Horizontal Distance');
    end;
    ylabel('Time');
else
    if ~exist('X','var');
        hIm=imagesc(1:size(Data,2),1:size(Data,1),Data,max(max(abs(Data)))*[-1,1]);
        xlabel('Trace #');
    else
        hIm=imagesc(X,1:size(Data,1),Data,max(max(abs(Data)))*[-1,1]);
        xlabel('Horizontal Distance');
    end;
    ylabel('Time Samples');
end;
hCB=InteractiveColorbar;hCB.Sym=true;