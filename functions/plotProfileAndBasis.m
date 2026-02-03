function plotProfileAndBasis(path_to_NMRdata_txt,path_to_PureSpectraTemplate)
%{
    Chris Esselman 1.30.26
    Edited Log


    plotProfileAndBasis

        Plot a combined figure that shows the profiling study and the basis
        elements. Will show basis elements that have multiplets in the area
        that was zoomed into. Will be used to see aread where the
        multiplets do not align well with the profile and will need
        refinement in BATMAN


    Inputs:

    path_to_NMRdata_txt - Path to the NMRdata.txt file that is created from
                          running makeBatmanFileToFit function

    path_to_PureSpectraTemplate - Path to the PureSprectraTemplate folder
                                  created from running makeBatmanOutput
                                  type functions
%}

% Read in NMRdata.txt
profiling_data = readtable(path_to_NMRdata_txt);

% Read in Basis Pure Spectra Template
filePattern = fullfile(path_to_PureSpectraTemplate, 'metabolite*.txt');
p_listing = dir(filePattern);
[~, reindex] = sort( str2double( regexp( {p_listing.name}, '\d+', 'match', 'once' )));
p_listing = p_listing(reindex);
p_spec = struct;
w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
for i = 1:size(p_listing,1)
    p_spec(i).tables = readtable(p_listing(i).name);
end
warning(w);
% Get the ppm for the profiling data and the basis
ppm_profile = profiling_data.ppm';
ppm_basis = p_spec(1).tables.Var1';

% Get the data into two matrices 
X_profile = zeros(size(profiling_data,2)-1,length(ppm_profile));
X_basis = zeros(length(p_spec),length(ppm_basis));

for i = 1:size(X_profile,1)
    X_profile(i,:) = profiling_data{:,i+1}';
end

for i = 1:size(X_basis,1)
    X_basis(i,:) = p_spec(i).tables.Var2';
end
clear profiling_data p_spec

% Change the range of values to visualize
X_shrunk = X_profile;
ppm_shrunk = ppm_profile;
X_shrunk=(X_shrunk-min(X_shrunk(:)))/(max(X_shrunk(:))-min(X_shrunk(:)));
rows = size(X_shrunk,1);
for i = 1:rows
    x(i,:) = ppm_shrunk;
end
y = X_shrunk;
for i = 1:rows
    z(i,:) = (zeros(1,length(ppm_shrunk)) + i);
end
colors = winter(rows);
fig1 = figure('Units','pixels');
ax1 = axes(fig1);
hold on
for i = 1:rows
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
set(gca,'XDir','rev');
view(0,45)
originalZLim1 = zlim(ax1);
ax1.Clipping = "off";
% Create a slider control for Figure 1 (Independent)
slider1 = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 50, 'Value', 1, ...
    'Position', [100, 50, 220, 20], ...
    'Callback', @(src, event) adjustVerticalZoom(ax1, src, originalZLim1));

% Create a text label
sliderLabel1 = uicontrol('Style', 'text', 'Position', [100, 70, 120, 20], ...
    'String', 'Vertical Zoom');

% Create selection button
selectionButton1 = uicontrol('Style', 'pushbutton', 'String', 'Change Clipping', ...
    'Position', [100, 100, 120, 20], ...
    'Callback', @(src, event) selectRegion(ax1));


% Change the range of values to visualize
X_shrunk = X_basis;
ppm_shrunk = ppm_basis;
X_shrunk=(X_shrunk-min(X_shrunk(:)))/(max(X_shrunk(:))-min(X_shrunk(:)));
rows = size(X_shrunk,1);
for i = 1:rows
    x(i,:) = ppm_shrunk;
end
y = X_shrunk;
for i = 1:rows
    z(i,:) = (zeros(1,length(ppm_shrunk)) + i);
end
colors = lines(rows);
fig2 = figure('Units','pixels');
ax2 = axes(fig2);
hold on
for i = 1:rows
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
set(gca,'XDir','rev');
view(0,45)
originalZLim2 = zlim(ax2);
ax2.Clipping = "off";

% Create a slider control for Figure 2 (Independent)
slider2 = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 50, 'Value', 1, ...
    'Position', [100, 50, 220, 20], ...
    'Callback', @(src, event) adjustVerticalZoom(ax2, src, originalZLim2));

% Create a text label
sliderLabel2 = uicontrol('Style', 'text', 'Position', [100, 70, 120, 20], ...
    'String', 'Vertical Zoom');

% Create selection button
selectionButton2 = uicontrol('Style', 'pushbutton', 'String', 'Change Clipping', ...
    'Position', [100, 100, 120, 20], ...
    'Callback', @(src, event) selectRegion(ax2));

% --- LINKING LOGIC FOR SIZE ONLY ---
fig1.SizeChangedFcn = @(src, event) syncSize(fig1, fig2);
fig2.SizeChangedFcn = @(src, event) syncSize(fig2, fig1);

linkaxes([ax1 ax2],'x')
% Create a listener to synchronize the view
addlistener(ax1, 'View', 'PostSet', @(src,evt) set(ax2, 'View', get(ax1, 'View')));
% Create a listener to synchronize the view
addlistener(ax2, 'View', 'PostSet', @(src,evt) set(ax1, 'View', get(ax2, 'View')));
%% Function Definitions

    function syncSize(source, target)
        % 1. FOCUS CHECK (The Fix)
        % If 'source' is not the currently active figure (gcf), it means 
        % you are not interacting with it (it's being updated by code).
        % We ignore it to prevent the loop/bounce.
        if source ~= gcf
            return;
        end
        
        % 2. Get positions
        sPos = source.Position;
        tPos = target.Position;

        % 3. Check Tolerance (Optimization)
        % Only update if the difference is real (> 1 pixel)
        if max(abs(sPos(3:4) - tPos(3:4))) > 1
            % Update target size, preserve target position
            target.Position(3:4) = sPos(3:4);
        end
    end

    function adjustVerticalZoom(ax, slider, ~)
        zoomFactor = slider.Value;
        zData = get(ax.Children, 'ZData');
        if iscell(zData)
            zData = cell2mat(zData);
        end
        zMax = max(zData(:));
        baseline = -zMax / zoomFactor/3; 
        zlim(ax, [baseline, zMax / zoomFactor]);
    end

    function selectRegion(ax)
        if ax.Clipping == "off"
            ax.Clipping = "on";
        else
            ax.Clipping = "off";
        end
    end
end