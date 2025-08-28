function adjustDriverFL(X,ppm,driver_struct)
%{
    Chris Esselman 4.29.25
    
    adjustDriverFL
        Create a 3D stacked plot of fraction libraries with traces for
        creating a basis set on there. When asked to close the 

    Inputs:
    
    X                   Data matrix ex/ 140x32k double

    ppm                 ppm values corresponding to data matrix ex/ 1x32k

    driver_struct       struct containing ppms and fracs. Can be created by
                        using the tracePeaksFL function

%}

% So I can pass optional arguments
arguments
    X double
    ppm (1,:) double
    driver_struct
end

% Make the 3D plot
x = zeros(size(X,1),length(ppm));
for i = 1:size(X,1)
    x(i,:) = ppm;
end
y = X;
for i = 1:size(X,1)
    z(i,:) = (zeros(1,length(ppm)) + i);
end
colors = winter(size(X,1));
fig = figure('Position',[10 10 1500 900]);
ax1 = axes(fig);
hold on
% Plot the spectra
for i = 1:size(X,1)
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
view([0 45])
% Plot the trace
for i = 1:size(driver_struct,2)
    plot3(ppm(driver_struct(i).ppms),driver_struct(i).fracs,diag(X(driver_struct(i).fracs,driver_struct(i).ppms)),'LineWidth',1.5,'Color','m','LineStyle','-','Marker','.','MarkerSize',18)
end
hold off
set(gca,'XDir','rev');
originalZLim = zlim(ax1);
ax1.Clipping = "off";
z_data_num = 1;

% Create a slider control
slider = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 50, 'Value', 1, ...
    'Position', [100, 50, 220, 20], ...
    'Callback', @(src, event) adjustVerticalZoom(ax1, src, z_data_num,originalZLim));
% slider = uicontrol('Style', 'slider', 'Min', 0.01, 'Max', 100, 'Value', 1, ...
%                    'Position', [100, 50, 220, 20], ...
%                    'Callback', @(src, event) adjustVerticalZoom(ax, src, originalZLim));

% Create a text label for the slider
sliderLabel = uicontrol('Style', 'text', 'Position', [100, 70, 120, 20], ...
    'String', 'Vertical Zoom');

% Create a button to delete traces
uicontrol('Style', 'pushbutton', 'String', 'Save Deletions', ...
    'Position', [50, 100, 120, 20], ...
    'Callback', @(src, event)selectRegion(X,ppm,driver_struct));

% Create a button to add traces
uicontrol('Style', 'pushbutton', 'String', 'Save Additions', ...
    'Position', [200, 100, 120, 20], ...
    'Callback', @(src, event)add_traces(X,ppm,driver_struct));
%% function defs
    function adjustVerticalZoom(ax, slider, z_data_num ,~)
        % Get the current slider value
        zoomFactor = slider.Value;

        % Adjust the upper limit of the z-axis based on the zoom factor
        zData = get(ax.Children, 'ZData');
        if iscell(zData)
            zData = cell2mat(zData(1+z_data_num:end));
        end
        zMax = max(zData(:));
        baseline = -zMax / zoomFactor/3; % Ensure the baseline is always visible
        zlim(ax, [baseline, zMax / zoomFactor]);
    end
    % Function for deleting 
    function selectRegion(X,ppm,driver_struct)
        fig2 = gcf;
        d = datacursormode(fig2);
        vals = getCursorInfo(d);
        for r = 1:size(vals,2)
            index = driver_struct.fracs == vals(r).Position(2);
            driver_struct.ppms(index) = [];
            driver_struct.fracs(index) = [];
        end
        adjustDriverFL(X,ppm,driver_struct)
        close(fig2)
        assignin("base","driver_struct",driver_struct)
    end
    % Function for adding
    function add_traces(X,ppm,driver_struct)
        fig2 = gcf;
        d = datacursormode(fig2);
        vals = getCursorInfo(d);
        for r = 1:size(vals,2)
            if sum(ismember(vals(r).Position(2),driver_struct.fracs)) == 0
                driver_struct.fracs(end + 1) = vals(r).Position(2);
                driver_struct.ppms(end + 1) = matchPPMs(vals(r).Position(1),ppm);
            end
        end
        [driver_struct.fracs,I] = sort(driver_struct.fracs,'ascend');
        driver_struct.ppms = driver_struct.ppms(I);
        adjustDriverFL(X,ppm,driver_struct)
        close(fig2)
        assignin("base","driver_struct",driver_struct)
    end
end

