function plotTracesFL(X,ppm,driver_struct,correlated_traces)
%{
    Chris Esselman 4.29.25
    
    plotTracesFL
        Create a 3D stacked plot of fraction libraries with traces for
        creating a basis set on there

    Inputs:
    
    X                   Data matrix ex/ 140x32k double

    ppm                 ppm values corresponding to data matrix ex/ 1x32k

    driver_struct       struct containing ppms and fracs. Can be created by
                        using the tracePeaksFL function

    Optional Inputs:

    correlated_traces   struct containing ppms of correlated traces to
                        driver

    
        

%}

% So I can pass optional arguments
arguments
    X double
    ppm (1,:) double
    driver_struct
    correlated_traces = 1
end


% Make so values between 0 and 1
X=(X-min(X(:)))/(max(X(:))-min(X(:)));
if ~isstruct(correlated_traces)
    % Plot how the tracing did
    x = zeros(size(X,1),length(ppm));
    for i = 1:size(X,1)
        x(i,:) = ppm;
    end
    y = X;
    for i = 1:size(X,1)
        z(i,:) = (zeros(1,length(ppm)) + i);
    end
    colors = winter(size(X,1));
    fig = figure('Position',[1000 10 900 750]);
    ax1 = axes(fig);
    hold on
    for i = 1:size(X,1)
        plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
    end
    view([0 45])
    for i = 1:size(driver_struct,2)
        plot3(ppm(driver_struct(i).ppms),driver_struct(i).fracs,diag(X(driver_struct(i).fracs,driver_struct(i).ppms)),'LineWidth',1.5,'Color','m')
    end
    hold off
    set(gca,'XDir','rev');
    originalZLim = zlim(ax1);
    ax1.Clipping = "off";
    z_data_num = 1;
else
    for i = 1:size(correlated_traces,2)
        if length(correlated_traces(i).ppms) == length(driver_struct.ppms)
            if correlated_traces(i).ppms == driver_struct.ppms
                correlated_traces(i) = [];
                break
            end
        end
    end
    % Plot how the tracing did
    x = zeros(size(X,1),length(ppm));
    for i = 1:size(X,1)
        x(i,:) = ppm;
    end
    y = X;
    for i = 1:size(X,1)
        z(i,:) = (zeros(1,length(ppm)) + i);
    end
    colors = winter(size(X,1));
    fig = figure('Position',[1000 10 900 750]);
    ax1 = axes(fig);
    hold on
    for i = 1:size(X,1)
        plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
    end
    view([0 45])
    for i = 1:size(driver_struct,2)
        plot3(ppm(driver_struct(i).ppms),driver_struct(i).fracs,diag(X(driver_struct(i).fracs,driver_struct(i).ppms)),'LineWidth',1.5,'Color','m')
    end
    for i = 1:size(correlated_traces,2)
        plot3(ppm(correlated_traces(i).ppms),correlated_traces(i).fracs,diag(X(correlated_traces(i).fracs,correlated_traces(i).ppms)),'LineWidth',1.5,'Color',[0.929 0.694 0.125])
    end
    hold off
    set(gca,'XDir','rev');
    originalZLim = zlim(ax1);
    ax1.Clipping = "off";
    z_data_num = size(correlated_traces,2) + 1;
end
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

end