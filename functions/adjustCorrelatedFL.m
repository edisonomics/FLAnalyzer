function adjustCorrelatedFL(X,ppm,driver_struct,correlated_traces)
%{
    Chris Esselman 4.30.25
    
    adjustCorrelatedFL
        Create a 3D stacked plot of fraction libraries with traces for
        creating a basis set on there. Have buttons for adding and deleting
        traces

    Inputs:
    
    X                   Data matrix ex/ 140x32k double

    ppm                 ppm values corresponding to data matrix ex/ 1x32k

    driver_struct       struct containing ppms and fracs. Can be created by
                        using the tracePeaksFL function

    correlated_traces   struct containing ppms of correlated traces to
                        driver      

%}

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
fig = figure('Position',[10 10 1500 900]);
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

% Doing stuff to make the restore view right
tb = axtoolbar(ax1,{'datacursor','rotate','pan','zoomin','zoomout'});
btn = axtoolbarbtn(tb, 'push');
btn.Icon = 'restoreview';
btn.ButtonPushedFcn = @(src, event)restoreview(ax1);


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
uicontrol('Style', 'pushbutton', 'String', 'Delete Trace', ...
    'Position', [50, 125, 120, 20], ...
    'Callback', @(src, event)selectRegion(X,ppm,driver_struct,correlated_traces));

% Create a button to add traces to existing trace
uicontrol('Style', 'pushbutton', 'String', 'Add to Existing Trace', ...
    'Position', [200, 125, 150, 20], ...
    'Callback', @(src, event)add_traces(X,ppm,driver_struct,correlated_traces));

% Create a button to add a new trace
uicontrol('Style', 'pushbutton', 'String', 'Add New Trace', ...
    'Position', [50, 100, 120, 20], ...
    'Callback', @(src, event)add_trace_new(X,ppm,driver_struct,correlated_traces));

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
    function selectRegion(X,ppm,driver_struct,correlated_traces)
        fig2 = gcf;
        d = datacursormode(fig2);
        vals = getCursorInfo(d);
        for r = 1:size(vals,2)
            index_peak = matchPPMs(vals(r).Position(1),ppm);
            closest_in_correlated = 0;
            closest_in_correlated_num = inf;
            for g = 1:size(correlated_traces,2)
                for d = 1:length(correlated_traces(g).ppms)
                    if abs(index_peak-correlated_traces(g).ppms(d)) < closest_in_correlated_num
                        closest_in_correlated = g;
                        closest_in_correlated_num = abs(index_peak-correlated_traces(g).ppms(d));
                    end
                end
            end
            index = correlated_traces(closest_in_correlated).fracs == vals(r).Position(2);
            correlated_traces(closest_in_correlated).ppms(index) = [];
            correlated_traces(closest_in_correlated).fracs(index) = [];
        end
        % Delete rows that are empty
        correlated_traces(all( cell2mat( arrayfun( @(x) structfun( @isempty, x ), correlated_traces, 'UniformOutput', false ) ), 1 ) ) = [];
        xlim_old = xlim(gca);
        ylim_old = ylim(gca);
        zlim_old = zlim(gca);
        [view_old1,view_old2] = view(gca);
        adjustCorrelatedFL(X,ppm,driver_struct,correlated_traces)
        close(fig2)
        xlim(gca, xlim_old);
        ylim(gca, ylim_old);
        zlim(gca, zlim_old);
        view(gca, [view_old1 view_old2]);  % Optional
        assignin("base","correlated_traces",correlated_traces)
    end
% Function for adding
    function add_traces(X,ppm,driver_struct,correlated_traces)
        fig2 = gcf;
        d = datacursormode(fig2);
        vals = getCursorInfo(d);
        for r = 1:size(vals,2)
            index_peak = matchPPMs(vals(r).Position(1),ppm);
            closest_in_correlated = 0;
            closest_in_correlated_num = inf;
            for g = 1:size(correlated_traces,2)
                for d = 1:length(correlated_traces(g).ppms)
                    if abs(index_peak-correlated_traces(g).ppms(d)) < closest_in_correlated_num
                        closest_in_correlated = g;
                        closest_in_correlated_num = abs(index_peak-correlated_traces(g).ppms(d));
                    end
                end
            end
            if sum(ismember(vals(r).Position(2),correlated_traces(closest_in_correlated).fracs)) == 0
                correlated_traces(closest_in_correlated).fracs(end + 1) = vals(r).Position(2);
                correlated_traces(closest_in_correlated).ppms(end + 1) = index_peak;
                [correlated_traces(closest_in_correlated).fracs,I] = sort(correlated_traces(closest_in_correlated).fracs,'ascend');
                correlated_traces(closest_in_correlated).ppms = correlated_traces(closest_in_correlated).ppms(I);
            end
        end
        xlim_old = xlim(gca);
        ylim_old = ylim(gca);
        zlim_old = zlim(gca);
        [view_old1,view_old2] = view(gca);
        adjustCorrelatedFL(X,ppm,driver_struct,correlated_traces)
        close(fig2)
        xlim(gca, xlim_old);
        ylim(gca, ylim_old);
        zlim(gca, zlim_old);
        view(gca, [view_old1 view_old2]);  % Optional
        assignin("base","correlated_traces",correlated_traces)
    end
    function add_trace_new(X,ppm,driver_struct,correlated_traces)
        fig2 = gcf;
        d = datacursormode(fig2);
        vals = getCursorInfo(d);
        if size(vals,2) > 0
            index_correlated = size(correlated_traces,2);
            for r = 1:size(vals,2)
                correlated_traces(index_correlated+1).fracs(r) = vals(r).Position(2);
                correlated_traces(index_correlated+1).ppms(r) = matchPPMs(vals(r).Position(1),ppm);
            end
            [correlated_traces(index_correlated+1).fracs,I] = sort(correlated_traces(index_correlated+1).fracs,'ascend');
            correlated_traces(index_correlated+1).ppms = correlated_traces(index_correlated+1).ppms(I);
        end
        xlim_old = xlim(gca);
        ylim_old = ylim(gca);
        zlim_old = zlim(gca);
        [view_old1,view_old2] = view(gca);
        adjustCorrelatedFL(X,ppm,driver_struct,correlated_traces)
        close(fig2)
        xlim(gca, xlim_old);
        ylim(gca, ylim_old);
        zlim(gca, zlim_old);
        view(gca, [view_old1 view_old2]);  % Optional
        % % Reset 'Mode' to 'auto' so Home goes back to auto zoom
        assignin("base","correlated_traces",correlated_traces)
    end
    % New callback function for restore view
    function restoreview(ax)
        % Add our own end, which is to reset the axes manually
        set(ax, 'XLimMode', 'auto');
        set(ax, 'YLimMode', 'auto');
        set(ax, 'ZLimMode', 'auto');  % For 3D
        view([0 45])
    end
end