function plotBasisFL(X_original,ppm_original,basis,sand_end_regions_exclude_dss)
%{
    Chris Esselman 5.17.25
    
    plotBasisFL
        Create a 3D stacked plot of fraction libraries with traces of basis
        set on there

    Inputs:
    
    X_original          Data matrix of orgininal data ex/ 140x32k double

    ppm_original        ppm values corresponding to data matrix ex/ 1x32k

    basis               struct obtained from creating basis set a struct of
                        structs containing correlated traces
                        

    
        

%}

% Cut down the ends of the matrix to get rid of dss
X_original = X_original(:,ppm_original < sand_end_regions_exclude_dss(2) & ppm_original > sand_end_regions_exclude_dss(1));
ppm_original = ppm_original(:,ppm_original < sand_end_regions_exclude_dss(2) & ppm_original > sand_end_regions_exclude_dss(1));
% Make so values between 0 and 1
X_original=(X_original-min(X_original(:)))/(max(X_original(:))-min(X_original(:)));
% Plot how the tracing did
x = zeros(size(X_original,1),length(ppm_original));
for i = 1:size(X_original,1)
    x(i,:) = ppm_original;
end
y = X_original;
for i = 1:size(X_original,1)
    z(i,:) = (zeros(1,length(ppm_original)) + i);
end
colors = winter(size(X_original,1));
fig = figure('Position',[1000 10 900 750]);
ax1 = axes(fig);
hold on
for i = 1:size(X_original,1)
    plot3(x(i,:),z(i,:),y(i,:),'Color',colors(i,:),'LineWidth',0.5)
end
view([0 45])
total_traces = 0;
colors2 = lines(size(basis,2));
for j = 1:size(basis,2)
    total_traces = total_traces + size(basis(j).correlated_traces,2);
    for i = 1:size(basis(j).correlated_traces,2)
        plot3(ppm_original(basis(j).correlated_traces(i).ppms),basis(j).correlated_traces(i).fracs,diag(X_original(basis(j).correlated_traces(i).fracs,basis(j).correlated_traces(i).ppms)),'LineWidth',1.5,'Color',colors2(j,:))
    end
end
hold off
set(gca,'XDir','rev');
originalZLim = zlim(ax1);
ax1.Clipping = "off";
z_data_num = total_traces + 1;

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