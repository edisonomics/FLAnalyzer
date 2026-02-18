function sort_order = plotBatmanWavelets(path_to_results,sort_percentQ,region_sort)
% %
% %{
%     Chris Esselman 2.18.26
%     Edited log
% 
%     plotBatmanWavelets - Function to view residuals of many batman fits
% 
%     Inputs- 
% 
%         path_to_results - the path to the finished batman run output folder
%          ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"
% 
%     Optional Inputs-
% 
%         sort_percentQ - a bool.
%             True - Sort residuals in descending order of percent quantified
%             False - Plot residuals in numerical order
% 
%         region_sort - 1x2 vector to sort a specific region of the
%         residuals.sort_percetQ must also be true
% 
%         ex/ [7.1 7.5] - sort the spectra in descending order of percent
%         quantified for the specific region. 
% 
%     Outpus-
% 
%         sort_order - a vector containing the order of the sorting
% %}
% %

% So I can pass optional arguments
arguments
    path_to_results
    sort_percentQ (1,1) logical = false
    region_sort double {mustBeVector, mustBeTwoElements} = [Inf, Inf]
end

% Go to directory where output is
current_folder = pwd;
cd(path_to_results);

% Read in the specFit_i_rr_j.txt file names
s_listing = dir("./specFit_*_rr_0.txt");
[~, reindex] = sort( str2double( regexp( {s_listing.name}, '\d+', 'match', 'once' )));
s_listing = s_listing(reindex);

% Put the L_listing and s_listing into structures
fits = struct;
w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
for i = 1:size(s_listing,1)
    fits(i).tables = readtable(s_listing(i).name);
end
warning(w);

% Go back to current folder
cd(current_folder)
X_wavelet = zeros(size(fits,2),size(fits(1).tables,1));
ppm_fits = fits(1).tables.ppm;
for i = 1:size(fits,2)
    X_wavelet(i,:) = fits(i).tables.WaveletFit';
end
sort_order = (1:size(X_wavelet,1))';
if sort_percentQ == false

    annotation_fraction_plot(X_wavelet,ppm_fits)

elseif sort_percentQ == true && all(isinf(region_sort))

    percent_quant = zeros(1,size(fits,2));
    for i = 1:size(fits,2)
        percent_quant(i) = (1-(sum(fits(i).tables.WaveletFit.^2)/sum(fits(i).tables.OriginalSpectrum.^2)))*100;
    end
    [~,sort_order] = sort(percent_quant,'descend');
    sort_order = sort_order';
    X_wavelet_sort = X_wavelet(sort_order,:);
    annotation_fraction_plot(X_wavelet_sort,ppm_fits)

elseif sort_percentQ == true && ~all(isinf(region_sort))

    region_sort = sort(region_sort(:).');
    percent_quant = zeros(1,size(fits,2));
    for i = 1:size(fits,2)
        wave = fits(i).tables.WaveletFit(ppm_fits > region_sort(1) & ppm_fits < region_sort(2));
        og = fits(i).tables.OriginalSpectrum(ppm_fits > region_sort(1) & ppm_fits < region_sort(2));
        percent_quant(i) = (1-(sum(wave.^2)/sum(og.^2)))*100;
    end
    [~,sort_order] = sort(percent_quant,'descend');
    sort_order = sort_order';
    X_wavelet_sort = X_wavelet(sort_order,:);
    annotation_fraction_plot(X_wavelet_sort,ppm_fits)
    % hold on
    % x1 = zeros(1,size(fits,2)) + region_sort(1);
    % x2 = zeros(1,size(fits,2)) + region_sort(2);
    % y = 1:size(fits,2);
    % z = zeros(1,size(fits,2));
    % plot3(x1,y,z,'Color',"red",'LineWidth',2)
    % plot3(x2,y,z,'Color',"red",'LineWidth',2)
    % hold off
end
end
function mustBeTwoElements(v)
    if numel(v) ~= 2
        error("Value must have exactly 2 elements.");
    end
end