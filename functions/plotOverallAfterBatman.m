function plotOverallAfterBatman(path_to_results,which_fit)
    %{
    Chris Esselman 8.28.25
    Edited log
        
    plotOverallAfterBatman - Function to inspect BATMAn fits. Will plot the
    original spectrum, metabolite fit, and wavelet fit
    
    Inputs- 

        path_output - the path to the finished batman run output folder
         ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"

    Optional-
        which_fit - an integer saying which fit to visualize. The default
        is 1. If fitting more than one spectrum, can change number
         ex/ If fit 10 spectra, an input of 8 will visualize the eighth
         spectra
    %}

% So I can pass optional arguments
arguments
    path_to_results
    which_fit = 1
end

% Go to directory where output is
current_folder = pwd;
cd(path_to_results);

% Read in the specFit_i_rr_j.txt file names
s_listing = dir("./specFit_*_rr_0.txt");
[~, reindex] = sort( str2double( regexp( {s_listing.name}, '\d+', 'match', 'once' )));
s_listing = s_listing(reindex);

if size(s_listing,1) < which_fit
    cd(current_folder)
    error('Error. \nThe second argument is greater than the number of spectra fit')
end
% Put the L_listing and s_listing into structures
fits = struct;
w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
for i = 1:size(s_listing,1)
    fits(i).tables = readtable(s_listing(i).name);
end
warning(w);

% Go back to current folder
cd(current_folder)
X_metab_fit = zeros(size(fits,2),size(fits(1).tables,1));
X_ref_batman = zeros(size(fits,2),size(fits(1).tables,1));
X_wavelet = zeros(size(fits,2),size(fits(1).tables,1));
ppm_fits = fits(1).tables.ppm;
for i = 1:size(fits,2)
    X_metab_fit(i,:) = fits(i).tables.MetabolitesFit';
    X_ref_batman(i,:) = fits(i).tables.OriginalSpectrum';
    X_wavelet(i,:) = fits(i).tables.WaveletFit';
end
figure
hold on
plot(ppm_fits,X_ref_batman(which_fit,:),"k",'LineWidth',1.5)
plot(ppm_fits,X_metab_fit(which_fit,:),"r",'LineWidth',1.5)
plot(ppm_fits,X_wavelet(which_fit,:),"g",'LineWidth',1.5)
hold off
set(gca,'XDir','reverse');
set(gca,'ytick',[])
% See the other fits
figure
tiledlayout(3,1)
ax1 = nexttile;
plot(ppm_fits,X_ref_batman(which_fit,:),"k",'LineWidth',1.5)
set(gca,'XDir','reverse');
set(gca,'ytick',[])
set(gca,'xtick',[])
ax2 = nexttile;
plot(ppm_fits,X_metab_fit(which_fit,:),"r",'LineWidth',1.5)
set(gca,'XDir','reverse');
set(gca,'ytick',[])
set(gca,'xtick',[])
ax3 = nexttile;
plot(ppm_fits,X_wavelet(which_fit,:),"g",'LineWidth',1.5)
set(gca,'XDir','reverse');
set(gca,'ytick',[])
linkaxes([ax1 ax2 ax3],'xy')
end