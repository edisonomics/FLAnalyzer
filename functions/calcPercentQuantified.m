function percent_quant = calcPercentQuantified(path_to_results)
    %{
    Chris Esselman 12.5.25
    Edited log
        
    calcPercentQuantified - Function to calculate the percent of the
    original spectrum quantified by the BATMAN fit
    
    Inputs- 

        path_output - the path to the finished batman run output folder
         ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"

    Output - 

        percent_quant - 1xn vector containing the percent quantified for
        each spectra fit in a single BATMAN run
    %}

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
% Calculate the percent quantified
percent_quant = zeros(1,size(fits,2));
for i = 1:size(fits,2)
    percent_quant(i) = (1-(sum(fits(i).tables.WaveletFit.^2)/sum(fits(i).tables.OriginalSpectrum.^2)))*100;
end