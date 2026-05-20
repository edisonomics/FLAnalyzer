function optimizeBatmanFits(path_to_results,path_to_chemShiftPerSpec)
% %
% %{
%     Chris Esselman 5.20.26
%     Edited log
%
%     optimizeBatmanFits - Function to update starting positions of
%     multiplets
%
%     Inputs-
%
%         path_to_results - the path to the finished batman run output folder
%          ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"
%
%         path_to_chemShiftPerSpec - the path to chemShiftPerSpec.csv
%         ex/ "~~/runBATMAN/BatmanInput/chemShiftPerSpec.csv"
%
%     Optional Inputs-
%
%         region_sort - 1x2 vector to sort a specific region of the
%
%         ex/ [7.1 7.5] - sort the spectra in descending order of percent
%         quantified for the specific region.
%
%     Outpus-
%
%         Edits chemShiftPerSpec.csv for the next round of optimization
% %}
% %

% So I can pass optional arguments
arguments
    path_to_results
    path_to_chemShiftPerSpec
end

% Read in the L.txt files names
current_folder = pwd;
cd(path_to_results);

% Read in the chemShiftPerSpec.csv
chemShiftPerSpec = readtable(path_to_chemShiftPerSpec);
% Read in Mult data
w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
mult_data = readtable("./multi_data.dat");

warning(w);

% Read in MultipletsPpmShifts.txt
multiplet_shifts = readtable("./MultipletsPpmShifts.txt");

% Read in Metabolites used file
% FID = fopen("./metabolitesListUsed.txt");
% data = textscan(FID,'%s');
% fclose(FID);
% metabolites_used = string(data{:});

% Read in the specFit_i_rr_j.txt file names
s_listing = dir("./specFit_*_rr_0.txt");
[~, reindex] = sort( str2double( regexp( {s_listing.name}, '\d+', 'match', 'once' )));
s_listing = s_listing(reindex);

% Put the L_listing and s_listing into structures
s_fits = struct;
w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
for i = 1:size(s_listing,1)
    s_fits(i).tables = readtable(s_listing(i).name);
end
warning(w);

% Go back to current folder
cd(current_folder)



% First get the ppm vector
ppm = s_fits(1).tables.ppm';

% Read in the original and wavelets
X_wavelet = zeros(size(s_fits,2),size(s_fits(1).tables,1));
X_original = zeros(size(s_fits,2),size(s_fits(1).tables,1));
for i = 1:size(s_fits,2)
    X_wavelet(i,:) = s_fits(i).tables.WaveletFit';
    X_original(i,:) = s_fits(i).tables.OriginalSpectrum';
end

% Do the bucketing of the data
[~,~,chem_interval] = opt_bucket(ppm,X_original,0.2,1);
figure,plot(ppm,X_original)
hold on
xline(chem_interval(:,1))
hold off
[multiplet_shifts, mult_data, rows_e] = sync_datasets(multiplet_shifts, mult_data,chemShiftPerSpec);




end

function [clean_shifts, clean_data, rows_to_edit] = sync_datasets(multiplet_shifts, mult_data, chemShiftPerSpec)
    % SYNC_DATASETS Synchronizes datasets and extracts row indices for chemShiftPerSpec.
    %
    % Inputs:
    %   multiplet_shifts: Table from MultipletsPpmShifts.txt 
    %                     (Assumes Col 1 is the ID like 'metabolite1_2.433')
    %   mult_data:        Table from multi_data.dat 
    %                     (Assumes Col 1 is Name, Col 2 is PPM Shift, Col 6 is User Shift)
    %   chemShiftPerSpec: Table containing 'multiplets' and 'pos_in_ppm' columns
    %
    % Outputs:
    %   clean_shifts:     The filtered and aligned multiplet_shifts variable
    %   clean_data:       The filtered and aligned mult_data variable
    %   rows_to_edit:     Column vector of exact row indices to edit in the original chemShiftPerSpec

    % 1. Extract Keys from multiplet_shifts (The Reference)
    if istable(multiplet_shifts)
        keys1_raw = multiplet_shifts{:, 1};
    else
        keys1_raw = multiplet_shifts(:, 1);
    end
    keys1 = cellstr(string(keys1_raw));

    % 2. Construct Keys for mult_data
    if istable(mult_data)
        names2 = mult_data{:, 1};
        ppms_default = mult_data{:, 2};
        ppms_user = mult_data{:, 6};
    else
        names2 = mult_data(:, 1);
        if iscell(mult_data(:, 2))
            ppms_default = cell2mat(mult_data(:, 2));
            ppms_user = cell2mat(mult_data(:, 6));
        else
            ppms_default = mult_data(:, 2);
            ppms_user = mult_data(:, 6);
        end
    end

    keys2 = cell(length(names2), 1);
    for i = 1:length(names2)
        if ppms_user(i) ~= -50
            active_ppm = ppms_user(i);
        else
            active_ppm = ppms_default(i);
        end
        keys2{i} = sprintf('%s_%.3f', string(names2{i}), active_ppm);
    end

    % 3. Construct Keys for chemShiftPerSpec
    if istable(chemShiftPerSpec)
        names3 = chemShiftPerSpec.multiplets;
        ppms3 = chemShiftPerSpec.pos_in_ppm;
    else
        names3 = chemShiftPerSpec(:, 1);
        if iscell(chemShiftPerSpec(:, 2))
            ppms3 = cell2mat(chemShiftPerSpec(:, 2));
        else
            ppms3 = chemShiftPerSpec(:, 2);
        end
    end

    keys3 = cell(length(names3), 1);
    for i = 1:length(names3)
        keys3{i} = sprintf('%s_%.3f', string(names3{i}), ppms3(i));
    end

    % 4. Find the Intersections across all three datasets
    % Match multiplet_shifts with mult_data
    [~, idx1_2, idx2] = intersect(keys1, keys2, 'stable');
    
    % Match multiplet_shifts with chemShiftPerSpec
    % original_chemShift_rows stores the row numbers in the untouched chemShiftPerSpec file
    [~, idx1_3, original_chemShift_rows] = intersect(keys1, keys3, 'stable');

    % Find the common master indices that exist in ALL THREE datasets
    [common_idx1, iA, iB] = intersect(idx1_2, idx1_3, 'stable');

    % 5. Filter the reference variables and extract the rows to edit
    clean_shifts = multiplet_shifts(common_idx1, :);
    clean_data = mult_data(idx2(iA), :);
    
    % Extract the specific row numbers to pass back to the parent function
    rows_to_edit = original_chemShift_rows(iB);

end