function percent_quant = calcPercentQuantified(path_to_results, options)
    % %{
    % calcPercentQuantified - Function to calculate the percent of the
    % original spectrum quantified by the BATMAN fit
    % 
    % Inputs- 
    %   path_to_results - the path to the finished batman run output folder
    % 
    % Name-Value Optional Inputs-
    %   include_region  - 1x2 double array [min_ppm, max_ppm]. 
    %                     Only calculate the fit within this specific window.
    %
    %   exclude_regions - Nx2 double matrix. Each row is a [min_ppm, max_ppm]
    %                     region to exclude from the calculation.
    % 
    % Output - 
    %   percent_quant   - 1xn vector containing the percent quantified
    % %}

    % Allow optional inputs using Name-Value pairs
    arguments
        path_to_results
        options.include_region double = [-Inf, Inf] % Default includes everything
        options.exclude_regions double = []         % Default excludes nothing
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
    w = warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    for i = 1:size(s_listing,1)
        fits(i).tables = readtable(s_listing(i).name);
    end
    warning(w);

    % Go back to current folder
    cd(current_folder)
    
    % Calculate the percent quantified
    percent_quant = zeros(1, size(fits, 2));
    
    for i = 1:size(fits, 2)
        % Extract vectors to keep the code clean
        ppm = fits(i).tables.ppm;
        wave = fits(i).tables.WaveletFit;
        og = fits(i).tables.OriginalSpectrum;

        % 1. Create a baseline mask of all true values
        valid_idx = true(size(ppm));

        % 2. Apply Inclusion Region (Accessed via options.)
        if ~all(isinf(options.include_region))
            inc_bounds = sort(options.include_region); 
            valid_idx = valid_idx & (ppm >= inc_bounds(1) & ppm <= inc_bounds(2));
        end

        % 3. Apply Exclusion Regions (Accessed via options.)
        if ~isempty(options.exclude_regions)
            for j = 1:size(options.exclude_regions, 1)
                ex_bounds = sort(options.exclude_regions(j, :)); 
                valid_idx = valid_idx & ~(ppm >= ex_bounds(1) & ppm <= ex_bounds(2));
            end
        end

        % 4. Filter the data arrays using the finalized logical mask
        wave_filtered = wave(valid_idx);
        og_filtered = og(valid_idx);

        % 5. Calculate percent quantified on the filtered data
        if isempty(og_filtered) || sum(og_filtered.^2) == 0
            percent_quant(i) = NaN;
        else
            percent_quant(i) = (1 - (sum(wave_filtered.^2) / sum(og_filtered.^2))) * 100;
        end
    end
end