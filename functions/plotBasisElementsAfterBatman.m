function [easy_plotting,ppm_new] = plotBasisElementsAfterBatman(path_to_results,final_files,increase_ppm_size)
%{
    Chris Esselman 7.15.25
    Edited log

    - 7.25.25 - Making it so do not need final files to make the image

    - 8.24.25 - Need to increase the resolution of the underlying spectra
    to account for the shifts. Added new input

    - 9.23.25 - Also add the wavelet fit to the output
        
    plotBasisElementsAfterBatman
    
        Make datastructure to easily plot the underlying metabolites for
        the fitting
 
    Inputs:
 
    path_to_results -  Path to individual output directory inside of BatmanOutput
                       directory

    final_files   -   Boolean. True if RelCon.txt file has been made yet

    Increase ppm size - If the ppm size is too small, add increases amount
    of resolution. Is a number and will multiply the current ppm size
                        
     
    
%}

%{
    Psuedocode:
    1. Read in each L.txt file in the directory
    2. Read in RelCon.txt, each specFit_i_rr_j.txt, multidata_user_csv
    MultipletsPpmsShifts.txt, and metabolitesListUsed.txt
    3. For each multiplet, make new vectors out of the L.txt. Just find the
    area where it is zero next to the multiplet
    4. Use circshift to shift each multiplet based on the
    multipletppmshifts.txt. Sum each vector then multiply by RelCon value
%}
if final_files == true
    % Read in the L.txt files names
    current_folder = pwd;
    cd(path_to_results);
    L_listing = dir("L_*.txt");
    if isempty(L_listing)
        cd(current_folder)
        chr = ['Incorrect Directory.'' Provide path to batman results' newline 'ex/ "/runBatman/BatmanOutput/14_Jul_17_34_37"'];
        error(chr)
    end
    % Read in Mult data
    w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    mult_data = readtable("./multi_data.dat");

    warning(w);
    % Read in RelCon.txt
    relative_concentrations = readtable("./RelCon.txt");

    % Read in MultipletsPpmShifts.txt
    multiplet_shifts = readtable("./MultipletsPpmShifts.txt");

    % Read in Metabolites used file
    FID = fopen("./metabolitesListUsed.txt");
    data = textscan(FID,'%s');
    fclose(FID);
    metabolites_used = string(data{:});

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

    % Get the L matrix of the fits
    l_matrix = readmatrix(L_listing(1).name);

    % Go back to current folder
    cd(current_folder)

    % Next make a matrix that will hold every multiplet in its individual
    % vector

    % First get the ppm vector
    ppm = s_fits(1).tables.ppm';

    % Preallocate a matrix to hold the original multiplet locations
    original_multiplets = zeros(size(multiplet_shifts,1),length(ppm));

    % Do a loop to add the multiplet vectors in original location
    counter = 0;
    % relative concentrations counter
    relative_concentrations_counter = zeros(size(multiplet_shifts,1),1);
    for i = 1:length(metabolites_used)
        % Find part of table that matches the metabolite
        table_loop = mult_data(strcmp(metabolites_used(i),mult_data.Var1),"Var4");
        for j = 1:height(table_loop)
            counter = counter + 1;
            % Split the string
            indices_mult = strsplit(string(table_loop{j,"Var4"}),',');
            % Get the original from L_matrix
            original_multiplets(counter,:) = l_matrix(:,i)';
            % Make everything outside matrix = 0
            original_multiplets(counter,ppm < str2double(indices_mult(1)) | ppm > str2double(indices_mult(2))) = 0;
            % add to relative concentrations counter
            relative_concentrations_counter(counter) = i;
        end
    end
    % Will contain the structure for all the plotting
    easy_plotting = struct;

    ppm_length = length(ppm);
    ppm_new = linspace(ppm(1),ppm(end),ppm_length*increase_ppm_size);
    % Will tell how far to shift each multiplet
    differences_ppm = diff(ppm_new);
    average_size = abs(mean(differences_ppm));

    for i = 1:size(s_fits,2)
        % Add original spectrum and the underlying multiplets original
        easy_plotting(i).original = interp1(ppm,s_fits(i).tables.OriginalSpectrum',ppm_new);
        easy_plotting(i).metabFit = interp1(ppm,s_fits(i).tables.MetabolitesFit',ppm_new);
        easy_plotting(i).waveletFit = interp1(ppm,s_fits(i).tables.WaveletFit',ppm_new);
        easy_plotting(i).metabolites = zeros(length(metabolites_used),length(ppm_new));
        original_multiplets_loop = zeros(size(original_multiplets,1),length(ppm_new));
        for j = 1:size(original_multiplets,1)
            original_multiplets_loop(j,:) = interp1(ppm,original_multiplets(j,:),ppm_new);
            % shift each multiplet
            shift = floor(multiplet_shifts{j,i+1}/average_size);
            original_multiplets_loop(j,:) = circshift(original_multiplets_loop(j,:),shift*(-1));
        end
        for j = 1:length(metabolites_used)
            % Multiply by the concentration
            easy_plotting(i).metabolites(j,:) = sum(original_multiplets_loop(relative_concentrations_counter==j,:),1)*relative_concentrations{j,i+1};
        end
    end
else
    % Read in the L.txt files names
    current_folder = pwd;
    cd(path_to_results);
    L_listing = dir("L_*.txt");
    if isempty(L_listing)
        cd(current_folder)
        chr = ['Incorrect Directory.'' Provide path to batman results' newline 'ex/ "/runBatman/BatmanOutput/14_Jul_17_34_37"'];
        error(chr)
    end
    % Read in Mult data
    w=warning('off','MATLAB:table:ModifiedAndSavedVarnames');
    mult_data = readtable("./multi_data.dat");

    warning(w);
    % Read in beta files
    B_listing = dir("beta_*_rr_0.txt");
    B_index = logical(zeros(size(B_listing,1),1));
    for i = 1:size(B_listing,1)
        if contains(B_listing(i).name,"sam")
            B_index(i) = true;
        end
    end
    B_listing(B_index) = [];
    [~, reindex] = sort( str2double( regexp( {B_listing.name}, '\d+', 'match', 'once' )));
    B_listing = B_listing(reindex);
    for i = 1:size(B_listing,1)
        if i == 1
            loop_column = readmatrix(B_listing(i).name);
            relative_concentrations = [zeros(length(loop_column),1) loop_column];
        else
            loop_column = readmatrix(B_listing(i).name);
            relative_concentrations = [relative_concentrations loop_column];
        end
    end
    relative_concentrations = array2table(relative_concentrations);
    %relative_concentrations = readtable("./RelCon.txt");

    % Read in MultipletsPpmShifts.txt
    % Read in the specFit_i_rr_j.txt file names
    M_listing = dir("./delta_draw_mean_*.txt");
    [~, reindex] = sort( str2double( regexp( {M_listing.name}, '\d+', 'match', 'once' )));
    M_listing = M_listing(reindex);
    for i = 1:size(M_listing,1)
        if i == 1
            loop_column = readmatrix(M_listing(i).name);
            multiplet_shifts = [zeros(length(loop_column),1) loop_column];
        else
            loop_column = readmatrix(M_listing(i).name);
            multiplet_shifts = [multiplet_shifts loop_column];
        end
    end
    multiplet_shifts = array2table(multiplet_shifts);
    %multiplet_shifts = readtable("./MultipletsPpmShifts.txt");

    % Read in Metabolites used file
    FID = fopen("./metabolitesListUsed.txt");
    data = textscan(FID,'%s');
    fclose(FID);
    metabolites_used = string(data{:});

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

    % Get the L matrix of the fits
    l_matrix = readmatrix(L_listing(1).name);

    % Go back to current folder
    cd(current_folder)

    % Next make a matrix that will hold every multiplet in its individual
    % vector

    % First get the ppm vector
    ppm = s_fits(1).tables.ppm';

    % Preallocate a matrix to hold the original multiplet locations
    original_multiplets = zeros(size(multiplet_shifts,1),length(ppm));

    % Do a loop to add the multiplet vectors in original location
    counter = 0;
    % relative concentrations counter
    relative_concentrations_counter = zeros(size(multiplet_shifts,1),1);
    for i = 1:length(metabolites_used)
        % Find part of table that matches the metabolite
        table_loop = mult_data(strcmp(metabolites_used(i),mult_data.Var1),"Var4");
        for j = 1:height(table_loop)
            counter = counter + 1;
            % Split the string
            indices_mult = strsplit(string(table_loop{j,"Var4"}),',');
            % Get the original from L_matrix
            original_multiplets(counter,:) = l_matrix(:,i)';
            % Make everything outside matrix = 0
            original_multiplets(counter,ppm < str2double(indices_mult(1)) | ppm > str2double(indices_mult(2))) = 0;
            % add to relative concentrations counter
            relative_concentrations_counter(counter) = i;
        end
    end
    % Will contain the structure for all the plotting
    easy_plotting = struct;

    ppm_length = length(ppm);
    ppm_new = linspace(ppm(1),ppm(end),ppm_length*increase_ppm_size);
    % Will tell how far to shift each multiplet
    differences_ppm = diff(ppm);
    average_size = abs(mean(differences_ppm));

    for i = 1:size(s_fits,2)
        % Add original spectrum and the underlying multiplets original
        easy_plotting(i).original = interp1(ppm,s_fits(i).tables.OriginalSpectrum',ppm_new);
        easy_plotting(i).metabFit = interp1(ppm,s_fits(i).tables.MetabolitesFit',ppm_new);
        easy_plotting(i).waveletFit = interp1(ppm,s_fits(i).tables.WaveletFit',ppm_new);
        easy_plotting(i).metabolites = zeros(length(metabolites_used),length(ppm_new));
        original_multiplets_loop = zeros(size(original_multiplets,1),length(ppm_new));
        for j = 1:size(original_multiplets,1)
            original_multiplets_loop(j,:) = interp1(ppm,original_multiplets(j,:),ppm_new);
            % shift each multiplet
            shift = floor(multiplet_shifts{j,i+1}/average_size);
            original_multiplets_loop(j,:) = circshift(original_multiplets_loop(j,:),shift*(-1));
        end
        for j = 1:length(metabolites_used)
            % Multiply by the concentration
            easy_plotting(i).metabolites(j,:) = sum(original_multiplets_loop(relative_concentrations_counter==j,:),1)*relative_concentrations{j,i+1};
        end
    end

end
end
