function [basis_final,multi_data,multiplet_struct] = makeBatmanOutputMultipletsFL(X,ppm,all_peak_tables,basis,ends_ppm,where_to_save_files,mag_freq_mHz,options,include_standard)
%{
    Chris Esselman 6.11.25
    Edited log
    - 7.10.25 - Add functionality to make basis element of DSS
    
    - 7.11.25 - Subtracted each basis element by the min of the element
    to make sure the baseline is at zero

    - 7.25.25 - Fixed bug on line 188. Linear index was wrong when picking
    largest fraction that has the most peaks. Also changed multiplet center
    from median to mean

    - 7.29.25 - Fixed so the center of the multiplet is the actual center
    and not an offset due to end regions of the peak

    - 8.5.25 - Fixed bug where on line 269 where it was ppm and it should
    have been ppm_full
    

    makeBatmanOutputFL
        Finalize basis set and output files for BATMAN using rough
        multiplet estimation
 
    Inputs:
 
    X                          Data matrix 

    ppm                        Chemical shift vector

    all_peak_tables            structure- fields below
                               1. each_table - contains the SAND csvs
                               2. ft_peaks - contains all the SAND
                               underlying peaks fourier transformed

    basis                      the finalized basis set

    ends_ppm                   vector with two elements. ppm values to were
                               used when creating the basis set

    where_to_save_files        path to folder where all the batman input
                               files will be saved

    mag_freq_mHz               frequency of NMR instrument to collect data
                               in mHz

    options                    different methods to make each basis element
                               'max' - makes basis element from fraction
                               with largest inensity
                               'least_peaks' - makes basis element from
                               fraction with least number of overall peaks
                               'common_ppm' - makes basis element from
                               fraction that has the most peaks in the
                               basis element
  
    include_standard           true or false that if true will find the sand peak
                               that is closest to 0.0 ppm and will add it 
                               to the basis elements
    Outputs:

    basis_final               matrix containing spectra of the final basis
                              elements

    mulit-data                Table in matlab, but is main csv file needed
                              for batman

    mulitplet_struct          structure containing the center and end ppm
                              values for each multiplet in the basis

    batman_input_files                

      
%}
% Find the mins of the data and get ready to index into sand peak tables
% Cut down the X and ppm to match the excluded DSS region
ppm_full = ppm;
X = X(:,ppm < ends_ppm(2) & ppm > ends_ppm(1));
ppm = ppm(:,ppm < ends_ppm(2) & ppm > ends_ppm(1));
mins = islocalmin(X,2);
maxs = islocalmax(X,2);
end_min_idx = matchPPMs(ends_ppm,ppm);
mins(:,ppm < ends_ppm(1) | ppm > ends_ppm(2)) = 0;
mins(:,end_min_idx) = 1;
index_all_peak_tables = struct;
for i = 1:size(all_peak_tables,2)
    index_all_peak_tables(i).indexers = zeros(1,length(all_peak_tables(i).ft_peaks));
end
% Cut down basis to just include the basis elements that were added
if isfield(basis,"removed_traces")
    basis = rmfield(basis,"removed_traces");
end
idx1 = zeros(1,size(basis,2));
idx2 = 1:size(basis,2);
for i = 1:size(basis,2)
    if isempty(basis(i).correlated_traces)
        idx1(i) = 1;
    end
end
basis(idx2(logical(idx1))) = [];

basis_final = zeros(size(basis,2),length(all_peak_tables(1).ft_peaks(1).ft));

% Get the current working directory and then cd to the directory where the
% batman output will be put
currentFolder = pwd;
cd(where_to_save_files);
mkdir PureSpectraTemplate

% Make a variable that will be used to write a file with all the variables
metabolites = strings(size(basis,2),1);

% Get the multi_data_user.csv file ready
multi_data = table;

% Calculate the ppm value to consider multiplet based on magnet frequency
multiplet_ppm = (10/(mag_freq_mHz*1000000))*10^6;

% Struct containing the multiplets
multiplet_struct = struct;
for i = 1:size(basis,2)
    % Find the unique fractions and the unique ppms values for each
    % fraction

    % First find the unique fractions
    unique_fracs = 0;
    unique_ppms = 0;
    compare_max = struct;
    for j = 1:size(basis(i).correlated_traces,2)
        if j == 1
            unique_fracs = basis(i).correlated_traces(j).fracs;
            unique_ppms = basis(i).correlated_traces(j).ppms;
        else
            unique_fracs = [unique_fracs basis(i).correlated_traces(j).fracs];
            unique_ppms = [unique_ppms basis(i).correlated_traces(j).ppms];
        end
    end
    [unique_fracs,~,ic] = unique(unique_fracs);
    for j = 1:length(unique_fracs)
        compare_max(j).ppms = unique_ppms(ic == j);
        compare_max(j).ppms = unique(compare_max(j).ppms);
        % Get rid of ppm values that are not maxes
        compare_max(j).ppms(maxs(unique_fracs(j),compare_max(j).ppms) == 0) = [];
    end
    % Next find the fractions that have no peaks
    indexer_unique_fracs = zeros(1,length(unique_fracs));
    for j = 1:size(compare_max,2)
        if isempty(compare_max(j).ppms)
            indexer_unique_fracs(j) = 1;
        end
    end
    unique_fracs(logical(indexer_unique_fracs)) = [];
    compare_max(logical(indexer_unique_fracs)) = [];

    % Next build each fraction and see which one is the biggest
    sums_fracs = zeros(length(unique_fracs),length(all_peak_tables(1).ft_peaks(1).ft));
    for j = 1:length(unique_fracs)
        for k = 1:length(compare_max(j).ppms)
            %Find the index of the peak in the ppm
            max_idx = compare_max(j).ppms(k);
            %Locate the nearest minima surrounding the peak
            min1_idx = 0;
            min2_idx = 0;
            for l = max_idx:-1:1
                if mins(unique_fracs(j),l) == 1
                    min1_idx = l;
                    break;
                end
            end
            for l = max_idx:length(ppm)
                if mins(unique_fracs(j),l) == 1
                    min2_idx = l;
                    break;
                end
            end
            % new structure field will be a structure. Each row corresponds to
            % each fraction the peak is in
            indexer2 = 1:size(all_peak_tables(unique_fracs(j)).each_table,1);
            multiple_sand = indexer2(all_peak_tables(unique_fracs(j)).each_table.freq_ppm > ppm(min1_idx) & all_peak_tables(unique_fracs(j)).each_table.freq_ppm < ppm(min2_idx));
            for l = 1:length(multiple_sand)
                sums_fracs(j,:) = sums_fracs(j,:) + all_peak_tables(unique_fracs(j)).ft_peaks(multiple_sand(l)).ft;
            end
        end
    end
    if strcmp(options,'max')
        % Find the fraction that has the largest value and add to basis_final
        [~,lin_idx] = max(sums_fracs,[],"all","linear");
        [max_frac,~] = ind2sub(size(sums_fracs),lin_idx);
        basis_final(i,:) = sums_fracs(max_frac,:) - min(sums_fracs(max_frac,:));
    elseif strcmp(options,'least_peaks')
        number_of_peaks = zeros(1,length(unique_fracs));
        X_norm = (X - min(X(:)))/(max(X(:))-min(X(:)));
        for j = 1:length(unique_fracs)
            peak_number = findpeaks(X_norm(unique_fracs(j),:),'MinPeakHeight',0.001);
            number_of_peaks(j) = length(peak_number);
        end
        [~,min_peak_idx] = min(number_of_peaks);
        basis_final(i,:) = sums_fracs(min_peak_idx,:) - min(sums_fracs(min_peak_idx,:));
    elseif strcmp(options,'common_ppm')
        % Start with just having the one with the most ppm values
        max_length = 0;
        max_length_idx = 1;
        for j = 1:length(unique_fracs)
            if length(compare_max(j).ppms) > max_length
                max_length = length(compare_max(j).ppms);
                max_length_idx = j;
            elseif length(compare_max(j).ppms) == max_length
                % Pick the one that is larger intensity
                [~,lin_idx] = max(sums_fracs([max_length_idx,j],:),[],"all","linear");
                [max_frac,~] = ind2sub(size(sums_fracs([max_length_idx,j],:)),lin_idx);
                if max_frac == 2
                    max_length_idx = j;
                end
            end
        end
        basis_final(i,:) = sums_fracs(max_length_idx,:) - min(sums_fracs(max_length_idx,:));
    else
        error('Last input must be: "max","least_peaks",or "common_ppm"')
    end

    % Normalize each basis element
    % Make in between 0 and 1
    %norm_base_1 =  basis_final(i,:)*(1/max(basis_final(i,:)));
    %norm_base_2 = norm_base_1 - mean(norm_base_1(ppm_full>ends_ppm(2)));
    %norm_base_2 = norm_base_1 - min(norm_base_1);
    % Now normalize the integral to 1
    %norm_base_1 = basis_final(i,:) - min(basis_final(i,:));
    norm_base_1 = basis_final(i,:)-min(basis_final(i,:));
    norm_base_2 = norm_base_1/trapz(norm_base_1);
    % Make the puretemplatespectra
    pure_spectra = strings(length(all_peak_tables(1).ft_peaks(1).ft),2);
    pure_spectra(:,1) = compose('%.6f',ppm_full');
    pure_spectra(:,2) = compose('%.6f',norm_base_2');
    writematrix(pure_spectra,sprintf('PureSpectraTemplate/metabolite%d.txt',i),"Delimiter",'tab')

    % For writing the metaboliteList.csv
    metabolites(i) = sprintf('metabolite%d',i);

    % Get the multi_data_user.csv file ready
    labels = ["Metabolite" "pos_in_ppm" "couple_code" "J_constant" "relative_intensity" "overwrite_pos" "overwrite_truncation" "Include_multiplet"];
    var_types = ["string" "string" "double" "string" "string" "string" "string" "double"];

    % If this does not work, could put the peaks in arbitrary locations and
    % then make the raster muliplets based on that?
    [pks,locs] = findpeaks(norm_base_2);
    mins_loop = islocalmin(norm_base_2);
    % Make a 2D matrix the size of the peaks
    cluster_matrix = zeros(length(locs),length(locs));
    for j = 1:length(locs)
        for k = 1:length(locs)
            % Within 10hz or less than 100 percentage difference
            if abs(ppm_full(locs(j))-ppm_full(locs(k))) <= multiplet_ppm && ((abs(pks(j)-pks(k)))/((pks(j)+pks(k))/2))*100 < 100
                cluster_matrix(j,k) = 1;
            end
        end
    end
    S = cluster_matrix;
    r = fliplr(symrcm(cluster_matrix));
    C = {r(1)};
    for j = 2:numel(r)
        if any(S(C{end}, r(j)))
            C{end}(end+1) = r(j);
        else
            C{end+1} = r(j);
        end
    end
    multi_data_loop = table('VariableNames',labels,'Size',[length(C) 8],'VariableTypes',var_types);
    for j = 1:length(C)
        % Next find the end regions of each multiplet
        downfield_peak = max(locs(cell2mat(C(j))));
        upfield_peak = min(locs(cell2mat(C(j))));
        min1_idx = 0;
        min2_idx = 0;
        for l = upfield_peak:-1:1
            if mins_loop(l) == 1 || l == matchPPMs(ppm_full(upfield_peak)-0.03,ppm_full)
                min1_idx = l;
                break;
            end
        end
        for l = downfield_peak:length(ppm_full)
            if mins_loop(l) == 1 || l == matchPPMs(ppm_full(downfield_peak)+0.03,ppm_full)
                min2_idx = l;
                break;
            end
        end
        multiplet_struct(i).centers(j) = mean(ppm_full(locs(cell2mat(C(j)))))-(mean(ppm_full(locs(cell2mat(C(j)))))-mean([ppm_full(min1_idx) ppm_full(min2_idx)]));
        multiplet_struct(i).ends{j} = [ppm_full(min1_idx) ppm_full(min2_idx)];
        multi_data_loop.J_constant(j) = sprintf('%.4f,%.4f',ppm_full(min1_idx),ppm_full(min2_idx));
        multi_data_loop.couple_code(j) = -2;
        multi_data_loop.Include_multiplet(j) = 1;
        multi_data_loop.overwrite_pos(j) = 'n';
        multi_data_loop.overwrite_truncation(j) = 'n';
        multi_data_loop.pos_in_ppm(j) = sprintf('%.4f',multiplet_struct(i).centers(j));
        multi_data_loop.Metabolite(j) = sprintf('metabolite%d',i);
        multi_data_loop.relative_intensity(j) = sprintf('%.4f',trapz(norm_base_2(ppm_full > ppm_full(min1_idx) & ppm_full < ppm_full(min2_idx))));
    end
    if i == 1
        multi_data = multi_data_loop;
    else
        multi_data = [multi_data;multi_data_loop];
    end
end
% Add Dss to end of basis final to check if the ppm is right
% Will add the sand peak that is closest to 0
%basis_final
if include_standard == true
    % Make a table that we will add to the end of the other table
    labels = ["Metabolite" "pos_in_ppm" "couple_code" "J_constant" "relative_intensity" "overwrite_pos" "overwrite_truncation" "Include_multiplet"];
    var_types = ["string" "string" "double" "string" "string" "string" "string" "double"];
    multi_data_loop = table('VariableNames',labels,'Size',[1 8],'VariableTypes',var_types);
    % Find which Dss is closest to zero
    stand_idx_row = 0;
    stand_idx_sand_full = 0;
    closest_dss = inf;
    for i = 1:size(all_peak_tables,2)
        [value,stand_idx_sand] = min(abs(all_peak_tables(i).each_table.freq_ppm));
        if abs(value) < abs(closest_dss)
            closest_dss = value;
            stand_idx_row = i;
            stand_idx_sand_full = stand_idx_sand;
        end
    end
    % Add Dss to basis final
    basis_final(end+1,:) = all_peak_tables(stand_idx_row).ft_peaks(stand_idx_sand_full).ft - min(all_peak_tables(stand_idx_row).ft_peaks(stand_idx_sand_full).ft);
    % Normalize so the integral is 1
    norm_base_2 = basis_final(end,:)/trapz(basis_final(end,:));
    % Get index to add to end of the basis set
    size_ms = size(multiplet_struct,2)+1;
    % Add to multiplet struct
    multiplet_struct(size_ms).centers(1) = closest_dss-(closest_dss-0);
    multiplet_struct(size_ms).ends{1} = [-0.04 0.04];
    % Add to the multi_data
    multi_data_loop.J_constant(1) = sprintf('%.4f,%.4f',-0.04,0.04);
    multi_data_loop.couple_code(1) = -2;
    multi_data_loop.Include_multiplet(1) = 1;
    multi_data_loop.overwrite_pos(1) = 'n';
    multi_data_loop.overwrite_truncation(1) = 'n';
    multi_data_loop.pos_in_ppm(1) = sprintf('%.4f',multiplet_struct(size_ms).centers(1));
    multi_data_loop.Metabolite(1) = sprintf('metabolite%d',size_ms);
    multi_data_loop.relative_intensity(1) = sprintf('%.4f',trapz(norm_base_2(ppm_full > -0.04 & ppm_full < 0.04)));
    multi_data = [multi_data;multi_data_loop];
    % Add the pure spectra
    pure_spectra = strings(length(all_peak_tables(1).ft_peaks(1).ft),2);
    pure_spectra(:,1) = compose('%.6f',ppm_full');
    pure_spectra(:,2) = compose('%.6f',norm_base_2');
    writematrix(pure_spectra,sprintf('PureSpectraTemplate/metabolite%d.txt',size_ms),"Delimiter",'tab')
    % Add metabolite list
    metabolites(size_ms) = sprintf('metabolite%d',size_ms);
end
    
% Write the metaboliteList.csv
writematrix(metabolites,'metabolitesList.csv')
writetable(multi_data,'multi_data_user.csv')
cd(currentFolder)

end
