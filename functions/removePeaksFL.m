function [X_removed,all_peak_tables_removed,X_added] = removePeaksFL(X,ppm,all_peak_tables,peaks_to_remove_struct,ends_ppm)
%{
    Chris Esselman 4.23.25
    Edited log
        
    removePeaksFL
    
        Remove peaks as they elute over Fraction Library data. Input should
        be peak maxima in the data
 
    Inputs:
 
    X                          Data matrix 

    ppm                        Chemical shift vector

    all_peak_tables            structure- fields below
                               1. each_table - contains the SAND csvs
                               2. ft_peaks - contains all the SAND
                               underlying peaks fourier transformed

    peaks_to_remove_struct     structure - fields below
                               1. ppms - ppm indices of the trace
                               2. fracs - fracs/rows of the trace from
                               matrix

    ends_ppm                   vector with two elements. ppm values where
                               the ends of the spectrum is. The data to
                               consider
  

    Outputs:

    all_peak_tables_removed   Will be the same as all_peak_tables but with
                              the traces removed

    X_removed                 New X_matrix with the traces removed. Will
                              not be normalized

      
%}
all_peak_tables_removed = all_peak_tables;
mins = islocalmin(X,2);
end_min_idx = matchPPMs(ends_ppm,ppm);
mins(:,ppm < ends_ppm(1) | ppm > ends_ppm(2)) = 0;
mins(:,end_min_idx) = 1;
index_all_peak_tables = struct;
sum_over_fractions = zeros(1,length(all_peak_tables(1).ft_peaks(1).ft));
peak_index_each_fraction = struct;
for i = 1:size(all_peak_tables,2)
    index_all_peak_tables(i).indexers = zeros(1,length(all_peak_tables(i).ft_peaks));
    [~, peak_index_each_fraction(i).Peaks_idx] = findpeaks(X(i,:));
    %peak_index_each_fraction(i).ppm = ppm(peak_index_each_fraction(i).Peaks_idx);
end
for j = 1:size(peaks_to_remove_struct,2)
    for k = 1:length(peaks_to_remove_struct(j).fracs)
        %Find the index of the peak in the ppm
        max_idx = peaks_to_remove_struct(j).ppms(k);
        [~,ispeak] = ismember(max_idx,peak_index_each_fraction(peaks_to_remove_struct(j).fracs(k)).Peaks_idx);
        if ispeak > 0
            %Locate the nearest minima surrounding the peak
            min1_idx = 0;
            min2_idx =0;
            for l = max_idx:-1:1
                if mins(peaks_to_remove_struct(j).fracs(k),l) == 1
                    min1_idx = l;
                    break;
                end
            end
            for l = max_idx:length(ppm)
                if mins(peaks_to_remove_struct(j).fracs(k),l) == 1
                    min2_idx = l;
                    break;
                end
            end
            % new structure field will be a structure. Each row corresponds to
            % each fraction the peak is in
            indexer2 = 1:size(all_peak_tables(peaks_to_remove_struct(j).fracs(k)).each_table,1);
            multiple_sand = indexer2(all_peak_tables(peaks_to_remove_struct(j).fracs(k)).each_table.freq_ppm > ppm(min1_idx) & all_peak_tables(peaks_to_remove_struct(j).fracs(k)).each_table.freq_ppm < ppm(min2_idx));
            for l = 1:length(multiple_sand)
                sum_over_fractions = sum_over_fractions + all_peak_tables(peaks_to_remove_struct(j).fracs(k)).ft_peaks(multiple_sand(l)).ft;
            end
            index_all_peak_tables(peaks_to_remove_struct(j).fracs(k)).indexers(multiple_sand) = 1;
        end
    end
end
X_added = sum_over_fractions;
X_removed = zeros(size(X,1),length(all_peak_tables(1).ft_peaks(1).ft));
% Get rid of the peaks in the index
for i = 1:size(all_peak_tables_removed,2)
    % X_loop = zeros(size(all_peak_tables(i).ft_peaks,2),length(all_peak_tables(1).ft_peaks(1).ft));
    % for j = 1:size(X_loop,1)
    %     X_loop(j,:) = all_peak_tables(i).ft_peaks(j).ft;
    % end
    % X_loop(logical(index_all_peak_tables(i).indexers),:) = [];
    % X_removed(i,:) = sum(X_loop);
    all_peak_tables_removed(i).each_table(logical(index_all_peak_tables(i).indexers),:) = [];
    all_peak_tables_removed(i).ft_peaks(logical(index_all_peak_tables(i).indexers)) = [];
    X_loop = zeros(size(all_peak_tables_removed(i).ft_peaks,2),length(all_peak_tables(1).ft_peaks(1).ft));
    for j = 1:size(X_loop,1)
        X_loop(j,:) = all_peak_tables_removed(i).ft_peaks(j).ft;
    end
    X_removed(i,:) = sum(X_loop);
    X_removed(i,:) = X_removed(i,:) - X_removed(i,1);
end
end
