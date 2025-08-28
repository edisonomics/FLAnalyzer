function [X_added] = sumPeaksFL(X,ppm,all_peak_tables,peaks_to_add_struct,ends_ppm)
%{
    Chris Esselman 4.24.25
    Edited log
        
    sumPeaksFL
    
        Sum peaks as they elute over Fraction Library data. Input should
        be peak maxima in the data
 
    Inputs:
 
    X                          Data matrix

    ppm                        Chemical shift vector for cut down matrix

    all_peak_tables            structure- fields below
                               1. each_table - contains the SAND csvs
                               2. ft_peaks - contains all the SAND
                               underlying peaks fourier transformed

    peaks_to_add_struct        structure - fields below
                               1. ppms - ppm indices of the trace
                               2. fracs - fracs/rows of the trace from
                               matrix

    ends_ppm                   vector with two elements. ppm values where
                               the ends of the spectrum is. The data to
                               consider
  

    Outputs:

    X_added                   1D vector of the basis element

      
%}


mins = islocalmin(X,2);
end_min_idx = matchPPMs(ends_ppm,ppm);
mins(:,ppm < ends_ppm(1) | ppm > ends_ppm(2)) = 0;
mins(:,end_min_idx) = 1;
sum_over_fractions = zeros(1,length(all_peak_tables(1).ft_peaks(1).ft));
for j = 1:size(peaks_to_add_struct,2)
    for k = 1:length(peaks_to_add_struct(j).fracs)
        %Find the index of the peak in the ppm2
        max_idx = peaks_to_add_struct(j).ppms(k);
        %Locate the nearest minima surrounding the peak
        min1_idx = 0;
        min2_idx =0;
        for l = max_idx:-1:1
            if mins(peaks_to_add_struct(j).fracs(k),l) == 1
                min1_idx = l;
                break;
            end
        end
        for l = max_idx:length(ppm)
            if mins(peaks_to_add_struct(j).fracs(k),l) == 1
                min2_idx = l;
                break;
            end
        end
        % new structure field will be a structure. Each row corresponds to
        % each fraction the peak is in
        indexer2 = 1:size(all_peak_tables(peaks_to_add_struct(j).fracs(k)).each_table,1);
        multiple_sand = indexer2(all_peak_tables(peaks_to_add_struct(j).fracs(k)).each_table.freq_ppm > ppm(min1_idx) & all_peak_tables(peaks_to_add_struct(j).fracs(k)).each_table.freq_ppm < ppm(min2_idx));
        for l = 1:length(multiple_sand)
            sum_over_fractions = sum_over_fractions + all_peak_tables(peaks_to_add_struct(j).fracs(k)).ft_peaks(multiple_sand(l)).ft;
        end
    end
end
X_added = sum_over_fractions;

end
