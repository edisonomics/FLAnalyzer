function [X_NoNoise,all_peak_tables_NoNoise] = rmvNoiseFL(X_sand,ppm_ft,all_peak_tables,intensity_cutoff)
%{
    Chris Esselman 10.8.25
    
    rmvNoiseFL
        Remove SAND peaks below a noise threshold that are not part of
        large peak

    Inputs:
    
    X_sand              Data matrix obtained from using readInSandFL fxn
                        Use this matrix to set the noise threshold
                        ex/ 140x32k double

    ppm_ft              ppm vector corresponding to columns of X_sand

    all_peak_tables     Structure obtained from using readInSandFL fxn

    intensity_cutoff    cutoff for SAND peaks to be excluded. double


    Outputs:

    all_peak_tables_NoNoise   Structure with noise SAND peaks removed

    X_NoNoise                 Data matrix with SAND peaks removed


%}

% First peak pick each fraction with and without intensity cutoff
mins = islocalmin(X_sand,2);
all_peak_tables_NoNoise = struct;
X_NoNoise = zeros(size(X_sand));
for i = 1:size(X_sand,1)
    [~,peaks_cutoff] = findpeaks(X_sand(i,:),"MinPeakHeight",intensity_cutoff);
    % Return an error if the cutoff is to high
    if isempty(peaks_cutoff)
        error('Error. Intensity Cutoff too high. No peaks!')
    end
    % Get an indexer for all the peaks
    indexer_main = zeros(1,length(all_peak_tables(i).ft_peaks));
    % Loop through each peak
    for j = 1:length(peaks_cutoff)
        %Find the index of the peak in the ppm
        max_idx = peaks_cutoff(j);
        %Locate the nearest minima surrounding the peak
        min1_idx = 0;
        min2_idx =0;
        for l = max_idx:-1:1
            if mins(i,l) == 1 || l == 1
                min1_idx = l;
                break;
            end
        end
        for l = max_idx:size(X_sand,2)
            if mins(i,l) == 1 || l == size(X_sand,2)
                min2_idx = l;
                break;
            end
        end
        % Add the indices to the main loop
        indexer_loop = 1:length(all_peak_tables(i).ft_peaks);
        indexer_loop = indexer_loop(all_peak_tables(i).each_table.freq_ppm > ppm_ft(min1_idx) & all_peak_tables(i).each_table.freq_ppm < ppm_ft(min2_idx));
        indexer_main(indexer_loop) = 1;
    end
    % Only add the SAND peaks that are above cutoff to the new struct
    all_peak_tables_NoNoise(i).each_table = all_peak_tables(i).each_table(logical(indexer_main),:);
    all_peak_tables_NoNoise(i).ft_peaks = all_peak_tables(i).ft_peaks(logical(indexer_main));
    % Sum the ft peaks to make the data matrix
    for j = 1:length(all_peak_tables_NoNoise(i).ft_peaks)
        X_NoNoise(i,:) = X_NoNoise(i,:) + all_peak_tables_NoNoise(i).ft_peaks(j).ft;
    end
    X_NoNoise(i,:) = X_NoNoise(i,:) - min(X_NoNoise(i,:));
end
end
