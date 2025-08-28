function [all_peak_tables,X_sand] = readInSandFL(path_csv,mag_freq,sweep_width,ppm_ft,sand_end_regions,vector_region_remove)
    %% Read in the SAND
    % find the max of a peak that is in ppm range of just noise. Any peak less
    % than this peak * 1.5 is deleted from the csvs
    files = dir([path_csv '/*.csv']);
    all_peak_tables = struct;
    for i = 1:size(files,1)
        all_peak_tables(i).each_table = readtable(fullfile(files(i).folder,files(i).name));
    end
    % Get rid of any SAND peaks that we do not want
    for i = 1:size(all_peak_tables,2)
        indexer = zeros(1,size(all_peak_tables(i).each_table,1));
        for k = 1:2:length(vector_region_remove)-1
            for j = 1:length(indexer)
                if all_peak_tables(i).each_table.freq_ppm(j) > vector_region_remove(k) && all_peak_tables(i).each_table.freq_ppm(j) < vector_region_remove(k+1)
                    indexer(j) = 1;
                end
            end
        end
        all_peak_tables(i).each_table(logical(indexer),:) = [];
    end
    operating_frequency_mhz = mag_freq;
    sw_hz = sweep_width;
    dw = 1/(sw_hz);
    % td_ = app.td;
    fid_time = length(ppm_ft)/2*dw;
    X_sand = zeros(size(files,1),length(ppm_ft));
    freq_hz = linspace(-sw_hz, sw_hz, length(ppm_ft));
    reconstructed_ppm = freq_hz / operating_frequency_mhz;
    for i = 1:size(all_peak_tables,2)
        t = linspace(0, fid_time, length(ppm_ft));
        % Pre-allocate array to hold all FFTs
        %reconstructedPeaks = zeros(size(all_peak_tables(i).each_table, 1), length(t));
        sumFID = 0;
        composite_spectrum = zeros(1,length(ppm_ft));
        for j = 1:size(all_peak_tables(i).each_table, 1)
            freq_hz = (all_peak_tables(i).each_table.freq_ppm(j) + 0) .* operating_frequency_mhz;
            decay_constant = 1 / (all_peak_tables(i).each_table.decay_hz(j));
            amplitude = all_peak_tables(i).each_table.amplitude(j);
    
            fid = amplitude * exp(-t / decay_constant) .* (cos(2 * pi * freq_hz * t) + 1i * sin(2 * pi * freq_hz * t));
            individual_peak = fftshift(fft(fid));
            all_peak_tables(i).ft_peaks(j).ft = real(individual_peak);
            all_peak_tables(i).ft_peaks(j).ft = interp1(reconstructed_ppm,all_peak_tables(i).ft_peaks(j).ft,ppm_ft,'linear');
            all_peak_tables(i).ft_peaks(j).ft = all_peak_tables(i).ft_peaks(j).ft - min(all_peak_tables(i).ft_peaks(j).ft);
            composite_spectrum = composite_spectrum + all_peak_tables(i).ft_peaks(j).ft;
            %sumFID = sumFID + fid; 
        end
        %composite_spectrum = fftshift(fft(sumFID));
        %X_sand(i,:) = real(composite_spectrum);
        %X_sand(i,:) = interp1(reconstructed_ppm,X_sand(i,:),ppm_ft,'linear');
        X_sand(i,:) = composite_spectrum;
        X_sand(i,:) = X_sand(i,:) - min(X_sand(i,:)); 
    end
    % Make sure the baselines are at the right levels
    %min_matrix = min(X_sand(:,ppm_ft > sand_end_regions(2)),[],"all");
    % for i = 1:size(X_sand,1)
    %     %min_row = min(X_sand(i,ppm_ft > sand_end_regions(2)));
    %     %diff = min_row - min_matrix;
    %     %X_sand(i,:) = X_sand(i,:) - diff;
    %     X_sand(i,:) = X_sand(i,:) - min(X_sand(i,:));
    % end
    


end