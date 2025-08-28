function basis_norm = makeBatmanOutputSyntheticFits(basis_final,ppm,where_to_save_files,mag_freq_mHz,norm_tf)
%{
    Chris Esselman 7.14.25
    Edited Log
    
    - 7.28.25 - Make it so the center is correct
    makeBatmanOutputSyntheticFits
        Make Batman output files from synthetic fits
 
    Inputs:
 
    
      
%}


% Get the current working directory and then cd to the directory where the
% batman output will be put
currentFolder = pwd;
cd(where_to_save_files);
mkdir PureSpectraTemplate

% Make a variable that will be used to write a file with all the variables
metabolites = strings(size(basis_final,2),1);

% Get the multi_data_user.csv file ready
multi_data = table;

% Calculate the ppm value to consider multiplet based on magnet frequency
multiplet_ppm = (10/(mag_freq_mHz*1000000))*10^6;

% Struct containing the multiplets
multiplet_struct = struct;

basis_norm = zeros(size(basis_final));
for i = 1:size(basis_final,1)
    norm_base_1 = basis_final(i,:)-min(basis_final(i,:));
    if norm_tf == true
        norm_base_2 = norm_base_1/trapz(norm_base_1);
    else
        [~,indx_max] = max(basis_final,[],"all");
        [row,~] = ind2sub(size(basis_final),indx_max);
        norm_trap = basis_final(row,:) - min(basis_final(row,:));
        norm_base_2 = norm_base_1/trapz(norm_trap);
    end
    basis_norm(i,:) = norm_base_2;
    % Make the puretemplatespectra
    pure_spectra = strings(length(ppm),2);
    pure_spectra(:,1) = compose('%.6f',ppm');
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
            if abs(ppm(locs(j))-ppm(locs(k))) <= multiplet_ppm && ((abs(pks(j)-pks(k)))/((pks(j)+pks(k))/2))*100 < 100
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
            if mins_loop(l) == 1 || l == matchPPMs(ppm(upfield_peak)-0.03,ppm)
                min1_idx = l;
                break;
            end
        end
        for l = downfield_peak:length(ppm)
            if mins_loop(l) == 1 || l == matchPPMs(ppm(downfield_peak)+0.03,ppm)
                min2_idx = l;
                break;
            end
        end
        multiplet_struct(i).centers(j) = mean(ppm(locs(cell2mat(C(j)))))-(mean(ppm(locs(cell2mat(C(j)))))-mean([ppm(min1_idx) ppm(min2_idx)]));
        multiplet_struct(i).ends{j} = [ppm(min1_idx) ppm(min2_idx)];
        multi_data_loop.J_constant(j) = sprintf('%.4f,%.4f',ppm(min1_idx),ppm(min2_idx));
        multi_data_loop.couple_code(j) = -2;
        multi_data_loop.Include_multiplet(j) = 1;
        multi_data_loop.overwrite_pos(j) = 'n';
        multi_data_loop.overwrite_truncation(j) = 'n';
        multi_data_loop.pos_in_ppm(j) = sprintf('%.4f',multiplet_struct(i).centers(j));
        multi_data_loop.Metabolite(j) = sprintf('metabolite%d',i);
        multi_data_loop.relative_intensity(j) = sprintf('%.4f',trapz(norm_base_2(ppm > ppm(min1_idx) & ppm < ppm(min2_idx))));
    end
    if i == 1
        multi_data = multi_data_loop;
    else
        multi_data = [multi_data;multi_data_loop];
    end
end
    
% Write the metaboliteList.csv
writematrix(metabolites,'metabolitesList.csv')
writetable(multi_data,'multi_data_user.csv')
cd(currentFolder)

end
