function [X_removed,all_peak_tables_removed] = resumeAnalysisFL(X_sand,ppm_ft,all_peak_tables,basis,ends_ppm)

%{
    Chris Esselman 11.12.25
    
   resumeAnalysisFL
        Create X_removed and all_peak_tables_removed just from the basis
        structure. Saves memory when saving the basis elements as you go
        To continue from previous basis will use X_removed and 
        all_peak_tables_removed as input to FLAnalyzer

    Inputs:
    
    X_sand              Data matrix of orgininal sand data ex/ 140x32k double

    ppm_ft              ppm values corresponding to data matrix ex/ 1x32k

    all_peak_tables     Original peak tables from using readInSandFL

    basis               struct obtained from creating basis set a struct of
                        structs containing correlated traces

    ends_ppm            a 1X2 array containing the parameters used on
                        FLAnalyzer gui
                        ex/ ppm_exclude_dss_start and ppm_exclude_dss_end
                        First element is start and second element is end
                        

    Outputs:

    X_removed                  Data matrix with basis traces removed

    all_peak_tables_removed    Peak tables with traces from basis removed
    
        

%}
if ~isfield(basis,"removed_traces")
    % Remove the end regions
    X_removed = X_sand;
    ppm_removed = ppm_ft(:,ppm_ft < ends_ppm(2) & ppm_ft > ends_ppm(1));
    all_peak_tables_removed = all_peak_tables;
    for i = 1:size(basis,2)
        X_removed = X_removed(:,ppm_ft < ends_ppm(2) & ppm_ft > ends_ppm(1));
        [X_removed,all_peak_tables_removed] = removePeaksFL(X_removed,ppm_removed,all_peak_tables_removed,basis(i).correlated_traces,ends_ppm);

    end
else
    basis_full = basis;
    % First getting rid of the actual basis elements
    basis = rmfield(basis,"removed_traces");
    if ~isempty(fieldnames(basis))
        idx1 = zeros(1,size(basis,2));
        idx2 = 1:size(basis,2);
        for i = 1:size(basis,2)
            if isempty(basis(i).correlated_traces)
                idx1(i) = 1;
            end
        end
        basis(idx2(logical(idx1))) = [];
        X_removed = X_sand;
        ppm_removed = ppm_ft(:,ppm_ft < ends_ppm(2) & ppm_ft > ends_ppm(1));
        all_peak_tables_removed = all_peak_tables;
        for i = 1:size(basis,2)
            X_removed = X_removed(:,ppm_ft < ends_ppm(2) & ppm_ft > ends_ppm(1));
            [X_removed,all_peak_tables_removed] = removePeaksFL(X_removed,ppm_removed,all_peak_tables_removed,basis(i).correlated_traces,ends_ppm);
        end
    end
    % Now get rid of the removed traces
    idx1 = zeros(1,size(basis_full,2));
    idx2 = 1:size(basis_full,2);
    for i = 1:size(basis_full,2)
        if isempty(basis_full(i).removed_traces)
            idx1(i) = 1;
        end
    end
    basis_full(idx2(logical(idx1))) = [];
    for i = 1:size(basis_full,2)
        X_removed = X_removed(:,ppm_ft < ends_ppm(2) & ppm_ft > ends_ppm(1));
        [X_removed,all_peak_tables_removed] = removePeaksFL(X_removed,ppm_removed,all_peak_tables_removed,basis_full(i).removed_traces,ends_ppm);
    end
end
end


