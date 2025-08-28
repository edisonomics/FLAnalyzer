function [trace_struct] = tracePeaksFL(X_matrix,ppm,vector_row_idx,vector_col_idx,ppm_range_trace_points,percentage_peak_min)
%{
    Chris Esselman 4.23.25
    Edited log
        
    tracePeaksFL
    
        Trace peaks as they elute over Fraction Library data. Input should
        be peak maxima in the data
 
    Inputs:
 
    X_matrix                   Data matrix

    ppm                        Chemical shift vector

    vector_row_idx             row indices for starting point for tracing peaks

    vector_col_idx             col indices for starting point for tracing peaks

    ppm_range_trace_points     how many ppm points to consider for next peak   

    percentage_peak_min        percentage for considering still part of trace
  

    Outputs:

    trace_struct               Contains 2 fields. 
                               1. The ppm idxs of the trace 
                               2. The rows the trace is in

      
%}
% Find all the peaks in the data
peak_index_each_fraction = struct;
for i = 1:size(X_matrix,1)
    [peak_index_each_fraction(i).Peaks, peak_index_each_fraction(i).Peaks_idx,~,peak_index_each_fraction(i).Prominence] = findpeaks(X_matrix(i,:));
    peak_index_each_fraction(i).ppm = ppm(peak_index_each_fraction(i).Peaks_idx);
end


nodes_down = struct;
nodes_up = struct;
trace_struct = struct;
big_counter = 1;
% Big loop for going through all of the nodes
for i = 1:length(vector_row_idx)
    % Create variables for holding ppm idx and rows
    ppm_holder = zeros(1,1);
    fraction_holder = zeros(1,1);
    % Variable for going in and out of while loop
    at_end = false;
    % Varibale for indexing down through the fractions
    location = vector_row_idx(i);
    counter = 1;
    while ~at_end
        % Go down one fraction
        location = location - 1;
        % Break if going down one fraction is below first fraction
        if location < 1
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % In future, just look at that column and model the shape ? %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the closest peak in the next fraction
        [~,closest_next_peak] = min(abs(peak_index_each_fraction(location).Peaks_idx - vector_col_idx(i)));
        % To start, just see if next peak is below a cutoff
        if (X_matrix(location,peak_index_each_fraction(location).Peaks_idx(closest_next_peak))) > ((X_matrix(vector_row_idx(i),vector_col_idx(i)))*((1/100)*percentage_peak_min)) && ...
                (abs(peak_index_each_fraction(location).Peaks_idx(closest_next_peak) - vector_col_idx(i)) < ppm_range_trace_points)
            ppm_holder(counter) = peak_index_each_fraction(location).Peaks_idx(closest_next_peak);
            fraction_holder(counter) = location;
            counter = counter + 1;
        else
            at_end = true;
        end
    end
    nodes_down(big_counter).ppms = ppm_holder;
    nodes_down(big_counter).fracs = fraction_holder;
    % Now going to go up
    % Create variables for holding ppm idx and rows
    ppm_holder = zeros(1,1);
    fraction_holder = zeros(1,1);
    % Variable for going in and out of while loop
    at_end = false;
    % Varibale for indexing down through the fractions
    location = vector_row_idx(i);
    counter = 1;
    while ~at_end
        % Go down one fraction
        location = location + 1;
        % Break if going down one fraction is below first fraction
        if location > size(X_matrix,1)
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % In future, just look at that column and model the shape ? %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Find the closest peak in the next fraction
        [~,closest_next_peak] = min(abs(peak_index_each_fraction(location).Peaks_idx - vector_col_idx(i)));
        % To start, just see if next peak is below a cutoff
        if (X_matrix(location,peak_index_each_fraction(location).Peaks_idx(closest_next_peak))) > ((X_matrix(vector_row_idx(i),vector_col_idx(i)))*((1/100)*percentage_peak_min)) && ...
                (abs(peak_index_each_fraction(location).Peaks_idx(closest_next_peak) - vector_col_idx(i)) < ppm_range_trace_points)
            ppm_holder(counter) = peak_index_each_fraction(location).Peaks_idx(closest_next_peak);
            fraction_holder(counter) = location;
            counter = counter + 1;
        else
            at_end = true;
        end
    end
    nodes_up(big_counter).ppms = ppm_holder;
    nodes_up(big_counter).fracs = fraction_holder;
    if (nodes_down(big_counter).ppms(1)) == 0 && (nodes_up(big_counter).ppms(1) ~= 0)
        trace_struct(big_counter).ppms = [vector_col_idx(i) nodes_up(big_counter).ppms];
        trace_struct(big_counter).fracs = [vector_row_idx(i) nodes_up(big_counter).fracs];
        big_counter = big_counter + 1;
    elseif nodes_up(big_counter).ppms(1) == 0 && (nodes_down(big_counter).ppms(1) ~= 0)
        trace_struct(big_counter).ppms = [flip(nodes_down(big_counter).ppms) vector_col_idx(i)];
        trace_struct(big_counter).fracs = [flip(nodes_down(big_counter).fracs) vector_row_idx(i)];
        big_counter = big_counter + 1;
    elseif nodes_up(big_counter).ppms(1) ~= 0 && (nodes_down(big_counter).ppms(1) ~= 0)
        trace_struct(big_counter).ppms = [flip(nodes_down(big_counter).ppms) vector_col_idx(i) nodes_up(big_counter).ppms];
        trace_struct(big_counter).fracs = [flip(nodes_down(big_counter).fracs) vector_row_idx(i) nodes_up(big_counter).fracs];
        big_counter = big_counter + 1;
    else 
        trace_struct(big_counter).ppms = vector_col_idx(i);
        trace_struct(big_counter).fracs = vector_row_idx(i);
        big_counter = big_counter + 1;
    end
end
% If a trace has less than 4 add to it the next closest ppm on up and down
for i = 1:size(trace_struct,2)
    if length(trace_struct(i).fracs) < 4
        if trace_struct(i).fracs(end) + 1 <= size(X_matrix,2)
            trace_struct(i).fracs(end + 1) = trace_struct(i).fracs(end) + 1;
            trace_struct(i).ppms(end + 1) = trace_struct(i).ppms(end);
        end
        if trace_struct(i).fracs(1) - 1 > 0
            trace_struct(i).fracs = [trace_struct(i).fracs(1)-1 trace_struct(i).fracs];
            trace_struct(i).ppms = [trace_struct(i).ppms(1) trace_struct(i).ppms];
        end
    end
end
end
