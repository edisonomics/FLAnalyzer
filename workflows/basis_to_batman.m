%% Script for outputting basis for BATMAN
%% Tell MATLAB which functions to use
%{
    **If the necessary functions were already downloaded from github
      and put into a folder, ignore following comments**
    
    1. Download repository from github: github.com/edisonomics/FLAnalyzer
    2. Open zip and put folder on local machine or NMRbox
    3. Below you will need to provide the path to the functions folder in
    the repository. ex/ "~/FLAnalyzer-main/functions"

%}

% Replace "path_to_FLAnalyzer_functions" with the path to the FLAnalyzer
% functions
functions_folder = "path_to_FLAnalyzer_functions";
addpath(genpath(functions_folder))
%% Load the necessary MATLAB variables for outputting to BATMAN
%{
    Load the .mat file with the variables needed for outputting to BATMAN

    Variables needed:
    1. X_sand - orgininal matrix of SANDed mFL
    2. ppm_ft - chemical shift vector corresponding to X_sand
    3. all_peak_tables - struct containing CSV SAND files
    ^^^^ these variables were acquired from using the readInSandFL function
    from the create_basis_set.m script
    
    4. basis - Struct containing the basis elements. This variable is the
    result of using the FLAnalyzer app

%}
% Change '/folder_path_to_save_batman_files/' with the path to the folder where you want the BATMAN input
% files to be saved
folder = '/folder_path_to_save_batman_files/';

% Change this to magnet frequency the data was acquired with. 
mag_freq = 799.713758637;

% Very Important!!
% Replace these values with the exact same values for the start ppm (Exclude DSS) 
% and End ppm overall parameters in the FLAnalyzer app
start_end_ppm = [0.1,9.5];

% Creating the files
%{
    There are two inputs to the below function to play with.

    1. Where the input 'common_ppm' is located, there are several input
    options

    a. 'common_ppm' - chooses the fraction of the basis set element with the most peaks in the element
    b. 'max' - chooses the fractions with the largest intensity
    c. 'least_peaks' - chooses the fraction with the least number of peaks in total

    Visually inspect the result of these options and choose the one that
    has the best results

   %%%%%%%%%%%%%%
    2. The last input, there are two input options

    a. true - include DSS in the basis set
    b. false - do not include DSS in the basis set

    If there is no DSS in the the spectrum that will be fit, change to
    false
    
%}
basis_final = makeBatmanOutputMultipletsFL(X_sand,ppm_ft,all_peak_tables,basis,start_end_ppm,folder,mag_freq,'common_ppm',true);
%% Plot the basis set
annotation_fraction_plot(basis_final,ppm_ft)
%% 
basis_synth_norm = makeBatmanOutputSyntheticFits(basis_final2,ppm_ft,folder,799.713758637,false);