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
clear functions_folder
%% Allow for core edisonomics metabolomics toolbox functions to be used
%{
    **If on NMRbox, ignore following comments**

    If not on NMRbox, will need to do the same steps but for 
    metabolomics_toolbox.

    1. Download repository from github:
    github.com/edisonomics/metabolomics_toolbox
    2. Open zip and put folder on local machine
    3. Delete activate_metabolimics_toolbox code below
    4. Replace with:
        
    core_functions_folder = "path_to_core_functions"
    addpath(genpath(core_functions_folder))

    5. Replace "path_to_core_functions" with the path to the core
    functions folder
    ex/"~/metabolomics_toolbox-master/code"

%}
activate_metabolomics_toolbox
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
batman_folder = '/folder_path_to_save_batman_files/';

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
basis_final = makeBatmanOutputMultipletsFL(X_sand,ppm_ft,all_peak_tables,basis,start_end_ppm,batman_folder,mag_freq,'common_ppm',true);
%% Plot the basis set
annotation_fraction_plot(basis_final,ppm_ft)
%% Optional Step
%{
    If the last input of the below function is true, the BATMAN files will
    be the exact same as the ones created in the
    makeBatmanOutputMultipletsFL function above. The integral of each 
    basis set element equal to one.

    If the last input is false, the BATMAN files created will be so that
    each basis set element will be divided by the integral of the first
    basis set element
%}
basis_synth_norm = makeBatmanOutputSyntheticFits(basis_final,ppm_ft,batman_folder,mag_freq,false);
%% Create file of spectra that will be fit
%{
    Section below will create the file containing the spectra that will be
    fit using the metabolite basis set 

%}
current_folder = pwd;

% Change "path_to_folder_of_ft1_files_to_fit" to path where ft1 files that
% should be fit is located
folder_ft1_files_to_fit = "path_to_folder_of_ft1_files_to_fit";

%
cd(folder_ft1_files_to_fit)
spectra = load1DPipe();
[X_fit,ppm_fit] = Setup1D(spectra);
makeBatmanFileToFit(X_fit,ppm_fit,batman_folder);
cd(current_folder)
clear spectra folder_ft1_files_to_fit current_folder
