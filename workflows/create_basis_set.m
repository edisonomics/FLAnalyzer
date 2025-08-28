%% Script for using main FLAnalyzer App
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

%% Read in ppm vector of ft1 files that were Sanded
%{
    This section will read in the chemical shift vector used for 
    the metabolite fraction library (mFL)

    **** Important ****

    Use the same chemical shift vector throughout all use of 
    FLAnalyzer. Using a different chemical shift vector will cause
    abnormalities. 

%}
current_folder = pwd;
% Change "path_to_folder_of_ft1_files" to the path of the folder where the
% ft1 files that were sanded is located
folder_ft1_files = "path_to_folder_of_ft1_files";
cd(folder_ft1_files)
spectra = load1DPipe();
[~,ppm_ft]=Setup1D(spectra);
cd(current_folder)
clear spectra current_folder folder_ft1_files
%% Read in the SANDed data

%{
   This section reads in the SANDed data. Before running this section, the
   csv file from each SANDed fraction must all be put into one folder. This folder
   will only contain the csv files.
%}
% Required parameters

% Path to csv files where csv files from the FL are located
path_csv = "path_to_SAND_csv_files";

% Magnet Frequency the data was collected with
mag_freq = 799.713758637;

% Sweep width the data was collected with
sweep_width = 12500;

%{ 
  1x2 vector. First value is start of sand region and second is end of sand
  region. 
  For example, if sand was run with:
       sand1D.com -in fraction1.ft1 -x1 9.5ppm -xn -0.5ppm -simAll
  
  the first value would be -0.5 and the second value would be 9.5
    ex/ sand_end_regions = [-0.5,9.5];
%}
sand_end_regions = [-0.5,9.625];

%{
  1x(even_number) vector. Pairs of numbers detailing regions for SAND peaks
  to not be considered.

  For example, if water region is only region wanting to be removed, the
  input would be something like this:
       vector_region_remove = [4.6,4.8];
  If wanted to add another region, the input would look something like
  this:
       vector_region_remove = [4.6,4.8,8.3,8.6];

%}
vector_region_remove = [4.66974,4.90586,8.3522,8.5845];

% Read in the data
[all_peak_tables,X_sand] = readInSandFL(path_csv,mag_freq,sweep_width,ppm_ft,sand_end_regions,vector_region_remove);

% Clear unnecessary variables
clear mag_freq norm_peak_region path_csv sand_end_regions sweep_width
clear vector_region_remove
%% Start the app
%{
    This next section gets the app running. The current folder in matlab
    must be where the app file (.mlapp) is located

    Follow the protocols.io on how to use the application

%}
% Change "path_to_folder_with_app" with the path to the folder were the
% create_basis.mlapp is located
folder_app = "path_to_folder_with_app";
cd(folder_app)
create_basis


