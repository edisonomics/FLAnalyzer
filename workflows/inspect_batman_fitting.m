%% Script for inspecting BATMAN fitting
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

%% Make plot of summed fits
% Replace "path_to_Batman_output" with the path to the output of the BATMAN
% run. ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"
output_folder = "path_to_Batman_output";

% Can replace second argument with larger number if fit more than one
% spectrum in a run. For example, if wanting to inspect the ninth
% spectra out of ten spectra that were fit in a single BATMAN run, second
% argument would be 9
plotOverallAfterBatman(output_folder,1)

%% Make plot of underlying basis sets
% Replace "path_to_Batman_output" with the path to the output of the BATMAN
% run. ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"
output_folder = "path_to_Batman_output";

%{
If fit more than one spectrum in a single run, can inspect a different
spectrum than the first one by changing a few values below. For example,
if wanting to inspect the third spectrum fit, replace all instances of
easy_mult(1) with easy_mult(3)
%}
[easy_mult,ppm] = plotBasisElementsAfterBatman(output_folder,true,5);
figure
hold on 
plot(ppm,easy_mult(1).original,'LineWidth',2,'Color',"k")
plot(ppm,easy_mult(1).metabFit,'LineWidth',2,'Color',"r")
plot(ppm,easy_mult(1).metabolites)
hold off
set(gca,'XDir','rev')
%% Calculate the percent quantified
% Replace "path_to_Batman_output" with the path to the output of the BATMAN
% run. ex/ "~/runBATMAN/BatmanOutput/11_Jun_14_06_43"
output_folder = "path_to_Batman_output";

percent_quantified = calcPercentQuantified(output_folder);