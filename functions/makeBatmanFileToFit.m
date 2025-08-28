function makeBatmanFileToFit(X,ppm,where_to_save_files)
%{
    Chris Esselman 6.24.25
    Edited log
        
    makeBatmanFileToFit
    
        Make NMRdata.txt. This is the file which BATMAN will fit the data
        to
 
    Inputs:
 
    X                          Data matrix 

    ppm                        Chemical shift vector


    where_to_save_files        path to folder where NMRdata.txt will be
                               saved
            
    
%}

X_norm = zeros(size(X));
labels = strings(1,(size(X_norm,1)+1));
var_types = strings(1,(size(X_norm,1)+1));
labels(1) = "ppm";
var_types(1) = "string";
for i = 1:size(X,1)
    X_norm(i,:) = X(i,:)/trapz(X(i,:));
    labels(i + 1) = sprintf('x%d',i);
    var_types(i + 1) = "string";
end
Nmrdata = table('VariableNames',labels,'Size',[length(ppm) length(labels)],'VariableTypes',var_types);
Nmrdata.ppm = compose('%.6f',flip(ppm)');
for i = 1:size(X_norm,1)
    Nmrdata(:,i+1) = compose('%.6f',flip(X_norm(i,:))');
end
current_folder = pwd;
cd(where_to_save_files);
writetable(Nmrdata,'NMRdata.txt','Delimiter','tab')
cd(current_folder)
end