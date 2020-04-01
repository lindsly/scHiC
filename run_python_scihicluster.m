%% python test

% addpath(genpath('HiCluster'))

commandStr = 'py "C:\Users\lindsly\Google Drive\UMICH\Research\Single Cell Hi-C\Matlab\HiCluster\HiCluster\schicluster\single_cell.py" Mymatrix.txt';
[status, commandOut] = system(commandStr);
if status==0
    disp('Imputation complete');
end

load('imputed_matrix.mat')

figure
imagesc(log(Q))