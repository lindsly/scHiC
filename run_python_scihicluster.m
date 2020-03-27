%% python test

commandStr = 'py "C:\Users\lindsly\Google Drive\UMICH\Research\Single Cell Hi-C\Matlab\HiCluster\HiCluster\schicluster\single_cell.py" Mymatrix.txt';
[status, commandOut] = system(commandStr);
if status==0
 fprintf('squared result is %d\n',str2num(commandOut));
end