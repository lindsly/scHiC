if (exist('chrStart','var')) == 0
    load('chrStart.mat')
end
if (exist('chrStart','var')) == 0
    load('H_indiv_cells_1_114.mat')
end


%% Pre-process
% Trimming
A = H_indiv_cells_1_114(:,:,1);
trimLocs = zeros(1,size(A,1));
trimLoc = [];

for j = 1:size(A,1)
    if A(j,j) == 0
        trimLocs(1,j) = 1;
        trimLoc = [trimLoc;j];
    end 
end 

A(trimLoc,:) = [];
A(:,trimLoc) = [];
A(end,:) = [];
A(:,end) = [];

% Save into txt file
fid = fopen('Mymatrix.txt','wt');

for ii = 1:size(A,1)
    fprintf(fid,'%g\t',A(ii,:));
    fprintf(fid,'\n');
end


% Get new gene position   
chrEnd = [];

for i = 2:length(chrStart)
    chrEnd = [chrEnd;chrStart(i)-1];
end 

chrEnd = [chrEnd;size(H_indiv_cells_1_114,1)];
genePos = [chrStart chrEnd];
genePosNew = ceil(genePos);
trimLocIdx = sort(find(trimLocs),'descend');

for iTrimLocs = 1:length(trimLocIdx)
    genePosNew(genePosNew >= trimLocIdx(iTrimLocs)) =...
        genePosNew(genePosNew >= trimLocIdx(iTrimLocs))-1;
end

genePosNew(end,:) = [];





%% Get single chromosome 
load('test')
chrStartnew = genePosNew(:,1);
Chr = {};
Q = mydata;

for i = 1:length(chrStartnew)-1
    Chr{end+1} = Q(chrStartnew(i):chrStartnew(i+1)-1,chrStartnew(i):chrStartnew(i+1)-1);
end 
Chr{end+1} = Q(chrStartnew(length(chrStartnew)):size(Q,1),chrStartnew(length(chrStartnew)):size(Q,1));


%%
chr1 = Chr{1};
imagesc(log(chr1))

D = diag(sum(chr1,1));
L = D^(-1/2)*(D - chr1)*D^(-1/2);
[V S] = eig(L);











