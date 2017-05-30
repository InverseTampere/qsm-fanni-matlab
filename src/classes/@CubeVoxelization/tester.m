clear VARIABLES;
addpath('../plotting/');

CylData = [0.5, 2.0, 0 0 1, 1 1 1, 0, 0, 1 1 1;...
           0.3, 1.0, 0.5 0 0.3, -1 0 0, 0, 0, 1 1 1];

for iCyl = 1:size(CylData,1)
   
    CylData(iCyl,6:8) = CylData(iCyl,6:8)./norm(CylData(iCyl,6:8));
    
end

BranchData = [1 0 1 1 0 0];
ModelData = {CylData, BranchData};

figure(1);

qsm = QSMCylindrical(ModelData);
BlockPartition = qsm.toVoxels(0.4,[-1,-1 0],[2 2 3]);

hold off;

figure(2);
cla;

VoxelMat=zeros(BlockPartition.size);
for i=1:size(VoxelMat,1)
    for j=1:size(VoxelMat,2)
        for k=1:size(VoxelMat,3)
            if ~isempty(BlockPartition.occupied_table{i,j,k})
                VoxelMat(i,j,k) = 1;
            end
        end
    end
end
[vol_handle]=VoxelPlotter(VoxelMat,0.01); 
%visual effects (I recommend using the FigureRotator function from MATLAB
%Centeral

hold on;
daspect([1,1,1]);
%-
hold off;