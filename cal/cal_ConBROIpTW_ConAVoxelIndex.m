function ConBROIpTW = cal_ConBROIpTW_ConAVoxelIndex(ConA_VoxelIndex,ConB_runA,ConB_runB)
for i = 1:size(ConB_runA,1)
    A(i,:) = sum(ConB_runA(:,ConA_VoxelIndex==i),2)/sum(ConA_VoxelIndex==i);
    B(i,:) = sum(ConB_runB(:,ConA_VoxelIndex==i),2)/sum(ConA_VoxelIndex==i);
end
ConBROIpTW = (A + B) / 2;
end