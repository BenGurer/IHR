function alignFunctional2HighResT2Star(xformFile,t2Star,filenames)

% We assume the s-form in this header transforms from (high res) t2Star to whole head space
t2Star = mlrImageReadNiftiHeader(t2Star);
t2starVoxelSizeXform = [t2Star.pixdim(2) 0 0 0; 0 t2Star.pixdim(3) 0 0;0 0 t2Star.pixdim(4) 0; 0 0 0 1];

% xform = transformation from the fMRI to the low resolution (same as fMRI) T2star 
% this is a transformation matrix output by FLIRT, so in world space (mm)
% we will have to take into account the voxel size
xform = textread(xformFile);
xform = xform(1:4,1:4);
if ieNotDefined('filenames')
    % Prompt user for filename(s)
    filenames = mlrGetPathStrDialog(pwd,'Choose one or more nifti files',{'*.img;*.nii','NIFTI Files'},'on');
end
if ~iscell(filenames)
    filenames = {filenames};
end

for p = 1:length(filenames)
    hdr = mlrImageReadNiftiHeader(filenames{p});
    fMRIvoxelSizeXform = [hdr.pixdim(2) 0 0 0; 0 hdr.pixdim(3) 0 0;0 0 hdr.pixdim(4) 0; 0 0 0 1];
    
    % compose transformation matrix
    % - Convert fMRI voxel indices to mm (fMRIvoxelSizeXform), this is because xform matrices output by FLIRT are in world coordinates
    % - align to low-res T2 star (xform), 
    % - convert from mm to high resolution T2* voxel size (inv(t2starVoxelSizeXform)) 
    % - apply transformation from high res T2* to whole head (T2StarHDR.sform44)
    hdr = cbiSetNiftiSform(hdr,t2Star.sform44 * inv(t2starVoxelSizeXform) * xform * fMRIvoxelSizeXform);
%     hdr = cbiSetNiftiSform(hdr,T2StarHDR.sform44 * (t2starVoxelSizeXform \ xform) * fMRIvoxelSizeXform);
    
    hdr.sform_code = 1;
    hdr = mlrImageWriteNiftiHeader(hdr,filenames{p});
end
