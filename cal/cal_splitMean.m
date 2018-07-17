function [splitMean, ROI_data, Voxel_data, totalROIpCF] = cal_splitMean(betas_A,betas_B)

% function [condition_splitMean, voxel_MeanPeak, voxel_Mean, condition_splitMean_max, condition_splitMean_peak, VoxelIndex_A, VoxelIndex_B, VoxelMax_A, VoxelMax_B] = cal_splitMean(betas_A,betas_B)

% add error bars?
% variation in betas estiamted and indiviusdal voxel beta estimation

% output
% ROI centriod and spread, pRF properties - fit to bins
% voxel centriod and spread, pRF properties - fit to voxels


[VoxelMax_A, VoxelIndex_A] = max(double(betas_A));
[VoxelMax_B, VoxelIndex_B] = max(double(betas_B));
% compute split mean
for n = 1:size(betas_A,1)
    splitA(n,:) = sum(betas_A(:,VoxelIndex_B==n),2)/sum(VoxelIndex_B==n);
    splitB(n,:) = sum(betas_B(:,VoxelIndex_A==n),2)/sum(VoxelIndex_A==n);
    
    totalROIpCF(n) = sum(VoxelIndex_A==n) + sum(VoxelIndex_B==n);
end
splitMean = (double(splitA) + double(splitB)) / 2;

%% ROI properties
ROI_data = struct;
[ROI_data.splitMean_max, ROI_data.splitMean_peak] = max(splitMean);

% calculate weighted mean and spread
% [ROI_data.Cntrd, ROI_data.Sprd] = cal_WeightedMean(splitMean);
% pRF
% [ROI_data.pRF_params, ROI_data.rfModel, ROI_data.pRF_cntrd, ROI_data.pRF_spread, ROI_data.pRF_peak] = cal_pRF(splitMean, ROI_data.splitMean_peak, ROI_data.splitMean_max);

%% Voxel properties
Voxel_data = struct;
% voxel_Mean_Index = (VoxelIndex_A + VoxelIndex_B) / 2;
voxel_Mean = (double(betas_A) + double(betas_B)) / 2;
[Voxel_data.Mean_Max, Voxel_data.Mean_Peak] = max(voxel_Mean);

% calculate weighted mean and spread
% [Voxel_data.Cntrd, Voxel_data.Sprd] = cal_WeightedMean(voxel_Mean);
% pRF
% [Voxel_data.pRF_params, Voxel_data.rfModel, Voxel_data.pRF_cntrd, Voxel_data.pRF_spread, Voxel_data.pRF_peak] = cal_pRF(voxel_Mean, Voxel_data.Mean_Peak, Voxel_data.Mean_Max);

%%
% roi_av = mean([mean(betas_A,2),mean(betas_A,2)],2);

% perform same as below on condition_splitMean
% create weighted mean and fit gaussian functions

% %% ROI properties
% % calculate weighted mean and spread
% f = repmat((1:size(condition_splitMean,1))',1,size(condition_splitMean,2));
% condition_wA = (condition_splitMean-min(condition_splitMean))./(max(condition_splitMean)-min(condition_splitMean));
% condition_Cntrd = sum(condition_wA.*f)./sum(condition_wA);
% condition_Sprd = sqrt(sum(condition_wA.*(f-condition_Cntrd).^2)/sum(condition_wA));
% 
% % pRF
% % define gaussian function
% fun = @(x,xdata) x(3) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(4);
% lb = 0;
% ub = size(condition_splitMean,1).*1.25;
% xdata = 1:size(condition_splitMean,1);
% condition_pRF = zeros(size(condition_splitMean,2),4);
% 
% for i = 1:size(condition_splitMean,2)
%     % fit gaussian
%     x0 = [condition_splitMean_peak(i),1,1,condition_splitMean_max(i)];
%     condition_pRF(i,:) = lsqcurvefit(fun,x0,xdata,condition_splitMean(i,:),lb,ub);
% end
% 
% %% Voxel properties
% % calculate weighted mean and spread
% f = repmat((1:size(voxel_Mean,1))',1,size(voxel_Mean,2));
% voxel_wA = (voxel_Mean-min(voxel_Mean))./(max(voxel_Mean)-min(voxel_Mean));
% Voxel_Cntrd = sum(voxel_wA.*f)./sum(voxel_wA);
% Voxel_Sprd = sqrt(sum(voxel_wA.*(f-Voxel_Cntrd).^2)/sum(voxel_wA));
% 
% % pRF
% % define gaussian function
% fun = @(x,xdata) x(3) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(4);
% xdata = (1:size(voxel_Mean,1))';
% Voxel_pRF = zeros(size(voxel_Mean,2),4);
% for i = 1:size(voxel_Mean,2)
%     % fit gaussian
%     x0 = [voxel_Mean_Peak(i),1,1,voxel_Mean_Max(i)];
%     Voxel_pRF(i,:) = lsqcurvefit(fun,x0,xdata,voxel_Mean(:,i),lb,ub);
%     Voxel_rfModel{i} = fun(Voxel_pRF(i,:),xdata);
%     % Voxel_pRF_cntrd(i) =
% end

    function [centriod, spread] = cal_WeightedMean(data)
        xData = repmat((1:size(data,1))',1,size(data,2));
        wA = (data-min(data))./(max(data)-min(data));
        centriod = sum(wA.*xData)./sum(wA);
        spread = sqrt(sum(wA.*(xData-centriod).^2)/sum(wA));
                
    end

    function [pRF_params, rfModel, pRF_cntrd, pRF_spread, pRF_peak] = cal_pRF(data, Peak, Max)
        fun = @(x,xdata) x(3) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(4);
        lb = [0,0,-inf,-inf];
        ub = [size(data,1),size(data,1),+inf,+inf];
        xdata = (1:size(data,1))';
        pRF_params = zeros(size(data,2),4);
        opts = optimset('Display','off');
        for i = 1:size(data,2)
            % fit gaussian
            x0 = [Peak(i),1,1,Max(i)];
            pRF_params(i,:) = lsqcurvefit(fun,x0,xdata,data(:,i),lb,ub,opts);
            rfModel{i} = fun(pRF_params(i,:),xdata);
            [pRF_cntrd(i) pRF_spread(i)] = cal_WeightedMean(rfModel{i});
           [pRF_max(i), pRF_peak(i)] = max(rfModel{i});
        end
    end
% fit gaussian pRF

% output params of both

% voxel by voxel and condition by condition
% weight mean
%      wA = betasScanA(:,i)';
%      wA = (wA-min(wA))./(max(wA)-min(wA));
%      VoxelCntrdA(i) = sum(wA.*f)./sum(wA);
%      VoxelSprdA(i) = sqrt(sum(wA.*(f-VoxelCntrdA(i)).^2)/sum(wA));
%           wB = betasScanB(:,i)';
%      wB = (wB-min(wB))./(max(wB)-min(wB));
%      VoxelCntrdB(i) = sum(wB.*f)./(sum(wB));
%      VoxelSprdB(i) = sqrt(sum(wB.*(f-VoxelCntrdB(i)).^2)/sum(wB));
%      VoxelCntrd(i) = mean([VoxelCntrdA(i) VoxelCntrdB(i)]);
%      VoxelSprd(i) = mean([VoxelSprdA(i) VoxelSprdB(i)]);
%
% % pRF
% for i = 1:nBins-1
% fun = @(x,xdata) x(4) .* exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))) + x(5);
% % fun = @(x,xdata) x(4) .* (exp(-(((xdata-x(1)).^2)/(2*(x(2)^2))))).^x(3) + x(5);
% % fun = @(x,xdata) max(x(4) .* 1 + x(3).*log10((exp(-(((xdata-x(1)).^2)/(2*(x(2)^2)))))) + x(5),0);
% xdata = 1:length(ROIbetasSum);
% y = ROIbetasSum(i,:);
% [m index] = max(ROIbetasSum(i,:));
%
% x0 = [index,1,2,1,1];
% fpRF(i,:) = lsqcurvefit(fun,x0,xdata,y);
% rfModel{i} = fun(fpRF(i,:),xdata);
%
% subplot(subIndexSum(1),subIndexSum(2),i)
% plot(xdata,y)
% hold on
% plot(xdata,rfModel{i})
% % c = c+1;
% end
%
% figure
% xpRFplot = -32:32;
% for i = 1:nBins-1
% plot(xpRFplot,fun([0 fpRF(i,2:end)],xpRFplot));
% axis tight
% hold on
% end

end