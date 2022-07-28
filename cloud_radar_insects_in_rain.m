%% Code to examine insect flight during rain
% This code accompanies the manuscript "Using cloud radar to investigate
% the effect of rainfall on migratory insect flight" by Charlotte E.
% Wainwright, Sabrina  N. Volponi, Phillip M. Stepanian, Don R. Reynolds, and David
% H. Richter

%% Step 1: Download cloud radar data. The data used in this paper can be
% accessed at https://adc.arm.gov/discovery/. Data are free to download but registration with the ARM program is required.
% Data are available in 1-hr increments. Here we use data
% from 30 July 2015, including the time period 15-22 UTC. Here we use
% co-polar data, which is listed as SITEkazrspeccmaskgecopolC1.a0.YYYYDDMM.TIME.cdf
% where SITE is a site ID (here SGP). Data are stored in a packed
% format to save space. This code reconstructs and processes the data.

close all
clear all
clc

%% List radar files
files = dir('sgpkazrspeccmaskgecopol*20150730.*cdf');

%% Set threshold radial velocity values to delineate insects/rain [in m/s]
% Here set a 1 m/s, 0.5 and 2 m/s are also demonstrated in the paper
vrThresholdMin = -1.0;
vrThresholdMax = 1.0;

% Process each hour-long file one by one
for fileIdx=1:length(files)
    ncid = netcdf.open(files(fileIdx,1).name);
    basetime = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'base_time'));       % Base time [s since 1 January 1970]
    time = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'time_offset'));         % Time offset [s]
    range = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'range'));              % Height bins [m]
    vr = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'velocity_bins'));         % Radial velocity bins [m/s]
    mask = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'locator_mask'));        % Data mask
    spectraRaw = netcdf.getVar(ncid,netcdf.inqVarID(ncid,'spectra'));          % Spectra data
    specScale = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'spectra'),'scale_factor'); % Spectra scale factor
    specOff = netcdf.getAtt(ncid,netcdf.inqVarID(ncid,'spectra'),'add_offset');     % Spectra offset
    netcdf.close(ncid)
    
    % Reconstruct the full date and time from components in file
    for t=1:length(time)
        datetime(t) = datenum([1970 1 1 0 0 double(basetime)+time(t)]);
    end
    nCt = find(mask>-9999);
    
    % Reconstruct spectra accounting for offset and scale factor
    spectraRaw = double(spectraRaw).*specScale+specOff;
    
    % Rearrange the spectra to match the time-height-radial velocity format
    spectra = nan(length(time),length(range),length(vr));
    for count = 1:(length(unique(mask))-2)
        [ind1 ind2] = find(mask==count);
        spectra(ind2,ind1,1:256) = spectraRaw(:,count);
    end
    
    % NaNs are set to -9999 during processing
    spectra(isnan(spectra)>0)=-9999;
    %     save(['kazrspec' files(fileIdx,1).name(30:37) '_' files(fileIdx,1).name(39:40) 'Z_XC.mat'],'vr','range','time','datetime','spectra')
    
    % Remove heights above level we care about (here, we use 3 km as the limit)
    % from each variable
    heightMax = find(range<3000,1,'last');
    range = range(1:heightMax);
    spectra = spectra(:,1:heightMax,:);
    
    % Apply crude velocity dealiasing (velocity specified for this case)
    spectra = cat(3,squeeze(spectra(:,:,157:end)),squeeze(spectra(:,:,1:156)));
    vrTest = [vr(157:end)-11.9268; vr(1:156)];

    % Calculate noise level at each time and height using Hildebrand-Sekhon
    % algorithm
    for tCount=1:size(spectra,1)
        for htCount=1:size(spectra,2)
            [SnR2] = radarNoise(10.^(spectra(tCount,htCount,:)./10),length(vr));
            SNRcrit(tCount,htCount) = 10*log10(SnR2);
        end
    end
    SNRcrit(find(isinf(SNRcrit)>0)) = NaN; % Remove invalid points
    SNRmat = repmat(SNRcrit,[1,1,256]);    % Replicate noise matrix to match spectra size
    
    % Copy the spectra and apply the blob labeling technique to the copy
    spectraBlob = spectra;
    
    % Remove data falling below the noise threshold
    spectraBlob(find(spectraBlob<SNRmat)) = NaN;
    spectraBlob(find(spectraBlob<-9990)) = NaN;  % Turn missing data back to NaN
    
    % Turn this into a binary 0/1 matrix (Fig. 2d in paper)
    spectraBlob(find(isnan(spectraBlob)>0)) = 0;
    spectraBlob(find(spectraBlob<0)) = 1;
    
    % Setup empty matrices ready for the blob labeling
    labels = nan(size(spectraBlob));
    
    % Setup empty matrices to save blob mean, min, max radial velocities
    blobVRmean = nan(size(spectraBlob));
    blobVRmin = nan(size(spectraBlob));
    blobVRmax = nan(size(spectraBlob));
    
    
    for tCount=1:size(spectraBlob,1)
        % Apply connected component labeling, number each blob sequentially
        % (Fig. 2e in paper)
        blobIndex = double(labelmatrix(bwconncomp(squeeze(spectraBlob(tCount,:,:)))));
        labels(tCount,:,:) = blobIndex;
        
        % Setup temporary variables used during processing
        blobVRmeantmp = nan(size(spectraBlob,2),size(spectraBlob,3));
        blobVRmintmp = nan(size(spectraBlob,2),size(spectraBlob,3));
        blobVRmaxtmp = nan(size(spectraBlob,2),size(spectraBlob,3));
        
        % Find mean, min, max radial velocity for each blob
        for blob = 0:nanmax(nanmax(blobIndex))
            [x,y] = find(blobIndex==blob);
            ind = sub2ind(size(blobIndex),x,y);
            blobVRmeantmp(ind) = nanmean(vr(y));
            blobVRmintmp(ind) = nanmin(vr(y));
            blobVRmaxtmp(ind) = nanmax(vr(y));
        end
        blobVRmean(tCount,:,:) = blobVRmeantmp;
        blobVRmin(tCount,:,:) = blobVRmintmp;
        blobVRmax(tCount,:,:) = blobVRmaxtmp;
        clear blobIndex blobVRmeantmp blobVRmintmp blobVRmaxtmp
    end
    
    % Setup matrix to store only the part of the signal classed as insects
    insectSpectra = spectra;
    % Remove any parts of the signal outside the velocity threshold
    insectSpectra(find(blobVRmin<vrThresholdMin)) = NaN;
    insectSpectra(find(blobVRmax>vrThresholdMax)) = NaN;
    
    % Calculate total reflectivity (including insects and rain)
    totalReflectivity = pow2db(nansum(db2pow(spectra),3))';
    
    % Calculate total reflectivity (including insects and rain)
    insectReflectivity = pow2db(nansum(db2pow(insectSpectra),3))';
    
    % Range correct to give reflectivity in dBZ
    htRefDB = 20*log10(range);
    totalReflectivity = totalReflectivity + repmat(htRefDB,[1,size(totalReflectivity,2)]);
    insectReflectivity = insectReflectivity + repmat(htRefDB,[1,size(totalReflectivity,2)]);
    
    %% the data can now be saved in any format you wish
    % key variables to save are datetime (time in native MATLAB format),
    % range, totalReflectivity, and insectReflectivity. Example to save MAT
    % file shown on line below. File will be named radar_YYYYMMDD_HH.mat
    % save(['radar_' datestr(datetime(1),'yyyymmdd') '_' datestr(datetime(1),'HH') '.mat'])
end

function [SNR] = radarNoise(signal,iters)
% Calculate signal-to-noise ratio based on Hildebrand and Sekhon algorithm
signal = sort(signal, 'descend');
SNR = signal(end); % Set signal-to-noise threshold to lowest value in spectrum.
for i = 1:iters
    noiseVar = nanmean(signal(i:end).^2) - (nanmean(signal(i:end)))^2;
    testThresh = (nanmean(signal(i:end)))^2/noiseVar;
    if (testThresh > 1)
        SNR = signal(i);
        break
    end
end
end

%% Sample lines of code to produce figures from the paper:
% Individual spectra-height plots as in Fig. 4b-e can be plotted using:
% pcolor(vr, range, spectra)