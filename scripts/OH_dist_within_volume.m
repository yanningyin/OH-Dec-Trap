% clear;clc;
% 
% % trajectories: [Ntime x 7 x Ntraj]
% % channels: 1=t(s), 2=x(m), 3=y(m), 4=z(m), 5=vx(m/s), 6=vy(m/s), 7=vz(m/s)
% S = load('../data/traj/traj_selected_x2mm_after2000us_compact_20260303_183911.mat');
% trajectory = S.trajectory_kept; 
% 
% % --- constants ---
% amu = 1.66053906660e-27;     % kg
% m = 17 * amu;                % 17 amu in kg
% e  = 1.602176634e-19;        % J per eV (optional)
% 
% % --- thresholds in SI ---
% us = 1e-6;       % microseconds in seconds
% um = 1e-6;       % micrometers in meters
% t0  = trajectory(1,1,1)+2000e-6;   % 2000 us in seconds = 0.002 s
% xlim = 25*um;    % 25 um in meters
% ylim = 25*um;    % 25 um in meters
% zlim = 75*um;    % 75 um in meters
% 
% % --- extract channels ---
% t  = trajectory(:,1,:);
% x  = trajectory(:,2,:);
% y  = trajectory(:,3,:);
% z  = trajectory(:,4,:);
% vx = trajectory(:,5,:);
% vy = trajectory(:,6,:);
% vz = trajectory(:,7,:);
% 
% % --- kinetic energy per sample ---
% v2   = vx.^2 + vy.^2 + vz.^2;
% KE_eV = 0.5 * m .* v2./ e;           % eV
% mask = (t >= t0) & (abs(x) <= xlim) & (abs(y) <= ylim) & (abs(z) <= zlim);
% 
% avg_all_eV   = mean(KE_eV, 'all', 'omitnan');
% avg_filt_eV  = mean(KE_eV(mask), 'omitnan');
% 
% fprintf('Average KE (all):     %.6g eV\n', avg_all_eV);
% fprintf('Matched samples:      %d\n', nnz(mask));
% fprintf('Average KE (filtered): %.6g eV\n', avg_filt_eV);

filtered = 0;
outdir = '../data/traj/hist/';  % change if you want, e.g. outdir = '/swdata/yin/...';

if filtered
    % --- filtered energies (vector) ---
    KEf_eV = KE_eV(mask);
    KEf_eV = KEf_eV(isfinite(KEf_eV));   % remove NaN/Inf
    
    % --- choose bins (pick ONE of the following) ---
    
    % Option A: fixed number of bins
    numBins = 100;
    edges = linspace(min(KEf_eV), max(KEf_eV), numBins+1);
    
    % Option B: fixed bin width (uncomment and set)
    % binWidth = 1e-7; % eV
    % edges = (floor(min(KEf_eV)/binWidth)*binWidth) : binWidth : (ceil(max(KEf_eV)/binWidth)*binWidth);
    
    % --- histogram counts (no normalization) ---
    [counts, edges] = histcounts(KEf_eV, edges);
    
    % Bin centers for KEf_eV vs counts
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    
    % --- plot (counts histogram) ---
    % figure;
    % histogram(KEf_eV, edges);  % uses your edges, plots counts
    % xlabel('Kinetic energy (eV)');
    % ylabel('Counts');
    % title('Filtered KE distribution (counts)');
    
    % --- save outputs ---
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    
    matfile = fullfile(outdir, ['KE_hist_counts_' ts '.mat']);
    txtfile = fullfile(outdir, ['KE_hist_counts_' ts '.txt']);
    
    % Save .mat (raw energies + binning + results)
    save(matfile, 'KEf_eV', 'counts', 'edges', 'centers');
    
    % Save text: two columns [KE_center_eV, counts]
    data_out = [centers(:), counts(:)];
    fid = fopen(txtfile, 'w');
    fprintf(fid, '# KE_center_eV\tcounts\n');
    fprintf(fid, '%.12e\t%d\n', data_out.');   % transpose for column-wise fprintf
    fclose(fid);
    
    fprintf('Saved MAT: %s\n', matfile);
    fprintf('Saved TXT: %s\n', txtfile);

else
    mask = (t >= t0);
    % --- unfiltered energies (vector over all time & traj) ---
    KEall_eV =  KE_eV(mask);
    KEall_eV = KEall_eV(isfinite(KEall_eV));   % remove NaN/Inf
    
    % --- choose bins (pick ONE) ---
    
    % Option A: fixed number of bins
    numBins = 100;
    edges = linspace(min(KEall_eV), max(KEall_eV), numBins+1);
    
    % Option B: fixed bin width (uncomment and set)
    % binWidth = 1e-7; % eV
    % edges = (floor(min(KEall_eV)/binWidth)*binWidth) : binWidth : (ceil(max(KEall_eV)/binWidth)*binWidth);
    
    % --- histogram counts (no normalization) ---
    [counts, edges] = histcounts(KEall_eV, edges);
    centers = 0.5 * (edges(1:end-1) + edges(2:end));
    
    % --- plot (counts histogram) ---
    figure;
    histogram(KEall_eV, edges);
    xlabel('Kinetic energy (eV)');
    ylabel('Counts');
    title('Unfiltered KE distribution (counts)');
    
    % --- save outputs ---
    ts = datestr(now, 'yyyymmdd_HHMMSS');
    matfile = fullfile(outdir, ['KEall_hist_counts_' ts '.mat']);
    txtfile = fullfile(outdir, ['KEall_hist_counts_' ts '.txt']);
    
    % Save .mat
    save(matfile, 'KEall_eV', 'counts', 'edges', 'centers');
    
    % Save text: two columns [KE_center_eV, counts]
    data_out = [centers(:), counts(:)];
    fid = fopen(txtfile, 'w');
    fprintf(fid, '# KE_center_eV\tcounts\n');
    fprintf(fid, '%.12e\t%d\n', data_out.');   % transpose for fprintf
    fclose(fid);
    
    fprintf('Saved MAT: %s\n', matfile);
    fprintf('Saved TXT: %s\n', txtfile);

end