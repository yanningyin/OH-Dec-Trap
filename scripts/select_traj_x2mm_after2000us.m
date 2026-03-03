clear;clc;

% --- settings ---
xlim = 2e-3;        % 2 mm in meters
dt0  = 2000e-6;     % 2000 us in seconds

% trajectory is [Ntime x 7 x Ntraj]
t = trajectory(:,1,:);   % [Ntime x 1 x Ntraj]
x = trajectory(:,2,:);   % [Ntime x 1 x Ntraj]

% --- per-trajectory threshold time: t_start + 2000 us ---
t0 = squeeze(t(1,1,:)) + dt0;      % [Ntraj x 1]
t0 = reshape(t0, 1, 1, []);        % [1 x 1 x Ntraj] for broadcasting

after = (t >= t0);                 % [Ntime x 1 x Ntraj]
viol  = after & isfinite(x) & (abs(x) > xlim);

keep = squeeze(~any(viol, 1));     % [Ntraj x 1] logical
idx_keep = find(keep);

fprintf('Kept trajectories: %d / %d\n', numel(idx_keep), numel(keep));
trajectory_kept = trajectory(:,:,keep);  % [Ntime x 7 x Nkeep]

outFile = sprintf('traj_selected_x2mm_after2000us_compact_%s.mat', datestr(now,'yyyymmdd_HHMMSS'));
save(outFile, 'trajectory_kept', 'keep', 'idx_keep', 'xlim', 'dt0', '-v7.3');

fprintf('Saved: %s\n', outFile);