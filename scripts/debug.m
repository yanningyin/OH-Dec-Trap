%debug

n_row = 500;
n_traj = 17267;
tt  = trajectory(end-n_row:end,1,1:n_traj );
xx  = trajectory(end-n_row:end,2,1:n_traj );
yy  = trajectory(end-n_row:end,3,1:n_traj );
zz  = trajectory(end-n_row:end,4,1:n_traj );
vvxx = trajectory(end-n_row:end,5,1:n_traj );
vvyy = trajectory(end-n_row:end,6,1:n_traj );
vvzz = trajectory(end-n_row:end,7,1:n_traj );
% 
% speed = sqrt(vx.^2 + vy.^2 + vz.^2);
% 
% max(max(speed))


thr = 17;  % m/s

% speed: [Ntime x 1 x Ntraj]
speed = sqrt(vvxx.^2 + vvyy.^2 + vvzz.^2);

% max speed per trajectory over time (dim 1)
maxSpeed = squeeze(max(speed, [], 1));   % -> [Ntraj x 1]

% indices of trajectories to plot
idx = find(maxSpeed > thr);

fprintf('Plotting %d trajectories with max speed > %.1f m/s\n', numel(idx), thr);

figure; hold on; grid on;
xlabel('x (m)'); ylabel('y (m)'); zlabel('z (m)');
title(sprintf('Trajectories with max speed > %.1f m/s (last %d samples)', thr, n_row+1));

% plot each selected trajectory
for k = 1:8%numel(idx)
    j = idx(k);
    % plot3( squeeze(xx(:,1,j)), squeeze(yy(:,1,j)), squeeze(zz(:,1,j)) );
    plot3(squeeze(vvxx(:,1,j)), squeeze(vvyy(:,1,j)), squeeze(vvzz(:,1,j)));
end
axis equal;
view(3);