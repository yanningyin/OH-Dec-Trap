% trajectories: [Ntime x 7 x Ntraj]
% channels: 1=t(s), 2=x(m), 3=y(m), 4=z(m), 5=vx(m/s), 6=vy(m/s), 7=vz(m/s)

% --- constants ---
amu = 1.66053906660e-27;     % kg
m = 17 * amu;                % 17 amu in kg
e  = 1.602176634e-19;        % J per eV (optional)

% --- thresholds in SI ---
t0  = trajectory(1,1,1)+2000e-6;   % 2000 us in seconds = 0.002 s
xlim = 25e-6;    % 25 um in meters
ylim = 25e-6;    % 25 um in meters
zlim = 75e-6;    % 75 um in meters

% --- extract channels ---
t  = trajectory(:,1,:);
x  = trajectory(:,2,:);
y  = trajectory(:,3,:);
z  = trajectory(:,4,:);
vx = trajectory(:,5,:);
vy = trajectory(:,6,:);
vz = trajectory(:,7,:);

% --- kinetic energy per sample ---
v2   = vx.^2 + vy.^2 + vz.^2;
KE_J = 0.5 * m .* v2;        % Joules
KE_eV = KE_J ./ e;           % eV (optional)

% --- filter ---
mask = (t >= t0) & (abs(x) <= xlim) & (abs(y) <= ylim) & (abs(z) <= zlim);

% --- average ---
% avgKE_J  = mean(KE_J(mask),  'omitnan');
% avgKE_eV = mean(KE_eV(mask), 'omitnan');
% 
% fprintf('Matched samples: %d\n', nnz(mask));
% fprintf('Average KE: %.6g J (%.6g eV)\n', avgKE_J, avgKE_eV);



% Matched samples: 15322
% Average KE: 2.41932e-24 J (1.51002e-05 eV)


avgKE_J  = mean(KE_J,  'omitnan');
avgKE_eV = mean(KE_eV, 'omitnan');

fprintf('Matched samples: %d\n', nnz(mask));
fprintf('Average KE: %.6g J (%.6g eV)\n', avgKE_J, avgKE_eV);