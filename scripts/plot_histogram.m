% plot histogram

% S = load("../data/traj/filtered_KE_hist_counts_20260303_154616.mat");
% S = load("../data/traj/KE_hist_counts_20260303_155208.mat");
% 
% 
% figure;
% bar(S.centers, S.counts, 1);  % 1 = full bin width
% xlabel('Kinetic energy (eV)');
% ylabel('Counts');
% title('Filtered KE histogram (from saved counts)');



% fname = '../data/traj/hist/KEall_hist_counts_20260303_191140.txt';  

fname = '../data/traj/hist/KE_hist_counts_20260303_190733.txt';  

A = readmatrix(fname, 'FileType','text', 'CommentStyle','#');

KE_center_eV = A(:,1);
counts       = A(:,2);

meanE_eV = sum(KE_center_eV .* counts) / sum(counts);
fprintf('Mean energy = %.12e eV\n', meanE_eV);

figure;
bar(KE_center_eV, counts, 1);   % 1 = full bin width look
xlabel('Kinetic energy (eV)');
ylabel('Counts');
title('KE histogram (counts)');