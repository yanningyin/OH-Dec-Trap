clear; clc; close all; 
addpath('src')

rng('default') % fix the seed of the random nunber generator to be removed afterwards
params = Parameters();

example_no = 7; %{1: only dec; 2: only trap; 3: both}

switch example_no
    case 1
        % Example: only simulate deceleration
        tic;
        beam = Beam(params);
        beam.createParticles();
        dec = Deceleration(params, beam);
        dec.propagateParticles_euler();
        dec.plotTOF();
        dec.propagateParticles_verlet();
        dec.plotTOF();
        toc;
    case 2
        % Example: only simulate trapping
        tic;
        beam = Beam(params); % change params.BEAM_avg_velocity_beam and params.CALC_vel_synch_mol
        beam.createParticles();
        trap = Trapping(params, beam); % time adjustment needed
        trap.propagateParticles_verlet(10000);
        trap.plotTOF();
        toc;
    case 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
        % Example: simulate the whole process
        tic;
        beam = Beam(params);
        beam.createParticles();
        nsteps = 10000;
        dec = Deceleration(params, beam);
        dec.propagateParticles_verlet();
        trap = Trapping(params, dec);
        %trap.testnrsteps(); 
        % trajectory = zeros(nsteps+1,7);
        trajectory = [];
        k = 1;
        for i=1:length(trap.xyzVxyz_0)  % 对从减速器出来的粒子进行遍历，求它们的轨迹。
            trap.trajectory(trap.xyzVxyz_0(i,:),nsteps)  %进行轨迹计算
            if isempty(trap.traj_xyzVxyz)      % jump up untrapped molecule
                continue
            else
                plot3(trap.traj_xyzVxyz(:,1), trap.traj_xyzVxyz(:,2),trap.traj_xyzVxyz(:,3))

            end
            if length(trap.traj_time) == nsteps + 1
                trajectory(:,1,k) = trap.traj_time;
                trajectory(:,2:7,k) = trap.traj_xyzVxyz;
                k = k + 1;
            end
            
            hold on
        end
        fprintf('%d particles trapped.', k-1);
        % trap.propagateParticles_verlet(nsteps)
        % trap.plotTOF();

        % data = [trap.traj_time, trap.traj_xyzVxyz];      % one of OH molecules trajectory, it isn't all of trapped OH trajectory
        [a, b] = uiputfile('*.mat', 'Save trajectory MAT file'); % 让用户选择文件名
        if a ~= 0   % 检查是否按了“取消”
            filename = fullfile(b, a);  % 拼接完整路径
            save(filename, 'trajectory','-v7.3');  % 保存变量 trajectory 到 .mat 文件
        end
        toc;

    case 7
        % Replot the trajectories during loading and trapping
        path_to_traj_mat_file = './test.mat';
        load(path_to_traj_mat_file);
        for k=1:length(trajectory(1,1,:))
            plot3(trajectory(:,2,k), trajectory(:,3,k), trajectory(:,4,k));hold on;
        end

    case 4
        % Example: optimizing the trap loading parameters with three
        % parameters
        us = 1e-6;
        %current = 100;
        coil_onset_time = 3700*us:20*us:3900*us;
        coil_duration = 200*us:50*us:800*us;
        onset2 = 0*us:50*us:500*us;

        max_num_trapped = 0;
        best_current = 0;
        best_onset_time = 0;

        num_trapped = zeros(length(onset2), length(coil_onset_time), length(coil_duration));

        beam = Beam(params);
        beam.createParticles();
        dec = Deceleration(params, beam);
        dec.propagateParticles_verlet();
        for i = 1:length(onset2)
            for j = 1:length(coil_onset_time)
                for k = 1:length(coil_duration)
                    params = Parameters('TRAP_coil_onset_time', coil_onset_time(j),...
                        'TRAP_coil_duration', coil_duration(k)); % ...
                        %'TRAP_coil_current', current(i), ...
                    trap = Trapping_2coils(params, dec);
                    trap.propagateParticles_verlet(10000, onset2(i));
                    num_trapped(i,j,k) = trap.num_trapped_particles;
                end          
            end
        end

        max_num = max(num_trapped(:));
        [rows,cols,pages] = ind2sub(size(num_trapped),find(num_trapped == max_num)); %[rows, cols, pages] = find(num_trapped == max_num);
        fprintf('Highest number %d reached at the following params:\n', max_num);
        for k = 1:length(rows)
            fprintf('onset 2 %d, onset %d\n, duration %d\n', onset2(rows(k)), coil_onset_time(cols(k)), coil_duration(pages(k)));
        end
    case 5
        %trajectories
        tic;
        us = 1e-6;
        beam = Beam(params);
        beam.createParticles();
        dec = Deceleration(params, beam);
        dec.propagateParticles_verlet();
        full = [];
        verl = [];
        coil_onset_time = 3765*us;
        coil_duration = 635*us;
        res = [];
        params = Parameters('TRAP_coil_onset_time', coil_onset_time, 'TRAP_coil_duration', coil_duration);
        trap = Trapping(params, dec);
        trap.trajectory([-0.01,0,0,24.6,0,0]);
        plot(trap.traj_xyzVxyz(:,1), trap.traj_xyzVxyz(:,2))
        full(end+1) = 0;
        for i=1:1:length(trap.xyzVxyz_0(:,1))
            i;
            trap.trajectory(trap.xyzVxyz_0(i,:));
            if length(trap.traj_xyzVxyz) >= 10000
                full(end) = full(end)+1;     
            end
            hold on
            plot(trap.traj_xyzVxyz(:,1), trap.traj_xyzVxyz(:,2))
            hold off

            res = [trap.traj_time, trap.traj_xyzVxyz];
        end
        %xlim([-0.015 0.005])
        %ylim([-0.003 0.003])
            %hold on
            %xline(trap.onsets, color = 'black') -- too many lines if done for whole distribution
            %hold off
        save('data/traj/onset1_' + string( coil_onset_time ) + 'us.mat', 'res')
        savefig('data/traj/onset1_' + string( coil_onset_time ) + 'us.fig')
        trap.propagateParticles_verlet(10000) % for comparison
        verl(end+1) = trap.num_trapped_particles;
        %full/length(trap.xyzVxyz_0)
        verl/length(trap.xyzVxyz_0);
        toc;    
    
    case 6
        % Example: optimizing the trap loading parameters with two
        % parameters
        us = 1e-6;
        current = 100;
        coil_onset_time = 3600*us:5*us:3900*us;
        coil_duration = 400*us:5*us:700*us;

        max_num_trapped = 0;
        best_current = 0;
        best_onset_time = 0;

        num_trapped = zeros(length(coil_onset_time), length(coil_duration));

        beam = Beam(params);
        beam.createParticles();
        dec = Deceleration(params, beam);
        dec.propagateParticles_verlet();
        for j = 1:length(coil_onset_time)
            for k = 1:length(coil_duration)
                params = Parameters('TRAP_coil_onset_time', coil_onset_time(j),...
                    'TRAP_coil_duration', coil_duration(k),  ...
                    'TRAP_coil_current', current);%, ...
                trap = Trapping(params, dec);
                trap.propagateParticles_verlet(20000);
                num_trapped(j, k) = trap.num_trapped_particles;
            end
        end

        max_num = max(num_trapped(:));
        [rows,cols] = ind2sub(size(num_trapped),find(num_trapped == max_num)); %[rows, cols, pages] = find(num_trapped == max_num);
        fprintf('Highest number %d reached at the following params:\n', max_num);
        for k = 1:length(rows)
            fprintf('onset %d\n, duration %d\n', coil_onset_time(rows(k)), coil_duration(cols(k)));
        end


end






