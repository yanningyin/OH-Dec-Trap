classdef Trapping < handle 
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        params % class containing all relavant parameters

        ax_mag, ax_magcoil  % acceleration matrices; magnets/magnets and coil combined
        ay_mag, ay_magcoil
        az_mag, az_magcoil

        % ax_mag_extended, ax_magcoil_extended  % add some zeros outside the import acc fields for extropolation
        % ay_mag_extended, ay_magcoil_extended  % necessary?
        % az_mag_extended, az_magcoil_extended
        
        ax_mag_interpl,ay_mag_interpl,az_mag_interpl
        ax_magcoil_interpl,ay_magcoil_interpl,az_magcoil_interpl
        
        num_particles
        num_trapped_particles
        xyzVxyz_0                       % the vector of initial pos&vel, number_of_particles * 6
        xyzVxyz                         % the vector of computed pos&vel, number_of_particles * 6
        time                            % continue time recording
        t_start                         % to continue where decelerator left off
        has_the_simulation_been_run     % boolean, default is False in constructor, will be switched to True at the end of the simulation
        num_trajectories_saved          % number of trajectories to be saved
        arrival_time                    % TODO; WTF? this variable is defined in some plotting function?? Maybe better just a lcoal variable
        output                          % save xyzVxyz , t and flag at certain times and display here
        traj_xyzVxyz                    % trajectories x,y,z, Vx, Vy, Vz
        traj_time                       % ?????
        ind_particles                 % index of particles created at begining saving the indexes of the ones that make it to the end
        TOF_xyzVxyz                   % Time of Flight: variable to save xyzVxyz at the time step/steps of detection
        TOF_save
        TOF_profile
        nstep
    end

    methods
        %% constructor of the class
        function obj = Trapping(params, dec)
        % The argument "dec" can also accept "beam" object

            obj.params = params;
            obj.xyzVxyz_0 = dec.xyzVxyz;
            try
                obj.t_start = dec.M_time_vec(end-1)+ obj.params.FLY_incoupling_time;
                obj.xyzVxyz_0(:,1) = obj.xyzVxyz_0(:,1)-obj.params.PHYS_length_dec - obj.params.PHYS_exit_to_detection;
            catch %if "beam" is given instead of "dec"
                obj.t_start = 0;
                obj.xyzVxyz_0(:,1) = obj.xyzVxyz_0(:,1) + obj.params.PHYS_valve_to_dec- obj.params.PHYS_exit_to_detection;
            end
%            obj.t_start = 0;
            %delete(dec); clear dec;
            obj.num_trapped_particles = 0;

            %TODO: all the other methods in this class, including
            %   - loadAccelerationFields(obj) -- done
            %   - interpolateAccField(obj) -- done without extended matrix
            %   - propagateParticles_euler(obj) -- done (probably) 
            % definition of trapped? maybe Parameters.m needs trap size
            % check for loss during flight not yet implemented
            %   - plotTOF(obj) -- errors
            %   - etc.

            % load the acceleration fields
            tic;obj.loadAccelerationFields();toc;
            obj.InterpolateAccField()
        end


        %% Load the acceleration fields
        % updated to load directly the .mat files
        function loadAccelerationFields(obj)
            fprintf('Loading acceleration fields ...\t')

            load('data/acc/acc_combined/acc_ringmags') % load them in local workspace
            obj.ax_mag = accx*0.86; obj.ay_mag= accy*0.86; obj.az_mag= accz*0.86;  %因为OH的磁偶极矩没有1.41，而是1.21所以1.21/1.41=0.86
            clearvars accx accy accz % ugly: clear them from local workspace

            load('data/acc/acc_combined/acc_' + string( obj.params.TRAP_coil_current ) + 'A') % load them in local workspace
            obj.ax_magcoil = accx; obj.ay_magcoil= accy; obj.az_magcoil= accz;
            clearvars accx accy accz % ugly: clear them from local workspace
            
            fprintf('\tloaded\n')
        end
        
        %% Interpolate acceleration field
        function InterpolateAccField(obj)
            num_grids_x = 166;
            cx = 116;
            num_grids_y = 81;
            cy = 41;
            num_grids_z = 81;
            cz = 41;
            gridded_x = linspace(1, num_grids_x, num_grids_x);
            gridded_x = (gridded_x-cx)*0.001/10; %*mm/gu, gu=10, conversion to m
            gridded_y = linspace(1, num_grids_y, num_grids_y);
            gridded_y = (gridded_y-cy)*0.001/10;
            gridded_z =	linspace(1, num_grids_z, num_grids_z);
            gridded_z = (gridded_z-cz)*0.001/10;
% 
%             obj.ax_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
%             obj.ay_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);
%             obj.az_norm_extended = zeros(num_grids_x, num_grids_y, num_grids_z);          
%             obj.ax_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ax_norm, flip(-obj.ax_norm(2:110,:,:),1)), (obj.params.PHYS_number_of_electrodes - 1)/2, 1, 1), obj.ax_norm);
%             obj.ay_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.ay_norm, flip(obj.ay_norm(2:110,:,:),1)), (obj.params.PHYS_number_of_electrodes - 1)/2, 1, 1), obj.ay_norm);
%             obj.az_norm_extended(3:num_grids_x-2,3:43,3:43) = cat(1, repmat(cat(1,obj.az_norm, flip(obj.az_norm(2:110,:,:),1)), (obj.params.PHYS_number_of_electrodes - 1)/2, 1, 1), obj.az_norm);
            
            obj.ax_mag_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_mag,'linear','linear');
            obj.ay_mag_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_mag,'linear','linear');
            obj.az_mag_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_mag,'linear','linear');
            obj.ax_magcoil_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ax_magcoil,'linear','linear');
            obj.ay_magcoil_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.ay_magcoil,'linear','linear');
            obj.az_magcoil_interpl = griddedInterpolant({gridded_x, gridded_y, gridded_z}, obj.az_magcoil,'linear','linear');
             
        end


        function propagateParticles_verlet(obj, nsteps)
            obj.xyzVxyz = [];
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz_0(:,1:3);
            obj.xyzVxyz(:,4:6) = obj.xyzVxyz_0(:,4:6);
            dt = 1e-6;
            %nsteps = 2000; % to be optimised
            obj.nstep = nsteps;
            obj.time = linspace(obj.t_start, obj.t_start + dt*nsteps, nsteps);
%             ax_mag_interpl = obj.ax_mag_interpl;
%             ay_mag_interpl = obj.ay_mag_interpl;
%             az_mag_interpl = obj.az_mag_interpl;
%             dxyzVxyz = {@(y) [y(:, 4:6), ax_interpl(y(:, 1), y(:, 2), y(:, 3)),...
%                                              [0;obj.ay_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))],...
%                                              [0;obj.az_interpl(y(2:end, 1), y(2:end, 2), y(2:end, 3))]]                   
%                         @(y) [y(:, 4:6), -obj.ax_interpl(obj.params.PHYS_exit_to_detection - y(:, 1), y(:, 3), -y(:, 2)),...
%                                              [0;-obj.az_interpl(obj.params.PHYS_exit_to_detection - y(2:end, 1), y(2:end, 3), -y(2:end, 2))],...
%                                              [0;obj.ay_interpl(obj.params.PHYS_exit_to_detection - y(2:end, 1), y(2:end, 3), -y(2:end, 2))]]};
%                 
             for t = obj.time
                %if size(obj.xyzVxyz,1) == 0
                %    error("No molecules remain for trapping.")
                %end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.params.num_particles,i);
                if (t > obj.params.TRAP_coil_onset_time) && (t < obj.params.TRAP_coil_off_time)
                    ax = obj.ax_magcoil_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                    ay = obj.ay_magcoil_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                    az = obj.az_magcoil_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                else
                    ax = obj.ax_mag_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                    ay = obj.ay_mag_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                    az = obj.az_mag_interpl(obj.xyzVxyz(:,1), obj.xyzVxyz(:,2), obj.xyzVxyz(:,3));
                end
                x = obj.xyzVxyz(:,1)+obj.xyzVxyz(:,4)*dt+0.5*dt*dt*ax;
                y = obj.xyzVxyz(:,2)+obj.xyzVxyz(:,5)*dt+0.5*dt*dt*ay;
                z = obj.xyzVxyz(:,3)+obj.xyzVxyz(:,6)*dt+0.5*dt*dt*az;
                if (t+dt > obj.params.TRAP_coil_onset_time) && (t+dt < obj.params.TRAP_coil_off_time)
                    ax2 = obj.ax_magcoil_interpl(x,y,z);
                    ay2 = obj.ay_magcoil_interpl(x,y,z);
                    az2 = obj.az_magcoil_interpl(x,y,z);
                else
                    ax2 = obj.ax_mag_interpl(x,y,z);
                    ay2 = obj.ay_mag_interpl(x,y,z);
                    az2 = obj.az_mag_interpl(x,y,z);
                end
                vx = obj.xyzVxyz(:,4)+0.5*(ax+ax2)*dt;
                vy = obj.xyzVxyz(:,5)+0.5*(ay+ay2)*dt;
                vz = obj.xyzVxyz(:,6)+0.5*(az+az2)*dt;
                obj.xyzVxyz = [x, y, z, vx, vy, vz];       % 当前时刻所有分子的位置和速度，不包括之前的位置和速度。
                    %obj.xyzVxyz = obj.xyzVxyz + dxyzVxyz{mod(i, 2)+2*(~mod(i, 2))}(obj.xyzVxyz) * dt;
                obj.xyzVxyz = obj.xyzVxyz(abs(obj.xyzVxyz(:,2)) < 0.002 & abs(obj.xyzVxyz(:,3)) < 0.002 & obj.xyzVxyz(:,1) < 0.005 & obj.xyzVxyz(:,1)>-0.015, :);  %筛选分子
                    %obj.xyzVxyz = obj.xyzVxyz(abs(obj.xyzVxyz(:,2)) < 0.001 & abs(obj.xyzVxyz(:,3)) < 0.001 & (abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);
                    %probably control of range needed, like obj.xyzVxyz = obj.xyzVxyz((obj.xyzVxyz(:,1) < 0.004), :);
                coil_entrance_mask = ~(obj.xyzVxyz(:,1) > -6.5e-3 & obj.xyzVxyz(:,1) < -2e-3 & obj.xyzVxyz(:,2).^2 + obj.xyzVxyz(:,3).^2 > obj.params.TRAP_coil_inner_radis^2);
                obj.xyzVxyz = obj.xyzVxyz(coil_entrance_mask, :); % clipped by the inner diameter of the front coil
                ion_trap_mask = ~(abs(obj.xyzVxyz(:,1)) > obj.params.TRAP_ion_trap_horiz_surf2surf/2 & abs(obj.xyzVxyz(:,1)) < obj.params.TRAP_ion_trap_horiz_surf2surf/2 + obj.params.TRAP_ion_trap_thickness & abs(obj.xyzVxyz(:,2)) > obj.params.TRAP_ion_trap_vert_surf2surf/2);
                obj.xyzVxyz = obj.xyzVxyz(ion_trap_mask, :);
                
             end
            obj.xyzVxyz = obj.xyzVxyz((abs(obj.xyzVxyz(:,1)) < 0.002), :); %control if trapped  输出最后被trap的分子的位置和速度
            obj.num_trapped_particles = size(obj.xyzVxyz,1);        % trap的分子个数
            fprintf("%d out of %d particles trapped using a current of %d A and an onset time of %d us.\n",size(obj.xyzVxyz,1), obj.params.num_particles, obj.params.TRAP_coil_current, obj.params.TRAP_coil_onset_time);
        end

        function plotTOF(obj)
            
            xyzVxyz=obj.xyzVxyz;
            xyzVxyz(:,1:3)=xyzVxyz(:,1:3)+xyzVxyz(:,4:6)*(obj.time(end)-obj.time(end-1));                                                           
            obj.arrival_time = obj.time(end) - (xyzVxyz(:,1))./(xyzVxyz(:,4)) + obj.params.FLY_incoupling_time;
            %             figure;histogram(obj.arrival_time); - obj.params.PHYS_length_dec - obj.params.PHYS_exit_to_detection
            % min_time = 0; max_time = 6e-3;
            % num_bins = 6000;
            % binsize = (max_time - min_time)/num_bins;
            % binsize = 1e-6;
            % tof_profile = zeros(num_bins,1);
            
%             obj.params.FLY_detection_laser_diameter = 1.4e-3;
            indices_detected = abs(xyzVxyz(:,2)) < obj.params.FLY_detection_laser_diameter/2 & abs(xyzVxyz(:,3)) < 2e-3;%1.3e-3;
            arrival_time = obj.arrival_time(indices_detected);
            xyzVxyz = xyzVxyz(indices_detected, :);
            % figure;histogram(arrival_time);
            min_time = 0; 
            max_time = max(arrival_time)+1e-3;
            binsize = 1e-6;
            num_bins = ceil((max_time - min_time)/binsize);

            bin_begin_each_particle = floor((arrival_time - obj.params.FLY_detection_laser_diameter/2./xyzVxyz(:,4) - min_time)/binsize + 0.5);
            bin_end_each_particle = floor((arrival_time + obj.params.FLY_detection_laser_diameter/2./xyzVxyz(:,4) - min_time)/binsize + 0.5);
            swap_idx = bin_begin_each_particle > bin_end_each_particle;
            temp = bin_begin_each_particle(swap_idx);
            bin_begin_each_particle(swap_idx) = bin_end_each_particle(swap_idx);
            bin_end_each_particle(swap_idx) = temp;

            % max_time = max(bin_end_each_particle)+1000;
            % binsize = 1e-6;
            % num_bins = (max_time - min_time)/binsize;
            tof_profile = zeros(num_bins,1);
            t = min_time:binsize:max_time;
            for i = 1: length(bin_begin_each_particle)
                tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) = tof_profile(bin_begin_each_particle(i):bin_end_each_particle(i)) + 1;
            end
            
            tof_save_path = './result/tof_NM_' + ...
                    strrep( string(obj.params.FLY_voltage_on_electrodes), '.', 'p') + ...
                    'kV_' + 'dec_' + string(obj.params.CALC_vel_synch_mol) ...
                    + '_' + string(obj.params.FLY_target_velocity)+ '_' + string(obj.params.TRAP_coil_current) ...
                    + '_' + string(obj.params.TRAP_coil_onset_time) + '_' + string(obj.params.TRAP_coil_duration) ...
                    + '.dat';
            % tof = [t(2:end)', tof_profile];
            tof = [t', tof_profile];
%             save(tof_save_path, "tof");
            fileID = fopen(tof_save_path, 'w');
            fprintf(fileID, '%f\t%d\n', tof');
            obj.TOF_profile = tof;
            
            figure;
            % plot(t(1:end-1)*1e6, tof_profile);
            plot(t*1e6, tof_profile);
            xlabel('arrival time (us)')
            ylabel('signal (arb. u)')

        end

        function testnrsteps(obj)
            %nrsteps = [1000,2000,3000,3100,3200,3300,3400,3500,3600,3700,3800,3900,4000,4100,4200,4300,4400,4500,4600,4700,4800,4900,5000,5100,5200,5300,5400,5500,5600,5700,5800,5900,6000,7000,8000,9000,10000];
            nrsteps = [1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000, 12000, 13000, 14000, 15000, 16000, 17000, 18000, 19000, 20000];
            for n = nrsteps
                obj.propagateParticles_verlet(n)
            end
        end


        function trajectory(obj, posv, nsteps)     % get one molecule trajcetory whether it can be trapped or not
            obj.traj_time = [];
            obj.traj_xyzVxyz = [];
            dt = 1e-6;
            % maxnsteps = 10000;
            maxnsteps = nsteps;
            % maxnsteps = obj.nstep;
            interval = linspace(obj.t_start, obj.t_start + dt*maxnsteps, maxnsteps);
            traj = [interval(1), posv(1:6)];
            for t = interval
%                 if size(obj.xyzVxyz,1) == 0
%                     error("No molecules remain for trapping.")
%                 end
%                     fprintf("%d/%d\t%d\n",size(obj.xyzVxyz,1), obj.params.num_particles,i);
                if (t > obj.params.TRAP_coil_onset_time) && (t < obj.params.TRAP_coil_off_time)
                    ax = obj.ax_magcoil_interpl(traj(end,2), traj(end,3), traj(end,4));
                    ay = obj.ay_magcoil_interpl(traj(end,2), traj(end,3), traj(end,4));
                    az = obj.az_magcoil_interpl(traj(end,2), traj(end,3), traj(end,4));
                else
                    ax = obj.ax_mag_interpl(traj(end,2), traj(end,3), traj(end,4));
                    ay = obj.ay_mag_interpl(traj(end,2), traj(end,3), traj(end,4));
                    az = obj.az_mag_interpl(traj(end,2), traj(end,3), traj(end,4));
                end
                x = traj(end,2)+traj(end,5)*dt+0.5*dt*dt*ax;
                y = traj(end,3)+traj(end,6)*dt+0.5*dt*dt*ay;
                z = traj(end,4)+traj(end,7)*dt+0.5*dt*dt*az;
                if (t+dt > obj.params.TRAP_coil_onset_time) && (t+dt < obj.params.TRAP_coil_off_time)
                    ax2 = obj.ax_magcoil_interpl(x,y,z);
                    ay2 = obj.ay_magcoil_interpl(x,y,z);
                    az2 = obj.az_magcoil_interpl(x,y,z);
                else
                    ax2 = obj.ax_mag_interpl(x,y,z);
                    ay2 = obj.ay_mag_interpl(x,y,z);
                    az2 = obj.az_mag_interpl(x,y,z);
                end
                vx = traj(end,5)+0.5*(ax+ax2)*dt;
                vy = traj(end,6)+0.5*(ay+ay2)*dt;
                vz = traj(end,7)+0.5*(az+az2)*dt;
                
                %%% Molecules are considered lost if they are 1) out of the
                %%% acceleration field range (x < 5 mm, y/z < 2 mm) or 2)
                %%% radially out of coil's hole when passing through coil
                if (~(abs(traj(end,3)) < 0.002 && abs(traj(end,4)) < 0.002 && traj(end,2) < 0.005 && traj(end,2) > -0.015)) ...
                        || (traj(end,2)>-0.0065 && traj(end,2)<-0.002 && traj(end,3)^2+traj(end,4)^2>obj.params.TRAP_coil_inner_radis^2)
                    break
                end
                %%% Molecules considered lost if their postions overlap with
                %%% ion trap electrodes
                if abs(traj(end,2)) > obj.params.TRAP_ion_trap_horiz_surf2surf/2 ...
                        && abs(traj(end,2)) < obj.params.TRAP_ion_trap_horiz_surf2surf/2 + obj.params.TRAP_ion_trap_thickness ...
                        && abs(traj(end,3)) > obj.params.TRAP_ion_trap_vert_surf2surf/2
                    break
                end

                traj(end+1,:) = [t+dt, x, y, z, vx, vy, vz]; 
            end
            if abs(traj(end,2)) <= 0.002  %  trap control 
                obj.traj_time = traj(:,1);
                obj.traj_xyzVxyz = traj(:, 2:7);% [traj(:,2), traj(:,3), traj(:,4), traj(:,5), traj(:,6), traj(:,7)];
            else
                obj.traj_time = [];
                obj.traj_xyzVxyz = [];
            end
        end
        
    end
end