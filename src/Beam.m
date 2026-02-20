classdef Beam < handle 
    % This class is for beam generation.

    properties
        params % class containing all relavant parameters
        xyzVxyz_0                       % the vector of initial pos&vel, number_of_particles * 6
        xyzVxyz
    end

    methods
        %% constructor of the class
        function obj = Beam(params)
            obj.params = params;
            obj.xyzVxyz_0 = [];
            obj.xyzVxyz = [];
            % obj.createParticles();
        end
        
        
          %% create particles
          function createParticles(obj)
          % This function can create particles with the same phase-space
          % distribution as the Fortran code, as longs as
          % obj.params.BEAM_trans_velocity_HWHM > 0.008*obj.params.BEAM_avg_velocity_beam.

            fprintf('Creating molecular beam...')
            rng('default') % fix the seed of the random nunber generator to be removed afterwards
            obj.xyzVxyz_0 = [];
            while size(obj.xyzVxyz_0, 1) < obj.params.num_particles
                trial_num_particles = obj.params.num_particles * 3;
                xVx_0 = randn(trial_num_particles, 2)/sqrt(8*log(2)).* [obj.params.BEAM_long_pos_spread, obj.params.BEAM_long_vel_spread] + [-obj.params.PHYS_valve_to_dec, obj.params.BEAM_avg_velocity_beam]; 
                %conversion factor to std from full width half max, below adjust normal dist. to represent particles %set v_x to avergae and set all particles to valve pos. (decc x=0)
                rAVrB = rand(trial_num_particles, 4).*[obj.params.BEAM_radius_of_nozzle.^2, 2*pi, (obj.params.BEAM_trans_velocity_HWHM).^2, 2*pi];
                % produce radom transverse positions and azimuthal angles, and then project to y and z axis
                yzVyz_0 = [sqrt(rAVrB(:,1)).* cos(rAVrB(:,2)), sqrt(rAVrB(:,1)).* sin(rAVrB(:,2)), sqrt(rAVrB(:,3)).* cos(rAVrB(:,4)), sqrt(rAVrB(:,3)).* sin(rAVrB(:,4))];
                xyzVxyz_0 = [xVx_0(:,1),yzVyz_0(:,1), yzVyz_0(:,2), xVx_0(:,2),yzVyz_0(:,3), yzVyz_0(:,4)];
                % 存储从nozzle喷射出的粒子的位置和速度
                xyzVxyz_0 = xyzVxyz_0((xyzVxyz_0(:,2) + xyzVxyz_0(:,5).*(-obj.params.PHYS_skimmer_to_dec - xyzVxyz_0(:,1))./xyzVxyz_0(:,4)).^2 + (xyzVxyz_0(:,3) + xyzVxyz_0(:,6).*(-obj.params.PHYS_skimmer_to_dec - xyzVxyz_0(:,1))./xyzVxyz_0(:,4)).^2 <= obj.params.PHYS_skimmer_radius^2, :);
                % 筛选可以通过skimmer的粒子
                obj.xyzVxyz_0 = [obj.xyzVxyz_0; xyzVxyz_0];
                % 把可以通过skimmer的粒子的位置和速度信息保存，添加到队尾。
%                 size(obj.xyzVxyz_0, 1)
            end
            obj.xyzVxyz_0 = obj.xyzVxyz_0(1:obj.params.num_particles, :);
            obj.xyzVxyz_0(1,:) = [-obj.params.PHYS_valve_to_dec, 0, 0, obj.params.CALC_vel_synch_mol, 0, 0]; % The first row is for a synchronous molecule
            % 手动调整第一个粒子为同步粒子
            obj.xyzVxyz = obj.xyzVxyz_0;

            fprintf('\tcreated\n')
          end

        
%% Propagate particles, intergate using time intervals in M_time_vector, propagating all the molecules together
        function propagateToEntranceOfDec(obj)
            obj.createParticles();
            obj.xyzVxyz = obj.xyzVxyz_0;
            obj.xyzVxyz(:,1:3) = obj.xyzVxyz(:,1:3) + obj.xyzVxyz(:,4:6) * obj.params.FLY_incoupling_time; % ekin at start of deacc.
            obj.xyzVxyz = obj.xyzVxyz(abs(obj.xyzVxyz(:,2)) < 0.001 & abs(obj.xyzVxyz(:,3)) < 0.001 & (abs(obj.xyzVxyz(:,1) - obj.xyzVxyz(1,1)) < 2* 5.5e-3), :);   %限定初始在减速器前的位置
        end
    end
end
