% Define all the parameters here
function params = Parameters(varargin)
    params = struct();

    % Initial molecular beams
    params.num_particles = 10000;
    params.BEAM_avg_velocity_beam = 480; % Average velocity incoming package (m/s) 
    params.BEAM_radius_of_nozzle = 0.25e-3; % radius of the valve nozzle (m)
    params.BEAM_long_pos_spread = 11.5e-3; %  longitudinal position spread (m) - along x aixs or beam propagation
    params.BEAM_long_vel_spread = 0.2*params.BEAM_avg_velocity_beam ; % relative velocity spread along beam axis 0.112 0.12
    params.BEAM_trans_velocity_HWHM= 0.01*params.BEAM_avg_velocity_beam; %as longs as it is greater than 0.008*params.BEAM_avg_velocity_beam, it will reproduce the same distribution.
    
    % Deceleration
    params.FLY_voltage_on_electrodes = 13.5; % choose between 10, 12.5, 13. Will fail if fields do not exists!
    params.FLY_focusing_mode_bool = true; % flag for focusing mode operation
    params.CALC_vel_synch_mol = 450; % vx velocity synchronous molecule (m/s)
    params.FLY_target_velocity = 25; % target velocity (mhhh, overlap with the phase angle, but is needed to load appropriate velocity sequence)
    
    params.CALC_phase_degrees = 45.34;
    params.verbose = false;
    params.fortran_seq_bool = false;
    params.always_generate_M_seq = false;
    params.autom_find_final_vel = false;
    
    % Simion array -- DO NOT CHANGE
    params.SIMION_ni = 151; % ni number of grid points x (along beam axis
    params.SIMION_nj = 41; % number of grid points y,z (perpendicular to beam axis)
    params.SIMION_nk = 41; % number of grid points y,z (perpendicular to beam axis)
    params.SIMION_nbegin = 21; % Used part of array along beam axis, rest is for fitting. Start at nbegin...
    params.SIMION_nend = 131; % ... ends at nend
    params.SIMION_grid_units_p_meter = 20000.0; %gu #gridunits/meter (2*6666.6667)
    
    % Parameters of the PHYSICAL DECELERATOR
    params.PHYS_valve_to_skimmer = 253.0e-3; % LA Nozzle to skimmer (m)
    params.PHYS_skimmer_radius = 1.5e-3; %	r_s radius skimmer (m)
    params.PHYS_skimmer_to_dec = 66.6e-3; % LB Skimmer to Hexapole or Decelerator (m)
    params.PHYS_exit_to_detection = 10.2e-3;%11.52e-3; %L4 last stage decelerator to detection (m) % set to 10.2 mm 
    params.PHYS_number_of_electrodes = 123; % nt (Number of electrodes in dec.-1) size!! 
    params.PHYS_distance_stages = (params.SIMION_nend - params.SIMION_nbegin)/params.SIMION_grid_units_p_meter; % distance between two stages, 5.5 mm
    params.PHYS_valve_to_dec = params.PHYS_valve_to_skimmer + params.PHYS_skimmer_to_dec;
    params.PHYS_length_dec = params.PHYS_distance_stages * params.PHYS_number_of_electrodes; % total length of decelerator
    params.PHYS_valve_to_detection = params.PHYS_valve_to_skimmer + params.PHYS_skimmer_to_dec + params.PHYS_length_dec + params.PHYS_exit_to_detection;
    params.PHYS_seperation_pins = (params.SIMION_nj - 1)/params.SIMION_grid_units_p_meter;
    
    % Parameters used in theoretical calculations
    params.CALC_phase_distance = params.CALC_phase_degrees * params.PHYS_distance_stages / 180;
    
    params.CALC_OH_mass_amu = 17; %	mass Mass of molecule in AMU
    params.CALC_dipole_moment = 1.668; % mu dipole moment (Debye)
    params.CALC_lambda = 1.649e9; % Lambda-doublet splitting (Hz) (it was in GHz before)
    params.CALC_Beff = 0.58; % B effective value of MK/J(J+1)
    
    % Parameters used in fly and in the simulation
    params.FLY_incoupling_time = params.PHYS_valve_to_dec/params.CALC_vel_synch_mol; %710.2e-6; % valve - decelerator incoupling time (s)
    params.FLY_detection_laser_diameter = 1e-3;
    params.FLY_simulated_target_vel = []; % to be assigned while loading or generating the Matlab sequence.
    % Stores the rounded value of the last velocity vector. Ideally equals params.FLY_target_velocity
    
    
    % Trapping
    params.TRAP_coil_inner_radis = 1.5e-3;
    params.TRAP_ion_trap_horiz_surf2surf = 1.5e-3;
    params.TRAP_ion_trap_vert_surf2surf = 2.0e-3;
    params.TRAP_ion_trap_thickness = 0.5e-3;
    params.TRAP_coil_current = 100; % A
    params.TRAP_coil_onset_time = 3805e-6;%3.950e-3;% s
    params.TRAP_coil_duration = 400e-6;%440e-6; % s
    params.TRAP_coil_off_time = params.TRAP_coil_onset_time + params.TRAP_coil_duration;


    % Extra functionality: name-value pair arguments
    % This function can accept optional name-value pair arguments, useful for, e.g. optimization
    % Usage outside this script: params = Parameters('Trap_coil_current', 88, ...)
    parser = inputParser;
    addParameter(parser, 'TRAP_coil_current', params.TRAP_coil_current, @isnumeric);
    addParameter(parser, 'TRAP_coil_onset_time', params.TRAP_coil_onset_time, @isnumeric)
    addParameter(parser, 'TRAP_coil_duration', params.TRAP_coil_duration, @isnumeric)
    addParameter(parser, 'help', false, @islogical);

    try
        parse(parser, varargin{:});
    catch ME
        fprintf('Error: %s\n', ME.message);
        fprintf('Acceptable name-value pair arguments are:\n');
        fprintf('  ''TRAP_coil_current'' (numberic)\n');
        fprintf('  ''TRAP_coil_onset_time'' (numeric)\n');
        fprintf('  ''TRAP_coil_duration'' (numeric)\n');
        return;
    end

    params.TRAP_coil_current = parser.Results.TRAP_coil_current; % A
    params.TRAP_coil_onset_time = parser.Results.TRAP_coil_onset_time;% s
    params.TRAP_coil_duration = parser.Results.TRAP_coil_duration; % s
    params.TRAP_coil_off_time = params.TRAP_coil_onset_time + params.TRAP_coil_duration;
    
    fprintf('Input parameters:\n');
    disp(params);

    % Save all parameters to a .mat file
    % save('params.mat', '-struct', 'params');