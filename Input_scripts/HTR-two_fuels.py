##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = True
correct = False
domain_decomposition = False
use_decnfy_lib = True

# Motion
discrete_motion = False
looping = True

# Restart
restart_calculation = False
read_first_compositions = False

# Output
plotting = True
saving = True
write_restart = True
write_restart_discharged = True
write_restart_discarded = True

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'HTR-10_two_fuels' # Name of the input folder
main_input_file = 'input.inp' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'pebbles.inp'
r_pebbles = 3 # cm
# pebbles_dict = {'u_fuel1_pebble':{'mat_name':'fuel1', 'pebbles_frac':0.5, 'r_fuel_kernel':0.025, 'Ntrisos':8335, 'threshold_type':'551370', 'threshold_dir':+1, 'threshold':4e-5},
#                 'u_fuel2_pebble':{'mat_name':'fuel2', 'pebbles_frac':0.2, 'r_fuel_kernel':0.025, 'Ntrisos':8335},
#                 'u_graph_pebble':{'pebbles_frac':0.3}}

pebbles_dict = {'u_fuel1_pebble':{'mat_name':'fuel1', 'pebbles_frac':0.513, 'r_fuel_kernel':0.025, 'Ntrisos':8335, 'threshold_type':'burnup', 'threshold_dir':+1, 'threshold':72},
                'u_fuel2_pebble':{'mat_name':'fuel2', 'pebbles_frac':0.057, 'r_fuel_kernel':0.025, 'Ntrisos':8335, 'threshold_type':'passes', 'threshold_dir':+1, 'threshold':1},
                'u_graph_pebble':{'pebbles_frac':0.43}}

# Others
pbed_universe_name = 'u_pb'
detector_names = ['flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']
detectors = {'flux_pebbles_thermal':{'E':[1e-11, 1.86e-6], 'normalized':True},
             'flux_pebbles_fast':{'E':[0.18, 20], 'normalized':True},
             'power_pebbles':{'extra_cards':['dr', -8, 'void']}}

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 10e6 # W
Nsteps = 10000
neutrons_per_cycle = [1000]*20 + [5000]*600+(Nsteps-600-20)*[100000]
decay_step = 0 # days

#%% Burnup cycle
inventory_names = ['551370']

#%% Motion

# Discrete motion
if discrete_motion:
    motion_direction = -1
    Nrows_to_move = [20]*20+[10]*20+[5]*(Nsteps-20-20)
    time_per_pass = 522/8 # days

# DEM
else:
    recirc_threshold = 200 # cm, absolute value
    positions_folder = '/global/scratch/users/co_nuclear/pebble_positions_larger/'
    DEM_step_increment = 1
    circulation_rate = 125 # pebbles/days
    DEM_circulation_step = 1350 # pebbles
    positions_scale = 100
    positions_translation = [0,0, -610]
    if looping:
    # Looper
        DEM_start = 60
        DEM_end = 244
        looper_Nr = 5
        looper_Nz = 10
        looper_method = 'xyz'

#%% Outputing
output_folder_name = 'Nonpro_10' # Name of the output folder
verbosity_mode = 0

#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 10
if write_restart_discharged:
    restart_discharged_write_every = 1
if write_restart_discarded:
    restart_discarded_write_every = 1

#%% Restart read
if restart_calculation:
    restart_step = 650
    restart_data = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/Data/data_650.csv'
    restart_binary = '/global/scratch/users/yvesrobert/HTR-10_latest/Cases/HTR-10_P1_T1/wrk_Serpent/input.inp.wrk_650'
elif read_first_compositions:
    restart_binary = './Tools/Create_initial_compositions/restart/first_compos.wrk'
    restart_data = './Tools/Create_initial_compositions/data_first_compos.csv'


if saving:
    write_global = True
    write_incore = True
    write_reinserted = True
    write_discarded = True
    write_discharged = True

# if plotting:
#     core_plot = [['id', 'pass_residence_time', 'residence_time', 'passes', 'burnup'], detector_names + [f'integrated_{name}' for name in detector_names]]

if plotting:
    plot_base = False
    plot_detectors = False
    plot_inventory = False
    plot_keff = True
    plot_cumulative = False
    plot_pass = True
    plot_geom = False