##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = False
correct = False
domain_decomposition = True

# Motion
discrete_motion = True
looping = False

# Restart
restart_calculation = True
read_first_compositions = False

# Output
plotting = True
saving = True
write_restart = True
write_restart_discharged = False
write_restart_discarded = False

#%% Case
path_to_case = './Models' # Path to input folder parent
case_name = 'gFHR' # Name of the input folder
main_input_file = 'input' # Name of the main input file

#%% Serpent data (from Serpent input)

# Pebbles
pbed_file = 'fpb_pos'
r_pebbles = 2 # cm
fuel_mat_name = 'fuel'
fuel_frac = 1
fuel_pebble_universe_name = 'u_fuel_pebble'
graph_pebble_universe_name = 'u_graph_pebble'

# Others
pbed_universe_name = 'u_pb'
detector_names = ['flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']
#detectors = {
#    "flux_pebbles_thermal": {
#        "E": "1E-11 0.625E-6",
#        "normalized":"normalized"
#    },
#    "flux_pebbles_fast": {
#        "E": "1E-4 20",
#        "normalized":"normalized"
#    },
#    "power_pebbles": {
#        "dl": "u_pb",
#        "extra_cards":"dr -8 void"
#    }
#}
### Used to define detectors that are not cell wisely 
detector_extra=False

#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 280e6 # W
Nsteps = 1000
neutrons_per_cycle = 200000
decay_step = 0 # days

#%% Burnup cycle
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':1, 'r_fuel_kernel':0.02125, 'Ntrisos':9022, 'threshold_type':'passes', 
'threshold_dir':+1, 'threshold':8}}

#%% Motion
motion_direction = +1
recirc_threshold = 1 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = 1
    time_per_pass = 522/8 # days

# DEM
else:
    positions_folder = '/global/scratch/users/yvesrobert/HTR-10/pebble_positions/'
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
output_folder_name = 'gFHR_ctherm_1_compo' # Name of the output folder
verbosity_mode = 0
inventory_names = ['551370','551340','551360','541350','541380','541370','531310','Ce','Eu','Pd','Ag','471101','Cs','Sr','Kr','360880','360851','360830','360840','922350','922360','922380','942390','942400','942410','952410','380900']


#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 1
if write_restart_discharged:
    restart_discharged_write_every = 1
if write_restart_discarded:
    restart_discarded_write_every = 100


#%% Restart read
if restart_calculation:
    restart_step = 578
    restart_data = '/global/scratch/users/ludovicjantzen/HxF/equilibirum_8_passes/core_578.csv'
    restart_binary = '/global/scratch/users/ludovicjantzen/HxF/equilibirum_8_passes/input.wrk_578'
elif read_first_compositions:
    restart_binary = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/restart/first_compos.wrk'
    restart_data = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/data_first_compos.csv'

if saving:
    write_global = True
    write_incore = True
    write_reinserted = False
    write_discarded = True
    write_discharged = True

if plotting:
    plot_base = True
    plot_detectors = True
    plot_inventory = True
    plot_keff = True
    plot_cumulative = True
    plot_pass = True
    plot_geom = True

extra_fields = ['burnup', 'fima']

