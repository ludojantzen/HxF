##### Input parameters #####

#%% Calculation options
# Transport/depletion
transport = True
resolve_first = True
correct = False
domain_decomposition = True

# Motion
discrete_motion = True
looping = False

# Restart
restart_calculation = False
read_first_compositions = False

# Output
plotting = True
saving = True
write_restart = True
write_restart_discharged = True

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

# Others
pbed_universe_name = 'u_pb'
detector_names = ['flux_pebbles_thermal', 'flux_pebbles_fast', 'power_pebbles']
detectors = {
    "flux_pebbles_thermal": {
        "E": "1E-11 0.625E-6",
        "normalized":"normalized"
    },
    "flux_pebbles_fast": {
        "E": "1E-4 20",
        "normalized":"normalized"
    },
    "power_pebbles": {
        "dl": "u_pb",
        "extra_cards":"dr -8 void"
    }
}
#%% Depletion steps
power_normalization_field = 'power'
power_normalization_value = 280e6 # W
Nsteps = 1000
# neutrons_per_cycle = [2000]*68*2+[500000]*(Nsteps-68*2)
neutrons_per_cycle = [20000]*Nsteps
decay_step = 0 # days

#%% Burnup cycle
pebbles_dict = {'u_fuel_pebble':{'mat_name':'fuel', 'pebbles_frac':1, 'r_fuel_kernel':0.02125, 'Ntrisos':9022, 'threshold_type':'burnup', 'threshold_dir':+1, 'threshold':168}}

#%% Motion
motion_direction = +1
recirc_threshold = 1 # cm, absolute value

# Discrete motion
if discrete_motion:
    Nrows_to_move = [6]*Nsteps
    time_per_pass = 522/8 # days

 

#%% Outputing
output_folder_name = 'gFHR_test' # Name of the output folder
verbosity_mode = 1
inventory_names = ['551370','551340','551360','541350','541380','541370','531310','Ce','Eu','Pd','Ag','471101','Cs','Sr','Kr','360880','360851','360830','360840','922350','922360','922380','942390','942400','942410','952410','380900']


#%% Domain decomposition
if domain_decomposition:
    allowed_decomposition_types = 'rs'
    max_domains = [6, 8]

#%% Restart write
if write_restart:
    restart_write_every = 10
if write_restart_discharged:
    restart_discharged_write_every = 1
    

#%% Restart read
if restart_calculation:
    restart_step = 240
    restart_data = '/global/scratch/users/co_nuclear/gFHR_equilibrium/fine_eq/core_240.csv'
    restart_binary = '/global/scratch/users/co_nuclear/gFHR_equilibrium/fine_eq/input.wrk_240'
elif read_first_compositions:
    restart_binary = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/restart/first_compos.wrk'
    restart_data = '/global/scratch/users/yvesrobert/HxF/Tools/Create_initial_compositions/data_first_compos.csv'

if saving:
    write_global = True
    write_incore = True
    write_reinserted = True
    write_discarded = True
    write_discharged = True

extra_fields = ['burnup', 'fima']
