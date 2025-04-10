import os
import shutil
import numpy as np
from glob import glob
import re
import matplotlib.pyplot as plt
import subprocess
import sys
import pandas as pd
from cerberus.solvers import CodeInput, Solver
import cerberus as cb
import time
import datetime

np.random.seed(12345)

def print_with_timestamp(message):
    timestamp = datetime.datetime.now().strftime("[%m/%d/%Y %H:%M:%S]")
    print(timestamp, message, flush=True)

def process_input(input_files, pbed_universe_name):
	# Clean input files
	if isinstance(input_files, str):
		input_files = rec_Serpent_file_search(input_files, verbose=False)
	print_with_timestamp("Cleaning input files...")
	patterns = (r"^set[ \t]depout[ \t]", r"^dep[ \t]", r"^set[ \t]inventory[ \t]", r"^det(?=.*dl[ \t]+{0}[ \t])".format(re.escape(pbed_universe_name)), r"^set[ \t]rfr[ \t]", r"^set[ \t]rfw[ \t]", r"^set[ \t]acelib[ \t]", r"^set[ \t]declib[ \t]", r"^set[ \t]nfylib[ \t]", r"^set[ \t]dd[ \t]", r"^set[ \t]relfactor[ \t]", r"^set[ \t]pcc[ \t]")
	for file_path in input_files:
		modified_lines = []
		removed_lines = []

		with open(file_path, "r") as f:
			lines = f.readlines()

		for line in lines:
			if any(re.match(pattern, line) for pattern in patterns):
				removed_lines.append(line)
			else:
				modified_lines.append(line)

		if removed_lines:
			print_with_timestamp(f"Lines removed from {file_path}:")
			for removed_line in removed_lines:
				print_with_timestamp(f'\t{removed_line.strip()}')
			with open(file_path, "w") as f: # Write the modified lines back to the file
				f.writelines(modified_lines)

	# Add to input important lines
	with open(input_files[0], 'a') as f:
		f.write('\nset depout 2 2\n') # removes .dep and dep.m files
		f.write('\ndep butot 1\n') # make depletion calculation

def start_Serpent(serpent_executable, ncores, input_files, nnodes=1, verbosity=1):
	cb.LOG.set_verbosity(verbosity)
	cb.LOG.set_compiled_table(True)
	if isinstance(input_files, str):
		input_files = rec_Serpent_file_search(input_files)
	
	if nnodes > 1:
		# It seems that on Savio, SLURM needs another process to be spawned to make Cerberus understand it needs to book the right number of nodes
		print_with_timestamp("Dummy Serpent")
		from mpi4py import MPI
		child_comm = MPI.COMM_SELF.Spawn(serpent_executable, args=["tuto1"], maxprocs=nnodes - 1) # replace by dummy command that does not crash

	print_with_timestamp("Initializing Serpent. Waiting...")
	serpent = Solver("Serpent", serpent_executable, f"-omp {ncores} -port".split())
	serpent.input = CodeInput(input_files, main_input_idx=0)
	if nnodes==1:
		serpent.initialize()
	else:
		serpent.initialize(use_MPI=True, n_MPI_procs_to_spawn=nnodes)

	return serpent, input_files

def plot_Serpent(tra_plot, plots_folder_path=None, output_suffix=None, output_folder=None, nplots=None, delay=60, serpent_instance=None):
	Serpent_set_values(tra_plot, 1, serpent_instance=serpent_instance)
	if plots_folder_path:
		plot_files = []
		while len(plot_files) < nplots:
			plot_files = [file for file in glob(f'{plots_folder_path}/*_geom*.png') if time.time()-os.path.getmtime(file)<=delay]
			time.sleep(0.1)
		plot_files.sort(key=lambda x: os.path.getmtime(x))
		plot_files = plot_files[-nplots:]
		if output_folder:
			new_plot_files = [output_folder for file in plot_files]
		else:
			new_plot_files = [plots_folder_path for file in plot_files]
		for i in range(len(plot_files)):
			if not output_suffix:
				output_suffix = ''
			new_plot_files[i] += f'/{os.path.basename(plot_files[i]).split(".png")[0]}{output_suffix}.png'
			shutil.move(plot_files[i], new_plot_files[i])

def count_plots(serpent_input_files):
	cnt = 0
	for file in serpent_input_files:
		try:
			with open(file) as f:
				lines = [line.strip()[:4] for line in f.readlines()]
			for line in lines:
				if line=='plot':
					cnt += 1
		except:
			pass
	return cnt
#%%

def rec_Serpent_file_search(main_input, verbose=True, level=0, log_delim='\t'):
	if verbose:
		print_with_timestamp(f'{log_delim*level}Reading: {main_input}')
	files_to_read = [main_input]
	if main_input.split('.')[-1] in ['stl','ifc']:
		return files_to_read
	elif '.wrk' in main_input:
		return files_to_read
	path = os.path.dirname(main_input)
	with open(main_input) as f:
		lines = f.readlines()
	for i_line in range(len(lines)):
		line = lines[i_line].strip().split('%')[0]
		if len(line) == 0:
			continue
		fields = line.split()
		cmd = fields[0]
		if cmd == 'include':
			field_spot = 1
		elif cmd == 'pbed':
			field_spot = 3
		elif cmd == 'file':
			field_spot = 2
		elif cmd == 'set' and len(fields)>1:
			if fields[1] == 'dd' and fields[2] =='5':
				field_spot = 3
			# elif fields[1] == 'rfr':
			# 	if fields[2] == 'idx':
			# 		field_spot = 4
			# 	else:
			# 		field_spot = 3
			else:
				continue
		elif cmd == 'ifc':
			field_spot = 1
		else:
			continue

		if isinstance(field_spot, int):
			field_spot = [field_spot]

		if len(fields) > max(field_spot):
			for i in field_spot:
				file_name = fields[i].replace('"', '')
				files_to_read += rec_Serpent_file_search(os.path.normpath(os.path.join(path, file_name)), verbose=verbose, level=level+1, log_delim=log_delim)
		else:
				spots = list(fields)
				while len(spots) <= max(field_spot):
					i_line += 1
					line = lines[i_line].strip().split('%')[0]
					while len(line) == 0:
						i_line += 1
						line = lines[i_line].strip().split('%')[0]
						if i_line > len(lines):
							raise Exception('Lost here')
					fields = line.split()
					for j in fields:
						spots.append(j)

				for j in field_spot:
					file_name = spots[j].replace('"', '')
					files_to_read += rec_Serpent_file_search(os.path.normpath(os.path.join(path, file_name)), verbose=verbose, level=level+1, log_delim=log_delim)
	if len(files_to_read) > 1:
		if verbose:
			print_with_timestamp(f'{log_delim*(level+1)}Additional files ({len(files_to_read)-1}): {files_to_read[1:]}')

	return files_to_read

def analyze_sbatch_file(path_sbatch_file):
	with open(path_sbatch_file) as f:
		lines = f.readlines()

	ncores = ''
	for line in lines:
		fields = line.strip().split()
		if len(fields) == 0:
			continue

		if fields[0] == 'mpirun':
			for i in range(1, len(fields)):
				if fields[i] == '-np':
					nnodes = fields[i+1]
				elif fields[i] == '-omp':
					ncores = fields[i+1]
		if fields[0] == '#SBATCH' and '--ntasks=' in fields[1]:
			ntasks = int(fields[1].split('--ntasks=')[-1])
		if fields[0] == '#SBATCH' and '--cpus-per-task=' in fields[1]:
			cpus_per_task = int(fields[1].split('--cpus-per-task=')[-1])

	# Replace symbols
	if nnodes == '$SLURM_JOB_NUM_NODES':
		nnodes = ntasks
	else:
		nnodes = int(nnodes)

	if ncores == '' or ncores == '$SLURM_CPUS_PER_TASK':
		ncores = cpus_per_task
	else:
		ncores = int(ncores)

	return nnodes, ncores

def use_disperser(geom_type, args_geom, N_particles, r_particles, particles_uni, output_file_path, cores, grow_and_shake=True, shake_factor=0.05, growth_factor=0.1, serpent_executable='sss2', print_output=True):
	if isinstance(args_geom, (float,int)):
		args_geom = [args_geom]
	if isinstance(N_particles, (float, int)):
		N_particles = [N_particles]
	N_types = len(N_particles)
	if isinstance(r_particles, (float,int)):
		r_particles = [r_particles for _ in range(N_types)]
	if isinstance(particles_uni, (float,int,str)):
		particles_uni = [particles_uni for _ in range(N_types)]

	geom_type = str(geom_type)
	if geom_type in ['1', 'sphere', 'sph']:
		N_args = 1
		geom_type = '1'
		geom_name = 'sphere'
	elif geom_type in ['2', 'cylinder', 'cyl']:
		N_args = 3
		geom_type = '2'
		geom_name = 'cylinder'
	elif geom_type in ['3', 'cube']:
		N_args = 1
		geom_type = '3'
		geom_name = 'cube'
	elif geom_type in ['4', 'annular cylinder']:
		N_args = 4
		geom_type = '4'
		geom_name = 'annular cylinder'
	elif geom_type in ['5', 'cuboid']:
		N_args = 3
		geom_type = '5'
		geom_name = 'cuboid'
	elif geom_type in ['6', 'parallelepiped', 'par', 'paral']:
		N_args = 6
		geom_type = '6'
		geom_name = 'parallelepiped'
	else:
		raise Exception(f'"{geom_name}" is a wrong geometry trype, please choose a number between 1 and 6 or one string in the following list: "sphere", "cylinder", "cube", "annular cylinder", "cuboid", "parallelepiped"')

	if len(args_geom)!= N_args:
		raise Exception(f'Type {geom_name} takes {N_args} arguments, and {len(args_geom)} are given')

	printf_fields = []
	printf_fields.append(geom_type)
	for arg in args_geom:
		printf_fields.append(f'{arg}')
	for i in range(N_types):
		for param in [N_particles[i], r_particles[i], particles_uni[i]]:
			printf_fields.append(f'{param}')
		if i+1 != N_types:
			printf_fields.append('y')
		else:
			printf_fields.append('n')
	printf_fields.append(output_file_path)
	if grow_and_shake:
		printf_fields.append('y')
		printf_fields.append(f'{shake_factor}')
		printf_fields.append(f'{growth_factor}')
	else:
		printf_fields('n')
	string_printf = '\\n'.join(printf_fields)
	if cores > 1:
		serpent_command = f'{serpent_executable} -disperse -omp {cores}'
	else:
		serpent_command = f'{serpent_executable}'

	full_command = f"printf '{string_printf}\\n' | {serpent_command}"
	proc = subprocess.Popen(full_command, shell=True, stdout=subprocess.PIPE ,stderr=subprocess.PIPE, universal_newlines=True)
	output = ''
	for line in proc.stdout:
		if print_output:
			sys.stdout.write(line)
		output += line
	proc.wait()
	data = pd.read_csv(output_file_path, header=None, delim_whitespace=True, names=['x','y','z','r','uni'])
	info = output.strip().splitlines()[-6].split()
	N_spheres = int(info[0])
	PF = float(info[-1])
	return data, PF, N_spheres, output
	# Nspheres =
	# n
	#print_with_timestamp(output.decode("utf-8") )

def convert_pos_to_simple_pbed(pos_file_path, radii, universes, path_out=None, translation=[0,0,0], scale_factors=[1,1,1], del_in=',', header=None):
	positions = pd.read_csv(pos_file_path, delimiter=del_in, names=['x','y','z'], header=header)
	positions = np.multiply(positions + translation, scale_factors)
	positions['r'] = radii
	positions['uni'] = universes
	if not isinstance(path_out, type(None)):
		positions.to_csv(path_out, header=False, index=False, sep='\t')
	return positions

def init_case(case_name, python_input, path_to_case, output_name=None):
	if not output_name:
			output_name = case_name
	try:
		shutil.rmtree(f'Cases/{output_name}/')
	except:
		pass
	shutil.copytree(f'{path_to_case}/{case_name}/', f'Cases/{output_name}/')

	print_with_timestamp(f'Creating folder at ./Cases/{output_name}/ and copying input')
	for path_folder in [f'Cases/{output_name}/Plots/', f'Cases/{output_name}/Data/', f'Cases/{output_name}/wrk_Serpent/']:
		os.makedirs(path_folder, exist_ok=True)
	for file in [python_input, 'Source/Default_Input.py']:
		if file[-3:] != '.py':
			file += '.py'
		shutil.copy(file, f'Cases/{output_name}/')
	shutil.copy('Utils/tuto1', f'Cases/{output_name}/')
	if os.path.exists(f'{path_to_case}/{case_name}/OF') and os.path.isdir(f'{path_to_case}/{case_name}/OF'):
		shutil.copytree(f'{path_to_case}/{case_name}/OF', f'Cases/{output_name}/wrk_Serpent/OF')

def reset_case(to_recreate='all'):
	if isinstance(to_recreate, str) and to_recreate=='all':
		to_recreate = ['./Plots/', './Data/', './wrk_Serpent/']
	elif isinstance(to_recreate, str):
		to_recreate = [to_recreate]
	for path_folder in ['./Plots/', './Data/', './wrk_Serpent/']:
		if os.path.exists(path_folder):
			shutil.rmtree(path_folder)
	for path_folder in to_recreate:
			os.makedirs(path_folder)

def erase_working_plots():
	for file in glob('./wrk_Serpent/*.png'):
		os.remove(file)

def get_transferrable(transferrable, serpent_instance=None, input_parameter=False):
	if isinstance(transferrable, str):
		if transferrable[:3] != 'sss':
			if input_parameter:
				transferrable = f'sss_iv_{transferrable}'
			else:
				transferrable = f'sss_ov_{transferrable}'
		tra = serpent_instance.get_transferrable(transferrable)
	else:
		tra = transferrable
	return tra

def Serpent_get(transferrable, serpent_instance=None, input_parameter=False):
	tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
	tra.communicate()
	return tra

def Serpent_get_values(transferrable, serpent_instance=None, input_parameter=False, return_singles=True, communicate=True):
	if communicate:
		tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
	else:
		tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
	#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
	values = tra.value_vec
	if return_singles and len(values) == 1:
		values = values[0]
	return values

def Serpent_get_multiple_values(transferrables_matrix, serpent_instance=None, input_parameter=False, return_singles=True, communicate=True, ignore_None=True):
	transferrables_matrix = np.array(transferrables_matrix)
	values_matrix = np.empty_like(transferrables_matrix)
	dim = transferrables_matrix.ndim
	if dim == 1:
		for i, transferrable in enumerate(transferrables_matrix):
			if not ignore_None or not isinstance(transferrable, type(None)):
				if communicate:
					tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
				else:
					tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
				#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
				values = tra.value_vec
				if return_singles and len(values) == 1:
					values = values[0]
				values_matrix[i] = np.array(values)
	elif dim == 2:
		for i, el0 in enumerate(transferrables_matrix):
			for j, transferrable in enumerate(el0):
				if not ignore_None or not isinstance(transferrable, type(None)):
					if communicate:
						tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
					else:
						tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
					#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
					values = tra.value_vec
					if return_singles and len(values) == 1:
						values = values[0]
					values_matrix[i, j] = np.array(values)
	elif dim == 3:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, transferrable in enumerate(el1):
					if not ignore_None or not isinstance(transferrable, type(None)):
						if communicate:
							tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
						else:
							tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
						#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
						values = tra.value_vec
						if return_singles and len(values) == 1:
							values = values[0]
						values_matrix[i, j, k] = np.array(values)
	elif dim == 4:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, el2 in enumerate(el1):
					for l, transferrable in enumerate(el2):
						if not ignore_None or not isinstance(transferrable, type(None)):
							if communicate:
								tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
							else:
								tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
							#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
							values = tra.value_vec
							if return_singles and len(values) == 1:
								values = values[0]
							values_matrix[i, j, k, l] = np.array(values)
	elif dim == 5:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, el2 in enumerate(el1):
					for l, el3 in enumerate(el2):
						for m, transferrable in enumerate(el3):
							if not ignore_None or not isinstance(transferrable, type(None)):
								if communicate:
									tra = Serpent_get(transferrable, serpent_instance, input_parameter=input_parameter)
								else:
									tra = get_transferrable(transferrable, serpent_instance, input_parameter=input_parameter)
								#print_with_timestamp(f'Getting transferrable values for "{tra.name}"')
								values = tra.value_vec
								if return_singles and len(values) == 1:
									values = values[0]
								values_matrix[i, j, k, l, m] = np.array(values)
	return values_matrix

def Serpent_get_material_wise(parent_name, parameter, serpent_instance, prefix='material', input_parameter=False):
	zones = natural_sort([name for name in Serpent_get_values('materials', serpent_instance) if f'{parent_name}z' in name])
	transferrables = [get_transferrable(f'{prefix}_{name}_{parameter}', serpent_instance, input_parameter=input_parameter) for name in zones]
	return np.array(transferrables)

def Serpent_set_values(transferrable, values, serpent_instance=None, communicate=True):
	tra = get_transferrable(transferrable, serpent_instance, input_parameter=True)
	#print_with_timestamp(f'Setting transferrable "{tra.name}" to values : {values}')
	if isinstance(values, (int, float, np.integer)):
		values = [values]
	tra.value_vec = np.array(values)
	if communicate:
		tra.communicate()
	return tra

def Serpent_set_multiple_values(transferrables_matrix, values_matrix, serpent_instance=None, communicate=True, ignore_None=True):
	transferrables_matrix = np.array(transferrables_matrix)
	values_matrix = np.array(values_matrix)
	dim = transferrables_matrix.ndim
	if dim == 1:
		for i, transferrable in enumerate(transferrables_matrix):
			if not ignore_None or not isinstance(transferrable, type(None)):
				tra = Serpent_set_values(transferrable, values_matrix[i], serpent_instance=serpent_instance)
			#print_with_timestamp(f'Setting transferrable values for "{tra.name}"')
	elif dim == 2:
		for i, el0 in enumerate(transferrables_matrix):
			for j, transferrable in enumerate(el0):
				if not ignore_None or not isinstance(transferrable, type(None)):
					tra = Serpent_set_values(transferrable, values_matrix[i,j], serpent_instance=serpent_instance)
					#print_with_timestamp(f'Setting transferrable values for "{tra.name}"')
	elif dim == 3:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, transferrable in enumerate(el1):
					if not ignore_None or not isinstance(transferrable, type(None)):
						tra = Serpent_set_values(transferrable, values_matrix[i,j,k], serpent_instance=serpent_instance)
					#print_with_timestamp(f'Setting transferrable values for "{tra.name}"')
	elif dim == 4:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, el2 in enumerate(el1):
					for l, transferrable in enumerate(el2):
						if not ignore_None or not isinstance(transferrable, type(None)):
							tra = Serpent_set_values(transferrable, values_matrix[i,j,k,l], serpent_instance=serpent_instance)
						#print_with_timestamp(f'Setting transferrable values for "{tra.name}"')
	elif dim == 5:
		for i, el0 in enumerate(transferrables_matrix):
			for j, el1 in enumerate(el0):
				for k, el2 in enumerate(el1):
					for l, el3 in enumerate(el2):
						for m, transferrable in enumerate(el3):
							if not ignore_None or not isinstance(transferrable, type(None)):
								tra = Serpent_set_values(transferrable, values_matrix[i,j,k,l,m], serpent_instance=serpent_instance)
							#print_with_timestamp(f'Setting transferrable values for "{tra.name}"')
	return values_matrix

def Serpent_set_same_values(list_transferrables, values, serpent_instance=None, communicate=True):
	if isinstance(values, (int, float)):
		values = np.array([values], dtype=type(values))
	for transferrable in list_transferrables:
		tra = get_transferrable(transferrable, serpent_instance, input_parameter=True)
		#print_with_timestamp(f'Setting transferrable "{tra.name}" to user-defined values') # : {values}')
		tra.value_vec = values
		if communicate:
			tra.communicate()
	return list_transferrables

def Serpent_set_option(transferrable, serpent_instance=None, turn_on=True, communicate=True):
	tra = get_transferrable(transferrable, serpent_instance, input_parameter=True)
	if turn_on:
		tra.value_vec[0] = 1
	else:
		tra.value_vec[0] = 0
		#print_with_timestamp(f'Setting transferrable "{tra.name}" to value : {tra.value_vec[0]}')
	if communicate:
		tra.communicate()
	return tra

def assign_random_array(source_list, possible_values, proportions='equal'):
	if isinstance(source_list, int):
		source_list = np.arange(source_list)

	if isinstance(proportions, str) and proportions=='equal':
		shuffled_list = np.array_split(np.random.permutation(source_list), len(possible_values))
	else:
		proportions = np.array(proportions) / np.sum(proportions)
		random_list = np.random.permutation(np.arange(len(source_list)))
		ind = np.add.accumulate(np.array(proportions) * len(random_list)).round(8).astype(int)
		shuffled_list = [x.tolist() for x in np.split(random_list, ind)][:len(proportions)]
	array = np.empty(len(source_list), dtype=np.array(possible_values).dtype)
	for i, val in enumerate(possible_values):
		for j in range(len(shuffled_list[i])):
			array[shuffled_list[i][j]] = val
	return array

def get_index_table(list1, list2):
	index_list = np.ones(len(list1))*np.nan
	for j in range(len(list1)):
		for k in range(len(list2)):
				if list1[j] == list2[k]:
					index_list[j] = k
					break
	return index_list

def estimate_burnup(z, zlim, direction, pass_number, max_passes, max_bu):
	if direction == -1:
		bu = (zlim[1] - z)/(zlim[1] - zlim[0])	* (pass_number/max_passes) * max_bu * (pass_number+1)/pass_number
	elif direction == +1:
		bu = (z - zlim[0])/(zlim[1] - zlim[0])	* (pass_number/max_passes) * max_bu * (pass_number+1)/pass_number
	return bu

def load_interpolator(interpolator_path, to_replace=[666], replace_value=[-1]):
	import pickle
	with open(interpolator_path, 'rb') as f:
		interpolator, zai_interpolator = pickle.load(f)
	zai_interpolator = np.array(zai_interpolator).astype(int)
	zai_interpolator[zai_interpolator == 666] = -1
	return interpolator, zai_interpolator

def interpolate_adens(interpolator, bu, index_table):
	adens_interpolated = interpolator(bu)
	adens = np.ones(len(index_table))*np.nan
	for j in range(len(index_table)):
		if not np.isnan(index_table[j]):
			adens[j] = adens_interpolated[int(index_table[j])]
		else:
			adens[j] = 0.0
	return adens

def clear_all_figures():
	for i in plt.get_fignums():
			plt.figure(i)
			plt.clf()
			plt.close("all")

def ZA_to_name(ZA):
		ZA = str(ZA)
		# Case where ZA (+ temperature card) -> name
		ZA = int(ZA.split('.')[0])
		Z = int(ZA/1000)
		A = int(ZA-Z*1000)
		if A>=300:
				suffix = 'm'
				if Z > 80:
						A -= 100
				else:
						A -= 200
		else:
				suffix = ''
		if A == 0:
				A = 'nat'
		element = Z_to_element(Z)
		name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
		return name

def ZAI_to_name(ZAI):
		# Case where ZAI is given -> name
		if str(ZAI)[-1] == '0':
				ZA = int(ZAI)/10
				suffix = ''
		else:
				ZA = (int(ZAI) - int(str(ZAI)[-1]))/10
				suffix = 'm'

		Z = int(ZA/1000)
		if Z > 120:
				print_with_timestamp('nan')
		element = Z_to_element(Z)
		A = int(ZA-Z*1000)
		if A == 0:
				A = 'nat'
		name = '{}{}{}'.format(element, str(A).zfill(3), suffix)
		return name

def name_to_ZAI(name, separated=False):
		# Case where name is given -> ZAI
		if name[-1] == 'm':
				name = name[:-1]
				I = '1'
		else:
				I = '0'

		if name[-3:] == 'nat':
				element = name[:-3]
				Z = str(element_to_Z(element))
				if not separated:
						return int(Z+'000')
				else:
						return int(Z), 0

		i = len(name)-1
		while name[i:].isdigit():
				i-=1
		A = name[i+1:].zfill(3)
		element = name[:i+1]
		Z = str(element_to_Z(element))
		ZAI = int(Z+A+I)
		if not separated:
				return int(ZAI)
		else:
				return int(Z),int(A),int(I)

def name_to_ZA(name, separated=False):
		# Case where name is given -> ZA
		if name[-1] == 'm':
				name = name[:-1]
				isomer = True
		else:
				isomer = False
		if name[-3:] == 'nat':
				element = name[:-3]
				Z = str(element_to_Z(element))
				if not separated:
						return int(Z+'000')
				else:
						return int(Z), 0

		i = len(name)-1
		while name[i:].isdigit():
				i-=1
		A = name[i+1:].zfill(3)
		element = name[:i+1]
		Z = str(element_to_Z(element))
		if isomer:
				if int(A) < 100:
						A = str(int(A)+200)
				else:
						A = str(int(A)+100)
		ZA = Z+A
		if not separated:
				return int(ZA)
		else:
				return int(Z),int(A)

def Z_to_element(Z):
		elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
		return elements[Z-1]

def element_to_Z(element):
		elements = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt']
		return elements.index(element)+1

def natural_sort(l):
	convert = lambda text: int(text) if text.isdigit() else text.lower()
	alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
	return sorted(l, key=alphanum_key)

def xyzr_to_array(xyzr):
	return xyzr.reshape((-1, 4))

def get_domain_candidates(nnodes, max_domains):
    candidates = []
    def generate_combinations(index, product, combination):
        if index == len(max_domains):
            if product == nnodes:
                candidates.append(combination)
            return
        for i in range(1, max_domains[index]+1):
            if product * i > nnodes:
                break
            generate_combinations(index+1, product*i, combination+[i])
    generate_combinations(0, 1, [])
    return candidates


def nodes_to_dd(nnodes, allowed_decomposition_types, max_domains):
	candidates = get_domain_candidates(nnodes, max_domains)
	print_with_timestamp(candidates)
	print_with_timestamp(max_domains)
	print_with_timestamp(nnodes)
	best_score = 0
	for c in candidates:
		score = np.sum(np.power(c, np.arange(len(c), 0, -1)))
		if score > best_score:
			best_score = float(score)
			best_candidate = list(c)

	if candidates==0:
		raise Exception(f"No configuration in {max_domains} can make exactly {nnodes} nodes.")

	decomposition_types = []
	decomposition_domains = []
	for i in range(len(allowed_decomposition_types)):
		if best_candidate[i] > 1:
			decomposition_types.append(allowed_decomposition_types[i])
			decomposition_domains.append(best_candidate[i])
	return decomposition_types, decomposition_domains

def process_inventory(inventory_names):
	if not inventory_names:
		return None
	if isinstance(inventory_names, str):
		if inventory_names == "all":
			return "all"
		elif inventory_names == "accident":
			return [541330, 531310, 551340, 551370, 521320, 380890, 380900, 561400, 400950, 420990, 441030, 441060, 581410, 581440, 932390, 942380, 942390, 942400, 942410, 962420, -1]
		elif inventory_names == "actinides" or inventory_names == "act":
			return [ 892250, 892260, 892270, 902270, 902280, 902290, 902300, 902320, 902330, 902340, 912310, 912320, 912330, 922320, 922330, 922340, 922350, 922360, 922370, 922380, 932350, 932360, 932370, 932380, 932390, 942360, 942370, 942380, 942390, 942400, 942410, 942420, 942430, 942440, 942460, 952410, 952420, 952421, 952430, 952440, 952441, 962400, 962410, 962420, 962430, 962440, 962450, 962460, 962470, 962480, 962490, 962500, 972470, 972490, 972500, 982490, 982500, 982510, 982520, 992530, 992540, 992550, 1002550, -1]
		elif inventory_names == "burnupcredit":
			return [922340, 922350, 922360, 922380, 942380, 942390, 942400, 942410, 942420, 932370, 952410, 952430, 962430, 962440, 962450, 551330, 601430, 601350, 621470, 621490, 621500, 621510, 621520, 631530, 641550, 420950, 430990, 441010, 451030, 471090, 481130, -1]
		elif inventory_names == "burnupindicators":
			return [551370, 571390, 601480, -1]
		elif inventory_names == "cosi6":
			return [922320, 922330, 922340, 922350, 922360, 922380, 942360, 942380, 942390, 942400, 942410, 942420, 952410, 952421, 952430, 932370, 932390, 962420, 962430, 962440, 962450, 962460, 962470, 962480, 902270, 902280, 902290, 902300, 902310, 902320, 912310, 912330, 551330, 551340, 551350, 551360, 551370, 531290, 461040, 461050, 461060, 461070, 461080, 461090, 461100, 461120, 340760, 340770, 340780, 340790, 340800, 340820, 501150, 501160, 501170, 501180, 501190, 501200, 501210, 501220, 501230, 501240, 501250, 501260, 430990, 400900, 400910, 400920, 400930, 400940, 400950, 400960, 400970, -1]
		elif inventory_names == "lanthanides":
			return [571380, 571390, 571400, 581400, 581410, 581420, 581430, 581440, 591410, 591420, 591430, 601420, 601430, 601440, 601450, 601460, 601470, 601480, 601500, 611470, 611480, 611481, 611490, 611510, 621440, 621470, 621480, 621490, 621500, 621510, 621520, 621530, 621540, 631510, 631520, 631530, 631540, 631550, 631560, 631570, 641520, 641540, 641550, 641560, 641570, 641580, 641600, 651590, 651600, 661600, 661610, 661620, 661630, 661640, 671650, 681640, 681660, 681670, 661680, 661700, 711750, 711760, -1]
		elif inventory_names == "longterm":
			return [10030, 40100, 60140, 170360, 180390, 180420, 190400, 200410, 250530, 250540, 260550, 270600, 280590, 280630, 300650, 340790, 360810, 360850, 370870, 380900, 400930, 410910, 410920, 410931, 410940, 420930, 430970, 430990, 441060, 461070, 471081, 471101, 481090, 481131, 501191, 501211, 501230, 501250, 501260, 521251, 521271, 531290, 551340, 551350, 551370, 561330, 571370, 571380, 581440, 611450, 611470, 621470, 621510, 631520, 631540, 631550, 641530, 671630, 671661, 691700, 691710, 711740, 711760, 721781, 721820, 781930, 812040, 822050, 822100, 832080, 832101, 842100, 882230, 882250, 882260, 882280, 892270, 902270, 902280, 902290, 902300, 902320, 902340, 912310, 912330, 922320, 922330, 922340, 922350, 922360, 922380, 932370, 942360, 942380, 942390, 942400, 942410, 942420, 952410, 952421, 952430, 962420, 952430, 952440, 952450, 952460, 952480, 962490, 962500, 962510, 962520, -1]
		elif inventory_names == "minoractinides":
			return [892250, 892260, 892270, 912310, 912320, 912330, 932350, 932360, 932370, 932380, 932390, 952410, 952420, 952421, 952430, 952440, 952441, 962400, 962410, 962420, 962430, 962440, 962450, 962460, 962470, 962480, 962490, 962500, 972470, 972490, 972500, 982490, 982500, 982510, 982520, 992530, 992540, 992550, 1002550, -1]
		elif inventory_names == "dp":
			return [802060, 812060, 812060, 812070, 812080, 812090, 812100, 822060, 822060, 822070, 822070, 822080, 822080, 822090, 822090, 822100, 822100, 822110, 822120, 822140, 832090, 832100, 832110, 832120, 832130, 832140, 832140, 842100, 842110, 842120, 842130, 842130, 842140, 842140, 842150, 842160, 842180, 852170, 852180, 862170, 862180, 862190, 862200, 862220, 872210, 872230, 882230, 882230, 882240, 882250, 882260, 882270, 882280, 892250, 892270, 892280, -1]
		else:
			return [inventory_names]
	else:
		return inventory_names

def getdict(dictionary, key=None):
	if key:
		return np.array(list(element.get(key) for element in dictionary.values()))
	else:
		return np.array(list(dictionary.keys()))
#data, PF, N, output = use_disperser('cylinder', (120.00, 60.00, 369.47), 0.6, 2, 1, 'fpb_pos', 64, grow_and_shake=True, shake_factor=0.05, growth_factor=0.05, print_output=True)

# %%
