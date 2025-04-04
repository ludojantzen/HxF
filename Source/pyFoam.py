########## 
#### 1D channel code for the gFHR. It reads the powerDensityNeutronics from Serpent 
#### tally in a mesh and solve for Coolant Velocity, Presure, Coolant Temperature, Pebble Temperature (surface, matrix, fuel)
#### Use MPI to speed up the process 
#### Return vtk files, tally in the same mesh tha the powerDensity
#### To do: check the possibility 





###module
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from mpi4py import MPI
from scipy.interpolate import interp1d
import joblib
from joblib import Parallel, delayed
from joblib import parallel_backend
import time, math
import psutil
import sys

### read the step to read the power and add the step size to have the file to write the resukts 
# step=int(sys.argv[1])
# step_size=25

def read_fields(path):
    # Read cell volumes
    with open(path, 'r') as f:
        lines = f.readlines()

    # Find the starting and ending indices of the block
    start_idx = next(i for i, line in enumerate(lines) if '(' in line)
    end_idx = next(i for i, line in enumerate(lines[start_idx:]) if ';' in line) + start_idx - 1

    # Extract the block of lines
    block = lines[start_idx + 1:end_idx]

    # Use NumPy to efficiently process and store the data
    return np.array([float(line) for line in block])
    
def read_data_mesh(step):
    data_mesh=pd.DataFrame()
    data_mesh['x']=read_fields('./OF/0/Cx')
    data_mesh['y']=read_fields('./OF/0/Cy')
    data_mesh['z']=read_fields('./OF/0/Cz')
    data_mesh['V']=read_fields('./OF/0/V')
    data_mesh['r_dist']=(data_mesh['x']**2+data_mesh['y']**2)**0.5
    data_mesh['powerDensityNeutronics']=read_fields(f'./OF/{step}/fluidRegion/powerDensityNeutronics')
    return data_mesh



#### Block to compute velocity and Presure for a given fluid temperature profiel ###
###################################
#### Has to be set up: U_inlet, P_outlet and solver parameters
#### Based on pimplefoam with inner iteration for the pressure correction 
# #### Can set up the physics you want to simulate  on_p = {
#         "advection": True,
#         "diffusion": True,
#         "pressure_calc": True,
#         "p_correction": True,
#         "gravity": True,
#         "pressure_drop": True
#     }

#### Pressure drop based on ergun formula, tuned to match INL pressure drop 


def compute_courant_number(u, dt, dx):
    return np.max(np.abs(u)) * dt / dx

def pimple(Tf,nx):
    # iterate through time
    velocity = []
    pressure = []
    rel_var_u = []
    rel_var_p = []
    time_step = []

    Lx=3.0947
    # nx = 200  # Number of cells
    nt = 5000  # nt is the number of timesteps
    dt = 1e-4  # Initial timestep (delta t)
    maxit = 1000  # Pressure solver maximum iteration
    tol = 1e-5  # Pressure solver maximum error
    tol_conv = 1e-3  # Velocity Solver maximum error
    plotfreq = nt // 50  # Output plot frequency
    porosity = 1 - 0.5988494068114064
    # porosity_list = np.array([0.4+0.1*np.sin(2*k*np.pi/nx) for k in range(nx+1)])
    porosity_list=np.array([porosity]*(nx+1))
    U_inlet = 1173 / ((2413.03 - 0.4884 * 823) * np.pi * 1.2 ** 2)
    P_outlet = 2e5  # Pa
    dp = 4e-2  # pebble diameter in m


    # Set constant of the problem
    g = 9.81  # m/s^2

    # CFL parameters
    cfl_target = 0.001
    dt_min = 1e-7
    dt_max = 10
    rho = -0.4884 * Tf + 2413.0
    eta = 1.16e-4 * np.exp(3755.0 / Tf)
    Cp = np.array([2416] * (nx + 1))
    k_th = 5.0e-4 * Tf + 0.63

    # Switch bits of the solver on
    on_p = {
        "advection": True,
        "diffusion": True,
        "pressure_calc": True,
        "p_correction": True,
        "gravity": True,
        "pressure_drop": True
    }

    # Initialise fields
    dx = Lx / (nx - 1.)  # Cell size
    u = np.zeros(nx + 1)
    p = np.zeros(nx)

    u[0] = U_inlet
    p[-1] = P_outlet

    divergence_limit = 1e10  # Divergence criterion

    for n in range(nt):
        # Compute the current Courant number and adjust the timestep
        courant_number = compute_courant_number(u, dt, dx)
        if n > 5:
            dt = np.clip(cfl_target * dx / np.max(np.abs(u)), dt_min, dt_max)
        time_step.append(dt)

        # Check for divergence
        if np.max(np.abs(u)) > divergence_limit or np.max(np.abs(p)) > divergence_limit:
            print("Divergence detected, stopping simulation.")
            break

        # Advection and diffusion to project pressure free velocity
        i = np.arange(1, nx)
        un = u.copy()
        if on_p["advection"]:
            advection_term = dt / dx * (un[i] * (un[i] / porosity_list[i] - un[i-1] / porosity_list[i-1]))
            u[i] -= advection_term
        if on_p["diffusion"]:
            diffusion_term = (dt / (rho[i])) * (eta[i+1] * (un[i+1] / porosity_list[i+1] - un[i] / porosity_list[i]) - eta[i] * (un[i] / porosity_list[i] - un[i-1] / porosity_list[i-1])) / dx ** 2
            u[i] += diffusion_term
        if on_p['gravity']:
            gravity_term = g * dt * porosity_list[i]
            u[i] -= gravity_term
        if on_p['pressure_drop']:
            if n > 50:
                Re = rho * un * dp / eta
                W =  (150 * (1 - porosity_list[i]) / Re[i] + 1.75) * ((1 - porosity_list[i]) / porosity_list[i]) * un[i] / dp
                pressure_drop_term = W * un[i] * dt / porosity_list[i]
                u[i] -= pressure_drop_term / porosity_list[i]

        # Get Pressure using incompressible continuity eqn
        if on_p["pressure_calc"]:
            for it in range(maxit):
                pn = p.copy()
                i = np.arange(1, nx-1)
                p[i] = 1/(porosity_list[i]+porosity_list[i-1])*(porosity_list[i]*p[i+1]+porosity_list[i-1]*p[i-1]-dx/dt*(rho[i] * u[i] - rho[i-1] * u[i-1]))

                # Pressure boundary conditions
                p[-1] = P_outlet
                p[0] = p[1]

                error = abs(pn - p).max()
                if error < tol:
                    break

        # Correct pressure free field
        if on_p["p_correction"]:
            i = np.arange(1, nx-1)
            correction_term = (dt*porosity_list[i]) / (rho[i] * dx) * (p[i+1] - p[i])
            u[i] -= correction_term

        # Velocity boundary conditions
        u[0] = U_inlet
        u[-1] = U_inlet*rho[0]/(rho[-1])
        u[-2] = U_inlet*rho[0]/(rho[-2])
        # u[-3] = U_inlet*rho[0]/(rho[-3])

        if n > 0:  # Skip the first step for relative variation calculation
            rel_var_u.append(np.max(np.abs((u - un) / (un + 1e-10))) * 100)
            if on_p['p_correction']==True:
                rel_var_p.append(np.max(np.abs ((p - pn) / (pn + 1e-10))) * 100)
        if n > 0 and rel_var_u[-1] < tol_conv:
            break

        if n % plotfreq == 0:
            velocity.append(u.copy())
            pressure.append(p.copy())



    # Plot relative variations
    # plt.figure()
    # plt.plot(rel_var_u, label='Relative variation of u')
    # plt.plot(rel_var_p, label='Relative variation of p')
    # plt.xlabel('Time step')
    # plt.ylabel('Relative variation %')
    # plt.legend()
    # plt.grid()
    # plt.yscale('log')
    # plt.title('Relative Variation of Fields at Each Time Step')

    # # Plot final velocity and pressure
    # plt.figure()
    # plt.plot(u, label='Final velocity')
    # plt.xlabel('Cell index')
    # plt.ylabel('Velocity (m/s)')
    # plt.legend()
    # plt.grid()
    # # plt.yticks(format='plain')
    # plt.title('Final Velocity Field')

    # plt.figure()
    # plt.plot(u/porosity_list, label='Final real velocity')
    # plt.xlabel('Cell index')
    # plt.ylabel('Velocity (m/s)')
    # plt.legend()
    # plt.grid()
    # # plt.yticks(format='plain')
    # plt.title('Final Velocity Field')

    # plt.figure()
    # plt.plot(p, label='Final pressure')
    # plt.xlabel('Cell index')
    # plt.ylabel('Pressure (Pa)')
    # plt.legend()
    # plt.grid()
    # plt.title('Final Pressure Field')

    # plt.figure()
    # plt.plot(time_step, label='timestep')
    # plt.xlabel('Time step index')
    # plt.ylabel('Time step (s)')
    # plt.legend()
    # plt.grid()
    # plt.yscale('log')

    print(f'Velocity fields converged with a maximum relative residual of {rel_var_u[-1]} %')
    print(f'Pressure fields converged with a maximum relative residual of {rel_var_p[-1]} %')
    return u, p 




#### Block to compute fluid and solid temperature for a given power and fluid velocity
###################################
#### Has to be set up: T_inlet
#### Based on explicit Euler scheme with inner iteration for the solid temperature (pebble surface temperature)
# Can set up the physics for both solid phase and liquid phase
    # on_f = {
    #     "advection": True,
    #     "diffusion": True,
    #     "interphase": True,
    #     "power_source": True,
    # }

    # on_s = {
    #     "diffusion": True,
    #     "interphase": True,
    #     "power_source": True,
    # }
    


### within the same loop:
def energy(u,power,nx):
    # Meshing
    Lx = 3.0947  # Domain size
    # nx = 100  # Number of cells
    nt = 50000  # Initial number of timesteps
    dt = 1e-2  # Initial timestep (delta t)
    dx = Lx / (nx - 1.)  # Cell size
    #Constant
    PF = 0.5988494068114064 # Example porosity fraction
    porosity = 1 - PF
    porosity_list = np.array([porosity] * (nx+1))
    # pow_dens=280e6/(3.0947*np.pi*1.2**2)
    dp = 4e-2  # pebble diameter in m
    #convergence
    tol_conv=1e-4
    # Set Initial and Boundary Conditions
    Tf_inlet = 823
    Tf = np.zeros(nx + 1)
    Tf[:] = Tf_inlet
    Ts_guess=950.0
    Ts = np.array([Ts_guess] * (nx+1))

    # qf=np.array([0]*(nx+1))
    qs=power
    qf=np.array([0]*(nx+1))
    
    # Switch bits of the solver on
    on_f = {
        "advection": True,
        "diffusion": True,
        "interphase": True,
        "power_source": True,
    }

    on_s = {
        "diffusion": True,
        "interphase": True,
        "power_source": True,
    }
    
    for n in range(nt):
        i = np.arange(1, nx)
        # Generate coolant properties
        rho = -0.4884 * Tf + 2413.0
        eta = 1.16e-4 * np.exp(3755.0 / Tf)
        Cp = np.array([2416] * (nx+1))
        k_th = 5.0e-4 * Tf + 0.63
        # Generate htc 
        Re = rho * u * dp / eta
        Pr = eta * Cp / k_th
        Nu= 0.357*Pr**0.33/porosity_list*Re**0.641
        alpha = 6 * (1 - porosity_list) / dp * Nu * k_th / dp
        # Advection and diffusion to project pressure free velocity
        Tfn = Tf.copy()
        if on_f["advection"]:
            advection_term = dt / (dx*porosity_list[i]) * (u[i] * (Tfn[i] - Tfn[i-1]))
            Tf[i] -= advection_term
        if on_f["diffusion"]:
            diffusion_term = (dt / (rho[i] * Cp[i] * porosity_list[i])) * (k_th[i+1] * Tfn[i+1] - 2. * k_th[i] * Tfn[i] + k_th[i-1] * Tfn[i-1]) / dx**2
            Tf[i] += diffusion_term
        if on_f['interphase']:
            interphase_term = (dt / (rho[i] * Cp[i] * porosity_list[i])) * alpha[i] * (Ts[i]-Tfn[i])
            Tf[i] += interphase_term
        if on_f['power_source']:
            source_term = (dt / (rho[i] * Cp[i] * porosity_list[i])) * qf[i]
            Tf[i] += source_term
        
        # Velocity boundary conditions
        Tf[0] = Tf_inlet
        Tf[-1] = Tf[-2]
        #### find Ts from this Tf solution  
        rho_f = -0.4884*Tf+2413.0
        rho_s=np.array([1632]*nx)
        Cp_s=4184*(0.5421-2.4e-6*Ts-90.273/Ts-43449/Ts**2+1.59e7/Ts**3-1.4369/Ts**4)
        eta = 1.16e-4*np.exp(3755.0/Tf)
        Cp_f=np.array([2416]*(nx+1))
        k_th_f=5.0e-4*Tf + 0.63
        k_th_s=-22.05679*np.log(Ts)+194.32788
        ### generate htc 
        Re = rho_f * u * dp/ eta
        Pr=eta*Cp_f/k_th_f
        ###Wakao
        Nu= 0.357*Pr**0.33/porosity_list*Re**0.641
        alpha=6*(1-porosity_list)/dp*Nu*k_th_f/dp
        # Advection and diffusion to project pressure free velocity
        Tsn = Ts.copy()
        if on_s['interphase']:
            interphase_term = (dt/(rho_s[i]*Cp_s[i]*(1-porosity_list[i])))*alpha[i]*(Tf[i]-Tsn[i])
            Ts[i] += interphase_term
        if on_s['power_source']:
            source_term=(dt/(rho_s[i]*Cp_s[i]*(1-porosity_list[i])))*qs[i]
            Ts[i] += source_term
        if on_s['diffusion']:
            diffusion_term = (dt/(rho_s[i]*Cp_s[i]*(1-porosity_list[i]))) * (k_th_s[i+1]*Tsn[i+1] - 2. *k_th_s[i]*Tsn[i] + k_th_s[i-1]*Tsn[i-1]) / dx**2
            Ts[i] += diffusion_term

        # Velocity boundary conditions
        Ts[0] = Ts[1]
        Ts[-1] = Ts[-2]

        if n > 0:  # Skip the first step for relative variation calculation
            var_s=np.max(np.abs((Ts - Tsn) / (Tsn + 1e-10))) * 100
            var_f=np.max(np.abs((Tf - Tfn) / (Tfn + 1e-10))) * 100
            if n>0 and var_s<tol_conv and var_f<tol_conv:
                print(f'Tf fields converged with a maximum residual of {np.max(Tfn-Tf)} K')
                print(f'Ts fields converged with a maximum residual of {np.max(Tsn-Ts)} K')
                break

    return Tf,Ts

# Example velocity field
# Tf,Ts=energy(u,interpolated_data)
# plt.plot(Tf,label='Coolant Temperature')
# plt.plot(Ts,label='Surface Temperature')
# plt.grid()
# plt.legend()



#### Combine pimple and energy solver 
#### Stops when u and Tf converges, can tune the tolerance. 
#### Needs the power input form powerDensityNeutronics
#### Guess Tf at the first step based on gFHR charactersitics 

def compute1D(power,nx):
    Tf=np.linspace(823,923,nx+1)
    U=[]
    P=[]
    TF=[]
    TS=[]
    tol_conv=1e-3
    for k in range(10):
        # print(f"Outter Iteration number {k+1}")
        # print(f"Runing Pimple Solver")
        u,p=pimple(Tf,nx)
        un=u.copy()
        pn=p.copy()
        Tf,Ts=energy(u,power,nx)
        Tfn=Tf.copy()
        Tsn=Ts.copy()
        if k>0:
            var_u=np.max(np.abs((u - un) / (un + 1e-10))) 
            var_p=np.max(np.abs((p - pn) / (pn + 1e-10))) 
            var_s=np.max(np.abs((Ts - Tsn) / (Tsn + 1e-10)))
            var_f=np.max(np.abs((Tf - Tfn) / (Tfn + 1e-10)))
        if k>0 and var_u<tol_conv:
            break
    return Tf,Ts

#### Mesh interpolation between powerDensityNeutronics and pimple/energy mesh 
def interp(data,len_out):
    original_indices = np.linspace(0, 1, len(data))
    new_indices = np.linspace(0, 1,len_out)
    interpolator = interp1d(original_indices, data, kind='cubic')
    interpolated_data = interpolator(new_indices)
    return interpolated_data



# def compute_3D(data_mesh,k,nx):
#     # print(f'Channel Number {k+1}')
#     z_list=list(set(data_mesh['z']))
#     x_list=data_mesh.loc[data_mesh['z']==z_list[0],'x'].values
#     y_list=data_mesh.loc[data_mesh['z']==z_list[0],'y'].values
#     data_c=data_mesh[(data_mesh['x']==x_list[k]) & (data_mesh['y']==y_list[k])]
#     power=data_c['powerDensityNeutronics']
#     idx=data_c.index
#     power_fin = interp(power, nx+1)
#     Tf, Ts = compute1D(power_fin,nx)
#     Tf_f = interp(Tf, len(z_list))
#     Ts_f = interp(Ts, len(z_list))
#     return idx, Tf_f, Ts_f

    
def write_field_file_T(field_name, num_cells, field_values, file_name):
    # Open the output file for writing
    with open(file_name, 'w') as f:
        # Write the header
        f.write('/*--------------------------------*- C++ -*----------------------------------*\\\n')
        f.write('| =========                 |                                                 |\n')
        f.write('| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n')
        f.write('|  \\\\    /   O peration     | Version:  v2012                                 |\n')
        f.write('|   \\\\  /    A nd           | Website:  www.OpenFOAM.org                      |\n')
        f.write('|    \\\\/     M anipulation  |                                                 |\n')
        f.write('\\*---------------------------------------------------------------------------*/\n')
        f.write('FoamFile\n')
        f.write('{\n')
        f.write('    version     2.0;\n')
        f.write('    format      ascii;\n')
        f.write('    class       volScalarField;\n')
        f.write('    location    "25/fluidRegion";\n')
        f.write('    object      ' + field_name + ';\n')
        f.write('}\n')
        ##### To be modiefied depending on what is the unit of your field 
        #### Power [1 -3 0 0 0 0 0]
        ##### Temperature [0 0 0 1 0 0 0]
        ##### density [1 -3 0 0 0 0 0]
        f.write('dimensions [0 0 0 1 0 0 0];\n')
        f.write('internalField   nonuniform List<scalar>\n')
        f.write(str(num_cells)+'\n')
        f.write('(\n')
        # Write the field values
        for i in range(num_cells):
            f.write(str(field_values[i]) + '\n')
        f.write(')\n')
        f.write(';\n')
        f.write('boundaryField\n')
        f.write('{\n')
        f.write('    top\n')
        f.write('    {\n')
        f.write('        type            zeroGradient;\n')
        f.write('    }\n')
        f.write('    wall\n')
        f.write('    {\n')
        f.write('        type            zeroGradient;\n')
        f.write('    }\n')
        f.write('    bottom\n')
        f.write('    {\n')
        f.write('        type            zeroGradient;\n')
        f.write('    }\n')
        f.write('}\n')

def write_data_mesh(data_mesh,step,step_size):
    #### Pebble characteristics for 1D conduction model within the pebble 
    PF=250190*4/3*np.pi*(2e-2)**3/(np.pi*1.2**2*3.0947)
    V_peb=4/3*np.pi*(2e-2)**3
    Dh=4e-2*(1-PF)/PF
    cpCoolant=2386
    kCoolant=1.1
    
    
    pebbleCoreRadius=1.38e-2
    pebbleMatrixRadius=1.8e-2
    pebbleShellRadius=2e-2
    trisoFuelRadius         =212.5e-6
    trisoBufferRadius       =312.5e-6
    trisoInnerPyCRadius     =352.5e-6
    trisoSiCRadius          =387.5e-6
    trisoOuterPyCRadius     =427.5e-6
    
    nTRISO                  =9022
    
    trisoFuelKCoeffs      =3.3073
    trisoBufferKCoeffs    =0.50
    trisoPyCKCoeffs       =4.00
    trisoSiCKCoeffs       =90.3
    pebbleGraphiteKCoeffs =41.68 
    
    trisoFuelDensityCoeffs      =11031.432 
    trisoBufferDensityCoeffs    =1000
    trisoPyCDensityCoeffs       =1900
    trisoSiCDensityCoeffs       =4.015E+02
    pebbleGraphiteDensityCoeffs =8.12E+01
    
    
    ktriso = (1/trisoFuelRadius-1/trisoOuterPyCRadius)/((1/trisoBufferKCoeffs)*(1/trisoFuelRadius-1/trisoBufferRadius)+(1/trisoPyCKCoeffs)*(1/trisoBufferRadius-1/trisoInnerPyCRadius)+(1/trisoSiCKCoeffs)*(1/trisoInnerPyCRadius-1/trisoSiCRadius)+1/trisoPyCKCoeffs*(1/trisoSiCRadius-1/trisoOuterPyCRadius))
    PF_triso=nTRISO*4/3*np.pi*trisoOuterPyCRadius**3/(4/3*np.pi*(pebbleMatrixRadius**3-pebbleCoreRadius**3))
    kappa=ktriso/pebbleGraphiteKCoeffs
    beta=(kappa-1)/(kappa+2)
    keff=pebbleGraphiteKCoeffs*(1+2*beta*PF_triso)/(1-beta*PF_triso)
    
    
    data_mesh['Q_p'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF
    data_mesh['Q_p_dens'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF / V_peb
    data_mesh['Q_p_surf'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF / (4 * np.pi * pebbleShellRadius**2)
    data_mesh['Tshell'] = data_mesh['Tps'] + data_mesh['Q_p'] / (np.pi * 4 * pebbleGraphiteKCoeffs) * (1 / pebbleMatrixRadius - 1 / pebbleShellRadius)
    data_mesh['Q_m_dens'] = data_mesh['Q_p'] / (4/3 * np.pi * (pebbleMatrixRadius**3 - pebbleCoreRadius**3))
    data_mesh['Tmav'] = data_mesh['Tshell'] + data_mesh['Q_m_dens'] / (6 * keff) * (pebbleMatrixRadius**2 - 3/5 * (pebbleMatrixRadius**5 - pebbleCoreRadius**5) / (pebbleMatrixRadius**3 - pebbleCoreRadius**3)) + data_mesh['Q_m_dens'] * pebbleCoreRadius**3 / (3 * keff) * (1 / pebbleMatrixRadius - 3/2 * (pebbleMatrixRadius**2 - pebbleCoreRadius**2) / (pebbleMatrixRadius**3 - pebbleCoreRadius**3))
    data_mesh['Q_triso'] = data_mesh['Q_p'] / nTRISO
    data_mesh['Q_fuel_dens'] = data_mesh['Q_triso'] / (4/3 * np.pi * (trisoOuterPyCRadius**3))
    data_mesh['Tfs'] = data_mesh['Tmav'] + data_mesh['Q_triso'] / (4 * np.pi) * ((1 / trisoBufferKCoeffs) * (1 / trisoFuelRadius - 1 / trisoBufferRadius) + (1 / trisoPyCKCoeffs) * (1 / trisoBufferRadius - 1 / trisoInnerPyCRadius) + (1 / trisoSiCKCoeffs) * (1 / trisoInnerPyCRadius - 1 / trisoSiCRadius) + 1 / trisoPyCKCoeffs * (1 / trisoSiCRadius - 1 / trisoOuterPyCRadius))
    data_mesh['Tfav'] = data_mesh['Tfs'] + data_mesh['Q_fuel_dens'] * trisoFuelRadius**2 / (15 * trisoFuelKCoeffs)
    data_mesh['Tfmax'] = data_mesh['Tfs'] + data_mesh['Q_fuel_dens'] * trisoFuelRadius**2 / (6 * trisoFuelKCoeffs)
    # data_mesh.to_csv('data_mesh.csv', index=False)
    # print('Write the results to the time 25/fluidRegion')
    write_field_file_T('T', len(data_mesh['T']), data_mesh['T'], f'./OF/{step+step_size}/fluidRegion/T')
    write_field_file_T('Tps', len(data_mesh['Tps']), data_mesh['Tps'], f'./OF/{step+step_size}/fluidRegion/Tps')
    write_field_file_T('Tmav.nuclearSteadyStatePebble', len(data_mesh['Tmav']), data_mesh['Tmav'], f'./OF/{step+step_size}/fluidRegion/Tmav.nuclearSteadyStatePebble')
    write_field_file_T('Tfav.nuclearSteadyStatePebble', len(data_mesh['Tfav']), data_mesh['Tfav'], f'./OF/{step+step_size}/fluidRegion/Tfav.nuclearSteadyStatePebble')
    write_field_file_T('Tfmax', len(data_mesh['Tfmax']), data_mesh['Tfmax'], f'./OF/{step+step_size}/fluidRegion/Tfmax')


def compute_task(data_mesh,channel):
    return compute_3D(data_mesh, channel, 200)
    
    
def compute_3D(data_mesh,k,nx):
    # print(f'Channel Number {k+1}')
    z_list=list(set(data_mesh['z']))
    x_list=data_mesh.loc[data_mesh['z']==z_list[0],'x'].values
    y_list=data_mesh.loc[data_mesh['z']==z_list[0],'y'].values
    data_c=data_mesh[(data_mesh['x']==x_list[k]) & (data_mesh['y']==y_list[k])]
    power=data_c['powerDensityNeutronics']
    idx=data_c.index
    power_fin = interp(power, nx+1)
    Tf, Ts = compute1D(power_fin,nx)
    Tf_f = interp(Tf, len(z_list))
    Ts_f = interp(Ts, len(z_list))
    return idx, Tf_f, Ts_f


# if __name__ == "__main__":
#     comm = MPI.COMM_WORLD
#     rank = comm.Get_rank()
#     size = comm.Get_size()
#     start_time = time.time()

#     # Assuming `data_mesh` is a preloaded DataFrame with the required structure
#     # data_mesh = pd.read_csv('./data_mesh.csv')  # Load your data mesh here
#     z_list = list(set(data_mesh['z']))
#     x_list = data_mesh.loc[data_mesh['z'] == z_list[0], 'x'].values
#     y_list = data_mesh.loc[data_mesh['z'] == z_list[0], 'y'].values
#     data_mesh['T']=823
#     data_mesh['Tps']=860

#     channel_list = list(range(len(x_list)))
#     # channel_list = list(range(390))
#     chunk_size = len(channel_list) // size
#     start_index = rank * chunk_size
#     end_index = start_index + chunk_size

#     if rank == size - 1:
#         end_index = len(channel_list)

#     local_channel_list = channel_list[start_index:end_index]

#     # Number of available CPUs on the node
#     num_cpus = psutil.cpu_count(logical=False)
    
#     print(f'Process {rank} is runing on {num_cpus} processors')

#     # Parallel processing with joblib
#     results = Parallel(n_jobs=num_cpus)(delayed(compute_3D)(k,200) for k in local_channel_list)

#     # Gather results from all nodes
#     all_results = comm.gather(results, root=0)

#     if rank == 0:
#         final_results = [result for sublist in all_results for result in sublist]
#         for k in range(len(final_results)):
#             data_mesh.loc[final_results[k][0], 'T'] = final_results[k][1]
#             data_mesh.loc[final_results[k][0], 'Tps'] = final_results[k][2]

#         print('Run 1D conduction model to solve for pebble temperature')
#         data_mesh['Q_p'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF
#         data_mesh['Q_p_dens'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF / V_peb
#         data_mesh['Q_p_surf'] = data_mesh['V'] * data_mesh['powerDensityNeutronics'] / PF / (4 * np.pi * pebbleShellRadius**2)
#         data_mesh['Tshell'] = data_mesh['Tps'] + data_mesh['Q_p'] / (np.pi * 4 * pebbleGraphiteKCoeffs) * (1 / pebbleMatrixRadius - 1 / pebbleShellRadius)
#         data_mesh['Q_m_dens'] = data_mesh['Q_p'] / (4/3 * np.pi * (pebbleMatrixRadius**3 - pebbleCoreRadius**3))
#         data_mesh['Tmav'] = data_mesh['Tshell'] + data_mesh['Q_m_dens'] / (6 * keff) * (pebbleMatrixRadius**2 - 3/5 * (pebbleMatrixRadius**5 - pebbleCoreRadius**5) / (pebbleMatrixRadius**3 - pebbleCoreRadius**3)) + data_mesh['Q_m_dens'] * pebbleCoreRadius**3 / (3 * keff) * (1 / pebbleMatrixRadius - 3/2 * (pebbleMatrixRadius**2 - pebbleCoreRadius**2) / (pebbleMatrixRadius**3 - pebbleCoreRadius**3))
#         data_mesh['Q_triso'] = data_mesh['Q_p'] / nTRISO
#         data_mesh['Q_fuel_dens'] = data_mesh['Q_triso'] / (4/3 * np.pi * (trisoOuterPyCRadius**3))
#         data_mesh['Tfs'] = data_mesh['Tmav'] + data_mesh['Q_triso'] / (4 * np.pi) * ((1 / trisoBufferKCoeffs) * (1 / trisoFuelRadius - 1 / trisoBufferRadius) + (1 / trisoPyCKCoeffs) * (1 / trisoBufferRadius - 1 / trisoInnerPyCRadius) + (1 / trisoSiCKCoeffs) * (1 / trisoInnerPyCRadius - 1 / trisoSiCRadius) + 1 / trisoPyCKCoeffs * (1 / trisoSiCRadius - 1 / trisoOuterPyCRadius))
#         data_mesh['Tfav'] = data_mesh['Tfs'] + data_mesh['Q_fuel_dens'] * trisoFuelRadius**2 / (15 * trisoFuelKCoeffs)
#         data_mesh['Tfmax'] = data_mesh['Tfs'] + data_mesh['Q_fuel_dens'] * trisoFuelRadius**2 / (6 * trisoFuelKCoeffs)
#         data_mesh.to_csv('data_mesh.csv', index=False)
#         print('Write the results to the time 25/fluidRegion')
#         write_field_file_T('T', len(data_mesh['T']), data_mesh['T'], './{step+step_size}/fluidRegion/T')
#         write_field_file_T('Tps', len(data_mesh['Tps']), data_mesh['Tps'], './{step+step_size}//25/fluidRegion/Tps')
#         write_field_file_T('Tmav', len(data_mesh['Tmav']), data_mesh['Tmav'], './{step+step_size}//25/fluidRegion/Tmav')
#         write_field_file_T('Tfav', len(data_mesh['Tfav']), data_mesh['Tfav'], './{step+step_size}//25/fluidRegion/Tfav')
#         write_field_file_T('Tfmax', len(data_mesh['Tfmax']), data_mesh['Tfmax'], './{step+step_size}//25/fluidRegion/Tfmax')

#         end_time = time.time()
#         elapsed_time = end_time - start_time
#         print(f"Execution time: {elapsed_time} seconds for {len(channel_list)} channels and {size} nodes")

#     cpu_load = psutil.cpu_percent()
#     print(f"Process {rank}: CPU load: {cpu_load}%")