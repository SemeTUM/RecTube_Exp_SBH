import numpy as np
import os
import matplotlib.pyplot as plt

##-----------------------------##
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
##-----------------------------##

def clean_complex_line(line):
    return line.strip().replace('(', '').replace(')', '')
    
num_trans_ref= np.loadtxt('./Simulation_Trans_Ref/TransCoefficient_asym_sym_data.txt',skiprows=5)
#print(os.getcwd())

exp_trans_ref_sym = np.loadtxt('./Exp_Trans_Ref/Transcoef_both_sym_glass.txt', dtype=complex, converters={0: lambda s: complex(clean_complex_line(s.decode()))})
sym = np.array(exp_trans_ref_sym)

exp_trans_ref_asym = np.loadtxt('./Exp_Trans_Ref/Transcoef_both_asym_glass.txt', dtype=complex, converters={0: lambda s: complex(clean_complex_line(s.decode()))})
asym = np.array(exp_trans_ref_asym)

plt.figure(figsize=(6,4), dpi=500)
plt.plot(np.arange(0,len(sym)), abs(sym), label='sym (Exp)', color='black', linewidth=2.5)
plt.plot(num_trans_ref[:,0], num_trans_ref[:,2], label='sym (num)', color='black',  marker='o',markeredgecolor='black', linewidth=2.5)
plt.plot(np.arange(0,len(sym)+1), abs(asym), label='asym (Exp)', color= 'lightgray', linewidth=2.5)
plt.plot(num_trans_ref[:,0], num_trans_ref[:,1], label='asym (num)', color='lightgray', marker='o',markeredgecolor='lightgray', linewidth=2.5)
plt.xlim((300, 3600))
plt.xlabel("Frequency [Hz]", fontsize=16)
plt.ylabel("$|T|$", fontsize=16)
plt.tick_params(labelsize=16, width=1.2)
plt.legend(fontsize=16)
plt.grid()
plt.tight_layout()
