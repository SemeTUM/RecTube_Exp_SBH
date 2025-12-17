import numpy as np
from scipy.signal import resample
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt

sym_real = np.loadtxt('./Num_center_pressure/sym_realPresdata.txt', skiprows=8)
sym_imag = np.loadtxt('./Num_center_pressure/sym_imagPresdata.txt', skiprows=8)

asym_real =np.loadtxt('./Num_center_pressure/asym_realPresdata.txt', skiprows=8)
asym_imag = np.loadtxt('./Num_center_pressure/asym_imagPresdata.txt', skiprows=8)

#----symmetric--------
sym_real_1500 = sym_real[:311,1]
sym_real_2270 = sym_real[311:622,1]
sym_real_3000 = sym_real[622:,1]

sym_imag_1500 = sym_imag[:311,1]
sym_imag_2270 = sym_imag[311:622,1]
sym_imag_3000 = sym_imag[622:,1]

sym_abs_1500 = np.sqrt(sym_imag_1500**2 + sym_real_1500**2)  
sym_abs_2270 = np.sqrt(sym_imag_2270**2 + sym_real_2270**2)  
sym_abs_3000 = np.sqrt(sym_imag_3000**2 + sym_real_3000**2)  

#----asymmetric--------
asym_real_1500 = asym_real[:337,1]
asym_real_2270 = asym_real[337:674,1]
asym_real_3000 = asym_real[674:,1]

asym_imag_1500 = asym_imag[:337,1]
asym_imag_2270 = asym_imag[337:674,1]
asym_imag_3000 = asym_imag[674:,1]

#print(asym_imag[:337,0])

asym_abs_1500 = np.sqrt(asym_imag_1500**2 + asym_real_1500**2)  
asym_abs_2270 = np.sqrt(asym_imag_2270**2 + asym_real_2270**2)  
asym_abs_3000 = np.sqrt(asym_imag_3000**2 + asym_real_3000**2)  


def adjust_resolution(newdata):
    x_num = np.linspace(0,5, len(newdata))
    f = interp1d(x_num, newdata, kind='linear')
    xf = np.linspace(x_num.min(), x_num.max(), 500)
    abs_interp = f(xf) 
    return abs_interp

abs22_numsym = adjust_resolution(sym_abs_2270)
abs30_numsym = adjust_resolution(sym_abs_3000)
abs22_numasym = adjust_resolution(asym_abs_2270)
abs30_numasym = adjust_resolution(asym_abs_3000)
##-------------------------##
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['text.latex.preamble'] = r"\boldmath"
##------------------------##
fig, ax = plt.subplots(2, 2, figsize=(13, 11),dpi=500)#, constrained_layout=True

# -------- Panel (0,0) — Symmetric 2270 Hz --------
ax[0, 0].plot(x_cen22sym, abs22_numsym, color='gray', lw=2.5, label=r"${\mathrm{{Num,sym}}}\,2270\,\mathrm{Hz}$")
ax[0, 0].plot(x_cen22sym, abs(p_cent22sym)/(2*4), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[0, 0].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[0, 0].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[0, 0].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[0, 0].get_xticklabels() + ax[0, 0].get_yticklabels():
    tick.set_fontweight('bold')
ax[0, 0].grid(alpha=0.3)
ax[0, 0].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# -------- Panel (0,1) — Symmetric 3000 Hz --------
ax[0, 1].plot(x_cen30sym, abs30_numsym, color='gray', lw=2.5,label=r"${\mathrm{{Num,sym}}}\,3000\,\mathrm{Hz}$")
ax[0, 1].plot(x_cen30sym, abs(p_cent30sym)/(2*7.26), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[0, 1].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[0, 1].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[0, 1].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[0, 1].get_xticklabels() + ax[0, 1].get_yticklabels():
    tick.set_fontweight('bold')
ax[0, 1].grid(alpha=0.3)
ax[0, 1].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# -------- Panel (1,0) — Asymmetric 2270 Hz --------
ax[1, 0].plot(x_cen22asym, abs22_numasym, lw=2.5, color='gray', label=r"${\mathrm{{Num,asym}}}\,2270\,\mathrm{Hz}$")
ax[1, 0].plot(x_cen22asym, abs(p_cent22asym)/(2*3.92), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[1, 0].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[1, 0].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[1, 0].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[1, 0].get_xticklabels() + ax[1, 0].get_yticklabels():
    tick.set_fontweight('bold')
ax[1, 0].grid(alpha=0.3)
ax[1, 0].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# -------- Panel (1,1) — Asymmetric 3000 Hz --------
ax[1, 1].plot(x_cen30asym, abs30_numasym, color='gray', lw=2.5, label=r"${\mathrm{{Num,asym}}}\,3000\,\mathrm{Hz}$")
ax[1, 1].plot(x_cen30asym, abs(p_cent30asym)/(2*9.26), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[1, 1].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[1, 1].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[1, 1].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[1, 1].get_xticklabels() + ax[1, 1].get_yticklabels():
    tick.set_fontweight('bold')
ax[1, 1].grid(alpha=0.3)
ax[1, 1].legend(prop={'weight':'bold', 'size':20}, frameon=False)


plt.show()
