import numpy as np
from scipy.signal import resample
from scipy.interpolate import interp1d
import matplotlib as mpl
import matplotlib.pyplot as plt


sym_abs_pml_hard = np.loadtxt('./Num_center_pressure/sym_absPresdata_PML_Hard.txt', skiprows=8)
sym_abs_pml_pml = np.loadtxt('./Num_center_pressure/sym_absPresdata_PML_PML.txt', skiprows=8)

asym_abs_pml_hard=np.loadtxt('./Num_center_pressure/asym_absPresdata_PML_Hard.txt', skiprows=8)
asym_abs_pml_pml = np.loadtxt('./Num_center_pressure/asym_absPresdata_PML_PML.txt', skiprows=8)


def adjust_resolution(newdata):
    x_num = np.linspace(0,5, len(newdata))
    f = interp1d(x_num, newdata, kind='linear')
    xf = np.linspace(x_num.min(), x_num.max(), 500)
    abs_interp = f(xf) 
    return abs_interp


abs22_numsym = adjust_resolution(sym_abs_pml_hard[:,1])#sym_abs_2270
abs30_numsym = adjust_resolution(sym_abs_pml_hard[:,2])
abs22_numasym = adjust_resolution(asym_abs_pml_hard[:,1])
abs30_numasym = adjust_resolution(asym_abs_pml_hard[:,2])

abs22_sym = adjust_resolution(sym_abs_pml_pml[:,1])#sym_abs_2270
abs30_sym = adjust_resolution(sym_abs_pml_pml[:,2])
abs22_asym = adjust_resolution(asym_abs_pml_pml[:,1])
abs30_asym = adjust_resolution(asym_abs_pml_pml[:,2])

xval = np.linspace(np.min(sym_abs_pml_pml[:,0]), np.max(sym_abs_pml_pml[:,0]), 500)

##-------------------------##
plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['text.latex.preamble'] = r"\boldmath"

fig, ax = plt.subplots(2, 2, figsize=(16, 14),dpi=500)#, constrained_layout=True

# ----
ax[0, 0].plot(xval, abs22_numsym, color='gray', lw=2.5, label=r"$\mathrm{sym\;PML{-}Hard\;2270\,Hz}$")
ax[0, 0].plot(xval, abs22_sym, color='blue', lw=2.5, label=r"$\mathrm{sym\;PML{-}PML}$")
ax[0, 0].plot(xval, abs(p_cent22sym)/(4), lw=2.5, color='black', label=r"${\mathrm{{Exp,sym}}}$")
ax[0, 0].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[0, 0].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[0, 0].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[0, 0].get_xticklabels() + ax[0, 0].get_yticklabels():
    tick.set_fontweight('bold')
ax[0, 0].grid(alpha=0.3)
ax[0, 0].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# ------
ax[0, 1].plot(xval, abs30_numsym, color='gray', lw=2.5,label=r"$\mathrm{sym\;PML{-}Hard\;3000\,Hz}$")
ax[0, 1].plot(xval, abs30_sym, color='blue', lw=2.5, label=r"$\mathrm{sym\;PML{-}PML}$")
ax[0, 1].plot(xval, abs(p_cent30sym)/(7.26), lw=2.5, color='black', label=r"${\mathrm{{Exp,sym}}}$")
ax[0, 1].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[0, 1].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[0, 1].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[0, 1].get_xticklabels() + ax[0, 1].get_yticklabels():
    tick.set_fontweight('bold')
ax[0, 1].grid(alpha=0.3)
ax[0, 1].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# -----
ax[1, 0].plot(xval, abs22_numasym, lw=2.5, color='gray', label=r"$\mathrm{asym\;PML{-}Hard\;2270\,Hz}$")
ax[1, 0].plot(xval, abs22_asym, color='blue', lw=2.5, label=r"$\mathrm{asym\;PML{-}PML}$")
ax[1, 0].plot(xval, abs(p_cent22asym)/(3.92), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[1, 0].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[1, 0].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[1, 0].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[1, 0].get_xticklabels() + ax[1, 0].get_yticklabels():
    tick.set_fontweight('bold')
ax[1, 0].grid(alpha=0.3)
ax[1, 0].legend(prop={'weight':'bold', 'size':20}, frameon=False)

# ----
ax[1, 1].plot(xval, abs30_numasym, color='gray', lw=2.5, label=r"$\mathrm{asym\;PML{-}Hard\;3000\,Hz}$")
ax[1, 1].plot(xval, abs30_asym, color='blue', lw=2.5, label=r"$\mathrm{asym\;PML{-}PML}$")
ax[1, 1].plot(xval, abs(p_cent30asym)/(9.26), lw=2.5, color='black', label=r"${\mathrm{{Exp,asym}}}$")
ax[1, 1].set_xlabel(r"$\mathbf{x\,[m]}$", fontsize=24)
ax[1, 1].set_ylabel(r"$\mathbf{|P|\,[Pa]}$", fontsize=24)
ax[1, 1].tick_params(labelsize=22, width=1.5, length=6)
for tick in ax[1, 1].get_xticklabels() + ax[1, 1].get_yticklabels():
    tick.set_fontweight('bold')
ax[1, 1].grid(alpha=0.3)
ax[1, 1].legend(prop={'weight':'bold', 'size':20}, frameon=False)


plt.show()
