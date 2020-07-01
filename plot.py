import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

def time_shape(t, tpulse):
    sigma = 4*np.log(2)/tpulse**2
    return np.sqrt(sigma/np.pi) * np.exp(-sigma*(t-3*tpulse)**2)

ttmpath = 'Saved\\TTM\\Au\\230.00_J_m_2_200.00_fs_400.00_nm'

cdata = np.loadtxt(os.path.abspath(os.path.join(ttmpath, 'conservation_10.dat')))
plt.figure(figsize=(8,6))
ax = plt.subplot()
l1 = ax.plot(cdata[:,0], cdata[:,1], label='F$_{abs}$')
l2 = ax.plot(cdata[:,0], cdata[:,2], label='E$_{el}$')
l3 = ax.plot(cdata[:,0], cdata[:,3], label='E$_{l}$')
l4 = ax.plot(cdata[:,0], cdata[:,2]+cdata[:,3], lw=2, ls='--', label='E$_{l}$+E$_{el}$')
ax1 = ax.twinx()
l5 = ax1.plot(cdata[:,0], cdata[:,4]*100, lw=0.5, ls=':', label='$\delta$E', color='black')
ax.set_xlabel('t (ps)')
ax.set_ylabel('F (J/m$^2$)')
ax1.set_ylabel('Error (%)')
# added these lines
lns = l1+l2+l3+l4+l5
labs = [l.get_label() for l in lns]
ax.legend(lns, labs, loc=0)
plt.savefig('conservation.png', dpi=150)
plt.close()

# My data
ttm_data = np.loadtxt(os.path.abspath(os.path.join(ttmpath, 'profiles_10.dat')))
time = np.unique(ttm_data[:, 0])
num_uniq = time.shape[0]
ind = np.where(ttm_data[:, 0] == time[1])[0][0]
time = ttm_data[:,0].reshape((num_uniq, ind)) #ps
coord = ttm_data[:,1].reshape((num_uniq, ind)) #nm
Te = ttm_data[:,2].reshape((num_uniq, ind))
Tl = ttm_data[:,3].reshape((num_uniq, ind))
# Hohlfeld data
Tedat = pd.read_csv('Hohlfeld_Te.csv', sep=' ', header=None, decimal=',')
Tldat = pd.read_csv('Hohlfeld_Tl.csv', sep=' ', header=None, decimal=',')
Tedat = Tedat.to_numpy()
Tldat = Tldat.to_numpy()

plt.figure(figsize=(8,6))
plt.suptitle("Temperature dynamics on the surface",fontsize = 16)
plt.title("Stack [100nm Au|5000 nm Si]")
ax = plt.subplot()
ax.set_xlabel("t (ps)",fontsize = 16)
ax.set_ylabel("T (K)",fontsize = 16)
ax.plot(Tedat[:,0], Tedat[:,1],label = f"T$_e$, Hohlfeld")
ax.plot(Tldat[:,0], Tldat[:,1],label = f"T$_l$, Hohlfeld")
ax.plot(time[:,0], Te[:,0],label='T$_e$, my code')
ax.plot(time[:,0], Tl[:,0],label='T$_l$, my code')
[y1, y2] = ax.get_ylim()
plt.plot(time[:,0], time_shape(time[:,0], 0.2)/np.amax(time_shape(time[:,0], 0.2))*y2/1.2, ls=':', label='pulse shape')
plt.tick_params(axis='x', labelsize=14)
plt.tick_params(axis='y', labelsize=14)
plt.legend(loc = "upper right",fontsize = 16)
plt.xlim(0, 60)
plt.grid()
plt.savefig('TTM.png')
plt.close()
