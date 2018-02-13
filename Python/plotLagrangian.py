# -*- coding: utf-8 -*-

# ==================================================================================== #
# ================================ Program plot.py =================================== #
# ==================================================================================== #
# Author and corrections:                                                              #
#    - Eliot AYACHE    : Oct 2017                                                      #
#                                                                                      #
# Description:                                                                         #
#    - This program takes the output file of the lh1 code and creates a movie of the   #
#      evolution of the simulation                                                     #


# ==================================================================================== #
# ===================================== Imports ====================================== #
# ==================================================================================== #
import os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as ascii
# import time
# import progressbar as pgbar

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=10) 
plt.rcParams['savefig.dpi'] = 200

# ==================================================================================== #
# ==================================== Constants ===================================== #
# ==================================================================================== #
cells = False
machSwitch = True



# ==================================================================================== #
# ================================ SCRIPT STARTS HERE ================================ #
# ==================================================================================== #

inputadr = "../C_Bondi/lh1.out"


# Input file:
print("\nReading data...")
inputf = open(inputadr, "r")
data = ascii.read(inputf)
print("Done!\n")


# Renaming columns:
colnames = ['istep',
            't', 
            'ix', 
            'fm(ix)',
            'r(ix)/r(nx)',
            'u(ix)',
            'rho(ix)',
            'p(ix)',
            'eps(ix)',
            'w(ix)',
            'ix=1',
            'nx']
for i in range(len(data.columns)):
    data.rename_column(data.colnames[i],colnames[i])


# In the case of values normalised by c_s and r_s:
if machSwitch:
    cs = np.sqrt(4./3. * data['p(ix)'] / (data['rho(ix)']+1e-15))
    condiSound = (cs != 0.)
    data['u(ix)'][condiSound] /= cs[condiSound]


# Plotting
os.system("rm -f shock_*")
finished = False
startPos = 0

i = 0
nFigs = len(np.unique(data['istep']))
print("Plotting %s figures...\n----------------" %nFigs)


# Variables to plot:
sx = 'r(ix)/r(nx)'
sy = 'u(ix)'
sy2 = 'rho(ix)'
xcol = data[sx]
ycol = data[sy]
ycol2 = data[sy2]

while not(finished):

    # initializing plot:
    f, ax = plt.subplots(2, sharex = True)

    if i%10 == 0:
        print("steps %s - %s --> plotting..." %(i, str(i+10)))

    istep = data['istep'][startPos]
    condi = (data['istep'] == istep)
    x = xcol[condi]
    y = ycol[condi]
    y2 = ycol2[condi]

    if(cells):
        for j in range(len(x)):
            if j%20 == 0:
                plt.plot((x[j],x[j]),(0,y[j]), 'k-', linewidth = 0.5)

    ymax = np.max(ycol)
    ymin = np.min(ycol)
    ymax2 = np.max(ycol2)
    ymin2 = np.min(ycol2)
    xmax = np.max(xcol)
    xmin = np.min(xcol)
    startPos += len(x)



    # Figure 1:
    # --------
    ax[0].set_xlim(xmin,xmax)
    ax[0].set_ylim(ymin,ymax)
    ax[0].set_ylabel("$v_r/c_s$")
    ax[0].plot(x,y)
    


    # Figure 2:
    # --------
    ax[1].set_ylim(0,ymax2)
    # ax[1].set_ylim(ymin2,ymax2)
    ax[1].set_xlabel("$r$")
    ax[1].set_ylabel("$\\rho / \\rho_\mathrm{ext}$")
    ax[1].plot(x,y2,color = 'orange')

    for a in ax:
        a.label_outer()

    f.savefig("figs/shock_%07d.png" %istep)
    plt.close(f)

    if startPos == len(data):
        finished = True

    i += 1

print('\nExporting as movie...\n--------------')

os.chdir("figs")
os.system("rm -f movie.mp4")
os.system("ffmpeg -f image2 -pattern_type glob -i 'shock_*.png' -r 50 movie.mp4")
# os.system("ffmpeg -f image2 -pattern_type glob -i 'shock_*.png' -r 50 movie_%s_vs_%s.mp4" %(sx,sy))
os.chdir("..")











