import numpy as np
import matplotlib.pyplot as plt
import sys

#add sys.argv 1 or 2 for diff plot options

nParticles = [2.86e5,1.14e6,2.86e6,5.72e6]

if sys.argv[1] in {"p","P","%"}:
    mus = [3.575,3.571,3.569,3.569] # %OutsideTA
    sigmas = [0.036,0.015,0.00944,7.59e-3] # %OutsideTA
elif sys.argv[1] in {"j","J"}:
    mus = [53.672,52.984,52.831,52.718]
    sigmas = [0.388,0.213,0.149,0.115]

fs=16
fig,ax = plt.subplots()
ax.scatter(nParticles,mus,c='m',label=r"$\mu$",s=100)
ax2 = ax.twinx()
ax2.scatter(nParticles,sigmas,c='mediumseagreen',label=r"$\sigma$")

ax.set_xlabel(r"Simulation Macro-Particles",fontsize=fs)

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc="upper right")

#plt.tight_layout()

if sys.argv[1] in {"p","P","%"}:
    ax.set_title("% Outside Target Area \nConvergence Study",fontsize=fs+2)
    ax.set_ylabel(r"$\mu$ [%]",fontsize=fs)
    ax2.set_ylabel(r"$\sigma$ [%]",fontsize=fs)
    plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png",bbox_inches='tight',dpi=500)
    print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png")
elif sys.argv[1] in {"j","J"}:
    ax.set_title(r"Peak $\langle \bf{J} \rangle$"+"\nConvergence Study",fontsize=fs+2)
    ax.set_ylabel(r"$\mu$ [$\mu$A/cm$^2$]",fontsize=fs)
    ax2.set_ylabel(r"$\sigma$ [$\mu$A/cm$^2$]",fontsize=fs)
    plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png",bbox_inches='tight',dpi=500)
    print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png")

plt.close()