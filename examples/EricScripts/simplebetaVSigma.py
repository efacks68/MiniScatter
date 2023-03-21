import numpy as np
import matplotlib.pyplot as plt
import sys

#add sys.argv 1 or 2 for diff plot options

nParticles = [2.86e5,1.14e6,2.86e6,5.72e6]

if sys.argv[1] in {"p","P","%"}:
    mus = [3.575,3.571,3.569,3.569] # %OutsideTA
    sigmas = [0.0375,0.015,0.00908,7.45e-3] # %OutsideTA
elif sys.argv[1] in {"j","J"}:
    mus = [53.669,52.984,52.830,52.716]
    sigmas = [0.384,0.206,0.143,0.109]

fig,ax = plt.subplots()
ax.scatter(nParticles,mus,c='m',label=r"$\mu$",s=65)
ax2 = ax.twinx()
ax2.scatter(nParticles,sigmas,c='g',label=r"$\sigma$")

ax.set_xlabel(r"Simulation Macro-Particles")

lines, labels = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax2.legend(lines + lines2, labels + labels2, loc="upper right")

#plt.tight_layout()

if sys.argv[1] in {"p","P","%"}:
    plt.setp(ax,title=r"% Outside Target Area Convergence Study",ylabel=r"$\mu$ [%]")
    ax2.set_ylabel(r"$\sigma$ [%]")
    plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png",bbox_inches='tight')
    print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png")
elif sys.argv[1] in {"j","J"}:
    plt.setp(ax,title=r"Peak Current Density Convergence Study",ylabel=r"$\mu$ [$\mu$A/cm$^2$]")
    ax2.set_ylabel(r"$\sigma$ [$\mu$A/cm$^2$]")
    plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png",bbox_inches='tight')
    print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png")

plt.close()