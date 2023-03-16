import numpy as np
import matplotlib.pyplot as plt

#add sys.argv 1 or 2 for diff plot options

#try to find plt.twinx use in jupyter notebook from openxal
nParticles = [2.86e5,1.14e6,2.86e6,5.72e6]
mus = [3.575,3.571,3.569,3.569] # %OutsideTA
sigmas = [0.0375,0.015,0.00908,7.45e-3] # %OutsideTA
#mus = [53.669,52.984,52.830,52.716]
#sigmas = [0.384,0.206,0.143,0.109]


fig,ax = plt.subplots()
ax.scatter(nParticles,sigmas,c='m',label=r"$\mu$",s=65)
ax2 = ax.twinx()
ax2.scatter(nParticles,mus,c='g',label=r"$\sigma$")
#plt.text(miniRatio[2],miniSizeX[2]+0.1,"3",fontsize=10)
ax.set_xlabel(r"Simulation Macro-Particles")
ax.set_ylabel(r"$\mu$ [%]")
ax2.set_ylabel(r"$\sigma$ [%]")
#ax.set_ylabel(r"$\mu$ [$\mu$A/cm$^2$]")
#ax2.set_ylabel(r"$\sigma$ [$\mu$A/cm$^2$]")
plt.title(r"% Outside Target Area Convergence Study")
#plt.title(r"Peak Current Density Convergence Study")
ax.legend()
ax2.legend()
plt.tight_layout()
plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png")
#plt.savefig("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png")
#plt.show()
plt.close()
print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/POutTAConvergence.png")
#print("/uio/hume/student-u52/ericdf/Documents/UiO/Forske/ESSProjects/PBWScattering/Pictures/PeakJConvergence.png")
