import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2 as chiSq

x = np.linspace(0,4500,4500)
fs=18
mus = [250,1000,2000,3000,4000]
colors = ['r','orange','y','g','b']
fig,ax = plt.subplots()
ax.set_title(r"$\chi^2$ Distributions for Various $\lambda$",fontsize=fs)
ax.set_xlabel(r"$\lambda$",fontsize=fs)
ax.set_ylabel("Probability Density",fontsize=fs)
i=0
for mu in mus:
    ax.plot(x, chiSq.pdf(x, mu),'r--', label=r'$\chi^2$ PDF, $\lambda$='+str(mu)+r", $\sigma$="+"{:.2f}".format(np.sqrt(mu)),color=colors[i])
    i+=1
ax.legend()
ylim = ax.get_ylim()
ax.set_ylim(0.0001,ylim[1])
plt.tight_layout()
plt.savefig("chiSqDistribs.png",dpi=1000)