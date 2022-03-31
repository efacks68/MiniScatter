#plotFit.py
#Eric Fackelman
#29 March 2022

#This is to plot a general case

def singleplot(x,y,titl,labl,xlm,ylm,leg,savename):
  import matplotlib.pyplot as plt

  plt.plot(x,y)
  if labl != 0:
      plt.plot(x,y,label=labl)
  if xlm != 0:
    plt.xlim([-xlm,xlm])
  if ylm != 0:
    plt.ylim([-ylm,ylm])
  plt.title(titl)
  if leg:
    plt.legend()
  plt.savefig(savename+".png")
  plt.close()

def histplot(x,bins,titl,xlabl,ylabl,xlm,ylm,leg,savename,logy):
  import matplotlib.pyplot as plt

  plt.hist(x,bins)
  if logy:
    plt.hist(x,bins,log=True)
  if xlabl != 0:
    plt.xlabel(xlabl)
  if ylabl != 0:
    plt.ylabel(ylabl)
  if xlm != 0:
    plt.xlim([0,xlm])
  if ylm != 0:
    plt.ylim([0,ylm])
  plt.title(titl)
  if leg:
    plt.legend()
  print("Energy plotting ",savename)
  plt.savefig(savename+".png")
  plt.close()


