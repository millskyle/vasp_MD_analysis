import matplotlib.pyplot as plt
import numpy as np


ntimesteps = 800;

fig = plt.figure()
ax = fig.add_subplot(111)
plt.ion()
print "Plotting frames"

xlim=ylim=[-5,20]


for t in xrange(ntimesteps):
   data = np.loadtxt("animated_wrapped"+str(t)+".data")
   ax.cla()
   ax.set_ylim(ylim)
   ax.set_xlim(xlim)
   print data[:,0]
   plt.scatter(data[:,0],data[:,1])
   plt.draw()
   print t

















