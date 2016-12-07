import glob
import commands
import matplotlib.pyplot as plt
import sys
from numpy import *

vals = []
maxvals = -1e100

for fl in glob.glob("samples.txt*"):
    lines = commands.getoutput("cut -c-100 " + fl).split('\n')
    
    lines = [item.split(",") for item in lines]
    for i in range(len(lines))[::-1]:
        try:
            lines[i] = float(lines[i][int(sys.argv[1])])
        except:
            del lines[i]
    vals.append(array(lines))
    maxvals = max(maxvals, max(lines))


for item in vals:
    plt.plot(maxvals - item)

plt.title("Max %f %f" % (maxvals, log10(abs(maxvals))))
plt.yscale('log')
plt.show()
    
