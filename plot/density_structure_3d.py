from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

fname = '/Users/Maya/Documents/GitHub/kernel-integration/src/for_3d_plotting.dat'

f = open(fname)

line = f.readline().split()
ind = int(line[0])

while ind!=0 :
    line = f.readline().split()
    x = [float(line[i]) for i in range(len(line))]
    line = f.readline().split()
    y = [float(line[i]) for i in range(len(line))]
    line = f.readline().split()
    z = [float(line[i]) for i in range(len(line))]
    
    tupleList = zip(x, y, z)
    vertices = []
        
    line = f.readline().split()
    while len(line)>=3 :
        vertices = vertices + [[int(line[i]) for i in range(len(line))]]
        line = f.readline().split()

    density = float(line[0])    
    line = f.readline().split()
    ind = int(line[0])
    
    poly3d=[]
    for ix in range(len(vertices)) :
        poly3d = poly3d + [[tupleList[vertices[ix][iy]] for iy in range(len(vertices[ix]))]]
        
    ax.set_xlim3d(-1, 1)
    ax.set_ylim3d(-1, 1)
    ax.set_zlim3d(-1, 1)
    ax.add_collection3d(Poly3DCollection(poly3d, facecolors='b', linewidths=1, alpha=density))
    ax.tick_params(axis='both', which='major', labelsize=15)
    ax.set_aspect('equal','box')

plt.savefig("/Users/Maya/Documents/GitHub/kernel-integration/plot/density-3d-50.pdf",format="pdf")

f.close()

plt.show()