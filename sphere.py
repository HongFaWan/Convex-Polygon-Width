from mpl_toolkits.mplot3d import Axes3D, axes3d
import matplotlib.pyplot as plt
import numpy as np

# center and radius
center = [0, 0, 0]
radius = 1

# data
u = np.linspace(0, 2 * np.pi, 200)
v = np.linspace(0, np.pi, 200)
x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]


with open('Japan_BL.dat','r') as f:
    data=f.readlines()
L=np.zeros((len(data),2))
print(len(data))
for i in range(len(data)):
    L[i,1],L[i,0]=data[i].split()
min0=np.min(L[:,0])
max0=np.max(L[:,0])
min1=np.min(L[:,1])
max1=np.max(L[:,1])
for i in range(len(data)):
    L[i,0]=1*(L[i,0]-min0)+min0
    L[i,1]=1*(L[i,1]-min1)+min1

# with open('coor2_s.txt','w') as f:
#     for i in range(len(data)):
#         f.write("{0}  {1}\n".format(L[i,0],L[i,1]))
L[:,0]*=np.pi/180.0
L[:,1]*=np.pi/180.0
X=(radius * np.cos(L[:,0])*np.cos(L[:,1]) + center[0])
Y = radius * np.sin(L[:,0])* np.cos(L[:,1]) + center[1]
Z = radius *np.sin(L[:,1]) + center[2]

X=np.reshape(X,(1,-1))
Y=np.reshape(Y,(1,-1))
Z=np.reshape(Z,(1,-1))
XYZ=np.vstack((X,Y,Z))
# print(Z)
# plot
fig = plt.figure()
# ax = fig.add_subplot(121, projection='3d')
#
# # surface plot rstride 值越大，图像越粗糙
# ax.plot_surface(x, y, z, rstride=4, cstride=4, color='b')

# wire frame
ax = fig.add_subplot(111, projection='3d')
ax.plot_wireframe(x, y, z, rstride=10, cstride=10)


# ax = fig.add_subplot(111, projection='3d')
# ax.plot3d(X, Y, Z)
ax.scatter(X,Y,Z,c='r',marker='.',s=10,linewidth=0,alpha=1,cmap='spectral')

# R=np.matrix([[0.058677,0.985515,-0.159114],
#  [-0.369626,0.169508,0.913588],
#   [0.927326,0.00520599,0.374218]])

# R=np.matrix([[0.545042,0.741249,-0.391764],
# [-0.305358,0.610683,0.730632],
#  [0.780824,-0.278597,0.559194]])



R=np.matrix([[0.13437,0.777729,-0.614071],
[-0.580652,0.563953,0.587197],
 [0.802987,0.27766,0.527368]])
#
# xi1,yi1=124.394,45.4394
# xi2,yi2=128.896,43.5401
# xj,yj=125.604,53.0761
# x1=np.linspace(0,2*np.pi,50)
# l=127.700562149253*np.pi/180.0
# fai=26.1181258772270*np.pi/180.0
# x2 =(radius * np.cos(l)*np.cos(fai) + center[0])
# y2 = radius * np.sin(l)* np.cos(fai) + center[1]
# z2 = radius *np.sin(fai) + center[2]
# xyz=np.matrix([
#     [x2],
#     [y2],
#     [z2]
# ])
# print(xyz)
# xyz=R@xyz
XYZ=R@XYZ
ax.scatter(XYZ[0,:],XYZ[1,:],XYZ[2,:],c='g',marker='.',s=10,linewidth=0,alpha=1,cmap='spectral')
L2=np.zeros((len(data),2))
for i in range(len(data)):
    L2[i,0]=np.arctan2(XYZ[1,i],XYZ[0,i])
    L2[i, 1]=np.arcsin(XYZ[2,i])


maxfai=np.max( L2[:,1])
minfai=np.min( L2[:,1])
print(maxfai,minfai)

# l2=np.arctan2(xyz[1,0],xyz[0,0])
# fai=np.arcsin(xyz[2,0])
# print(l2*180.0/np.pi,fai*180.0/np.pi)

# print(xyz)
# show
plt.show()