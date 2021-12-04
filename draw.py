import matplotlib.pyplot as plt
import numpy as np
import cv2

# with open('CN-border-La.gmt','r') as f:
#     data=f.readlines()

with open('Japan_BL.dat','r') as f:
    data=f.readlines()

# for i in range(1,127026):
#     if(data[i]=='>\n'):
#         print(i)

# x=np.zeros((5784,2))
# for i in range(1,5785):
#     x[i-1,0],x[i-1,1]=data[i].split()
#
# x=np.zeros((5677,2))
# for i in range(5786,11463):
#     x[i-5786,0],x[i-5786,1]=data[i].split()

x=np.zeros((922,2))
for i in range(922):
    x[i,1],x[i,0]=data[i].split()

#######################################
# y=np.zeros((34,2))
# with open('out.txt','r') as f:
#     data=f.readlines()
#
#
# for i in range(0,34):
#     y[i,0],y[i,1]=data[i].split()
#######################################

# y=np.zeros((29,2))
# with open('out2.txt','r') as f:
#     data=f.readlines()
#
# for i in range(0,29):
#     y[i,0],y[i,1]=data[i].split()

########################################

y=np.zeros((16,2))
with open('Japan_BL_out.txt','r') as f:
    data=f.readlines()

for i in range(0,16):
    y[i,0],y[i,1]=data[i].split()



############################################

# with open('coor.txt','w') as w:
#     for i in range(0, 5784):
#         w.write('{0}\t{1}\n'.format(x[i , 0], x[i , 1]))

fig = plt.figure()
# 将画图窗口分成1行1列，选择第一块区域作子图
ax1 = fig.add_subplot(1, 1, 1)

ax1.scatter(x[:,0], x[:,1], c='k', marker='.')

ax1.scatter(y[:,0], y[:,1], c='#FF0000', marker='.',s=300)

x_max=max(y[:,0])
x_min=min(y[:,0])
# xi1,yi1=124.394,45.4394
# xi2,yi2=128.896,43.5401
# xj,yj=125.604,53.0761

xi1,yi1=141.651,45.44
xi2,yi2=129.335,34.6339
xj,yj=140.335,35.1314

ax1.scatter([xi1,xj], [yi1,yj], c='#008000', marker='*',s=300)


ax1.scatter([xi2], [yi2], c='#FFFF00', marker='*',s=300)

ax1.scatter(135.836,35.9583, c='r', marker='*',s=1000)


x1=np.linspace(x_min,x_max,50)

a=(yi2-yi1)/(xi2-xi1)

b1=yi1-a*xi1

b2=yj-a*xj

y1=a*x1+b1

y2=a*x1+b2

ax1.plot(x1, y1, color='green', linewidth=1.0, linestyle='--')

ax1.plot(x1, y2, color='green', linewidth=1.0, linestyle='--')

print(abs(b1-b2)/np.sqrt(1+a*a))



plt.show()
