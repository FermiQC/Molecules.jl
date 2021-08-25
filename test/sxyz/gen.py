import numpy as np
from scipy.optimize import minimize
name = "dodecahedron.xyz"
strang = ""

#n = 12
#inc = 2*np.pi / n
#for i in range(n):
#    theta = i*inc
#    x = np.cos(theta)
#    y = np.sin(theta)
#    z = 0
#    strang += f"H   {x:10.8f}   {y:10.8f}   {z:10.8f}\n"
#r2 = np.sqrt(2)
#lads = [[1,1,2],[1,-1,2],[-1,1,2],[-1,-1,2],[1,1,-1],[1,-1,-1],[-1,1,-1],[-1,-1,-1]]
#lads = [[1,1,2],[1,-1,2],[-1,1,2],[-1,-1,2],[0,r2,-1],[r2,0,-1],[-r2,0,-1],[0,-r2,-1]]
#for i in range(len(lads)):
#    strang += f"H   {lads[i][0]:10.8f}   {lads[i][1]:10.8f}   {lads[i][2]:10.8f}\n"



def ico(n):
    name = "icosahedron.xyz"
    strang = ""
    theta = [n, 180.0-n]
    phistart = np.asarray([0.0, 36.0])
    phi = np.zeros((2,5))
    
    for i in range(5):
        for j in range(2):
            phi[j,i] = phistart[j] + 72.0 * i
    
    mat = np.zeros((12,3))
    
    for i in range(2):
        for j in range(5):
            th = np.deg2rad(theta[i])
            ph = np.deg2rad(phi[i,j])
            x = np.sin(th)*np.cos(ph)
            y = np.sin(th)*np.sin(ph)
            z = np.cos(th)
            mat[5*i+j,:] = [x,y,z]
    
    mat[10,:] = [0.0,0.0,1.0]
    mat[11,:] = [0.0,0.0,-1.0]

    return mat

def dist(a,b):
    return np.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def fxna(n):
    mat = ico(n)
    return abs(dist(mat[0],mat[1]) - dist(mat[0],mat[5]))


n = minimize(fxna, [50.0]).x[0]
print(n)

mat = ico(n)

def cent(a,b,c):
    m = np.asarray([(a[0]+b[0]+c[0])/3, (a[1]+b[1]+c[1])/3, (a[2]+b[2]+c[2])/3])
    return m

l = []
for i in range(12):
    for j in range(12):
        if i < j:
            for k in range(12):
                if j < k:
                    if i != j != k:
                        a = cent(mat[i,:], mat[j,:], mat[k,:])
                        nor = np.linalg.norm(a)
                        if nor > 0.7:
                            l.append(a)
l = np.asarray(l)
print(l)
print(len(l))
for i in range(20):
    strang += f"H   {l[i,0]:10.8f}   {l[i,1]:10.8f}   {l[i,2]:10.8f}\n"

with open(name, 'w') as fn:
    fn.write(strang)

