import numpy as np
import matplotlib.pyplot as plt
import HqEvo, h5py, sys

def projector(p):
	norm2 = (p**2).sum()
	return np.diag([1,1,1]) - np.outer(p,p)/norm2 
	

M = 1.3
sampler = HqEvo.HqEvo(options={'transport': {'2->2': True, '2->3': True, '3->2': False},
						 'mass': M,
						 'Nf': 3,
						 'Kfactor': 1.,
						 'Tc': 0.154,
						 'mD': {'mD-model': 0, 'mTc': 4.0, 
								'slope': 0.25, 'curv':-0.8}}
				)

E = 100  # GeV
T = 0.2 # GeV
dt = 10 # GeV^-1
y0 = np.log(E/M + np.sqrt(1.+E*E/M/M))
sqrts = E + np.sqrt(E**2-M**2)
s = sqrts**2
pmax = (sqrts**2-M*M)/2./sqrts
p3 = []
p4 = []
k = []
"""
for i in range(10000):
	sampler.sample_final(2, s, T, dt, 0, 0)
	p3.append(sampler.FS[0])
	p4.append(sampler.FS[1])
	k.append(sampler.FS[2])
p3 = np.array(p3)
p4 = np.array(p4)
k = np.array(k)
f = h5py.File(sys.argv[1], 'w')
f.create_dataset('p3', data=p3)
f.create_dataset('p4', data=p4)
f.create_dataset('k', data=k)
f.close()
#"""

def plot(fname, c):
	f = h5py.File(fname, 'r')
	p3 = f['p3'].value
	p4 = f['p4'].value
	k = f['k'].value
	E4 = p4.T[0]
	k0 = k.T[0]
	kz = k.T[3]
	eta = 0.5*np.log((k0+kz)/(k0-kz))
	costheta4 = p4.T[3]/p4.T[0]
	kt = np.sqrt(k.T[1]**2 + k.T[2]**2)
	print(100.-np.mean(p3.T[0]))
	plt.subplot(2,2,1)
	plt.hist(np.log(k0/pmax), bins=100, range=[-10,0], histtype='step', color=c)
	plt.subplot(2,2,2)
	plt.hist(np.log(1.-E4/pmax), bins=50, range=[-15,0], histtype='step', color=c)
	plt.subplot(2,2,3)
	plt.hist(0.5*np.log((1.+costheta4)/(1.-costheta4)) + 0.5*np.log((E+M)/(E-M)), bins=50, range=[-10,10], histtype='step', color=c)

	plt.subplot(2,2,4)
	Tp4 = np.array([projector(ip[1:]) for ip in p4])
	direction = np.array([[0,0,-1] for ip in p4])
	vec1 = np.array([np.dot(PP, vv) for PP, vv in zip(Tp4, direction)])
	vec2 = np.array([np.dot(PP, vv[1:]) for PP, vv in zip(Tp4, k)])
	c4k = np.array([np.dot(v1, v2)/np.sqrt(np.dot(v1, v1)*np.dot(v2, v2)) for v1, v2 in zip(vec1, vec2)])

	plt.hist(np.arccos(c4k), range=[0,np.pi], bins=50, histtype='step', color=c)

plot(sys.argv[1], 'r') # 20
plot(sys.argv[2], 'g') # 50
plot(sys.argv[3], 'b') # 100
plt.show()



