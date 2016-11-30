import numpy as np
import matplotlib.pyplot as plt
import HqEvo as HQ
import Transform as Tf
import h5py

def dfdE(e, T, M):
	return np.exp(-e/T)*np.sqrt(e**2-M**2)*e
def dfdp(p, T, M):
	x = np.sqrt(p**2+M**2)/T
	return (x+1.)*np.exp(-x)

def corner(ds, ranges, bins=100):
	N = ds.shape[0]
	for i in range(N):
		for j in range(i+1):
			plt.subplot(N, N, i*N+j+1)
			if i==j:
				plt.hist(ds[i], bins=bins, range=[ranges[i,0], ranges[i,1]], histtype='step', normed=True)
				plt.xlim(ranges[i,0], ranges[i,1])
			else:
				plt.hist2d(ds[j], ds[i], range=[[ranges[j,0], ranges[j,1]],[ranges[i,0], ranges[i,1]]], bins=bins)
				plt.xlim(ranges[j,0], ranges[j,1])
				plt.ylim(ranges[i,0], ranges[i,1])

Sampler = HQ.HqEvo(inelastic=False)

def update(p1, v3cell, T, M):
	# Boost to cell frame and take down orientation of p1_cell
	p1_cell = Tf.boost4_By3(p1, v3cell)
	E1_cell = p1_cell[0]
	alpha1_cell = np.arctan2(p1_cell[2], p1_cell[1]) + np.pi/2.
	beta1_cell = np.arctan2(np.sqrt(p1_cell[1]**2+p1_cell[2]**2), p1_cell[3])
	gamma1_cell = 0.0

	# Sample channel in cell
	channel, dt = Sampler.sample_channel(E1_cell, T)

	if channel < 0:
		return p1
	
	# Sample E2_cell, s in cell
	E2_cell, s = Sampler.sample_initial(channel, E1_cell, T)
	
	# Imagine rotate p1_cell to align with z-direction, construct p2_cell_align
	p1z_cell_align = np.sqrt(E1_cell**2 - M**2)
	costheta2 = (M**2 + 2.*E1_cell*E2_cell - s)/2./p1z_cell_align/E2_cell
	sintheta2 = np.sqrt(1. - costheta2**2)
	phi2 = np.random.rand()*2.*np.pi
	cosphi2 = np.cos(phi2)
	sinphi2 = np.sin(phi2)
	p1_cell_align = np.array([E1_cell, 0.0, 0.0, p1z_cell_align])
	p2_cell_align = np.array([1.0, sintheta2*cosphi2, sintheta2*sinphi2, costheta2])*E2_cell


	# Center of mass frame of p1_cell_align and p2_cell_align, and take down orientation of p1_com
	Pcom = p1_cell_align + p2_cell_align
	v3com = Pcom[1:]/Pcom[0]
	p1_com = Tf.boost4_By3(p1_cell_align, v3com)
	alpha1_com = np.arctan2(p1_com[2], p1_com[1]) + np.pi/2.
	beta1_com = np.arctan2(np.sqrt(p1_com[1]**2+p1_com[2]**2), p1_com[3])
	gamma1_com = 0.0

	# Sample final state momentum in Com frame, with incoming paticles on z-axis
	final_ps = Sampler.sample_final(channel, s, T)
	p1_new_com_aligen = final_ps[0]
	
	# Rotate final states back to original Com frame (not z-axis aligened)
	p1_new_com = Tf.rotate_ByEuler(p1_new_com_aligen, -gamma1_com, -beta1_com, -alpha1_com)
	
	# boost by to cell frame z-aligned
	p1_new_cell_align = Tf.boost4_By3(p1_new_com, -v3com)
	
	# rotate back to original cell frame
	p1_new_cell = Tf.rotate_ByEuler(p1_new_cell_align, -gamma1_cell, -beta1_cell, -alpha1_cell)

	# boost back to lab frame
	p1_new = Tf.boost4_By3(p1_new_cell, -v3cell)
	
	# return updated momentum of heavy quark
	return p1_new
		
	
T = 0.4
v3cell = np.array([0.0, 0.0, 0.0])
M = 1.3
E = 10.0


f = h5py.File("HQ.hdf5", 'w')
ds = []
for i in range(10000):
	ds.append([])	
	p1 = np.array([E, 0.0, 0.0, (E**2-M**2)**0.5])
	for j in range(1000):
		p1 = update(p1, v3cell, T, M)
		ds[-1].append(p1)
ds = np.array(ds).T

ranges = np.array([[0, 11], [-5, 5], [-5, 5], [-5, 10]])
e = np.linspace(M, 11, 1000)
de = e[1]-e[0]
dpde = dfdE(e, T, M)
dpde = dpde/np.sum(dpde)/de
p = np.linspace(-5, 5, 1000)
dp = p[1]-p[0]
dpdp = dfdp(p, T, M)
dpdp = dpdp/np.sum(dpdp)/dp
for i in range(0, 1000, 50):
	plt.clf()	
	corner(ds[:,i,:], ranges)
	plt.subplot(4,4,1)
	plt.plot(e, dpde, 'r-')
	plt.subplot(4,4,6)
	plt.plot(p, dpdp, 'r-')
	plt.subplot(4,4,11)
	plt.plot(p, dpdp, 'r-')
	plt.subplot(4,4,16)
	plt.plot(p, dpdp, 'r-')
	plt.pause(0.01)
plt.show()















