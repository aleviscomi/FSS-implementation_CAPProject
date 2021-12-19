import numpy as np;
import numpy.random;
import array;

def dsload(filename, dtype = 'f'):
	with open(filename,'rb') as datafile:
		size = array.array('I');
		size.fromfile(datafile, 2);
		data = array.array(dtype);
		data.fromfile(datafile, size[0]*size[1]);   
	data = np.array(data);
	if(size[0]>1):
		return data.reshape(size[1],size[0]);
	return data;

def dssave(ds, filename, dtype = 'f'):
	if(ds.ndim == 2):
		n, d = ds.shape;
	else:
		n, d = len(ds), 1;
	with open(filename,'wb') as datafile:
		array.array('I', [d, n]).tofile(datafile);
		array.array(dtype, ds.reshape(np.prod(ds.shape))).tofile(datafile);

d = 8; # 8 125 256
nx = 64; # 64 735 768
iters = np.min([int(250*np.round(75*np.sqrt(d)/250,0)), int((2**31-1)/((nx*(d+1)+d)*8)/4)*4]); # 250 750 1250
scale = np.power(10,-1*np.log10(d));

if d == 8:
        si = 0.1;
        sv = 0.0125;
        w = 1.25;

if d == 125:
        si = 0.006;
        sv = 0.0008;
        w = 0.08;

if d == 256:
        si = 0.003;
        sv = 0.0004;
        w = 0.04;
        
print("PARAMETERS:");
print("\tnx:", nx, "d:", d, "iters:", iters);

print("FILE SIZE f:",(iters*(nx*(d+1))+d)*4/1024**3,"G")
print("FILE SIZE d:",(iters*(nx*(d+1))+d)*8/1024**3,"G")
if (iters*(nx*(d+1))+d)*8 > 2**31-1:
    print("Troppo grande!");
    quit();

print("Genero il file rand.");
r = np.random.rand(nx*d + iters*(nx*(d+1))+d);
dssave(r,'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/rand32_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='f');
dssave(r,'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/rand64_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='d');
del r;

print("Genero il file coeff.");
c=-1*np.random.rand(d)*scale;
dssave(c, 'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/coeff32_'+str(d)+'.ds2', dtype='f');
dssave(c, 'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/coeff64_'+str(d)+'.ds2', dtype='d');

a = -1*scale;
b = 1*scale;
print("Genero il file x init.");
ds = (np.random.rand(nx*d)*(b-a) + a).reshape((nx,d));
dssave(ds, 'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/x32_'+str(d)+'_'+str(nx)+'.ds2', dtype='f');
dssave(ds, 'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/x64_'+str(d)+'_'+str(nx)+'.ds2', dtype='d');
del ds;

def f(x): return np.exp(np.sum(x**2,axis=x.ndim-1))+ np.sum(x**2-c*x,axis=x.ndim-1);

import scipy as sc;
import scipy.optimize;
xok = sc.optimize.minimize(f,a/2*np.ones(d)).x;

print("\nLa soluzione è:");
for i in range(0, d):
        print("xsol[",i,"] =", xok[i]);
print("\nf(xok) =", f(xok));
print('-c data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/coeff32_'+str(d)+'.ds2 -r data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/rand32_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2 -x data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/x32_'+str(d)+'_'+str(nx)+'.ds2 -np '+str(nx)+' -si '+str(si)+' -sv '+str(sv)+' -w '+str(w)+' -it '+str(iters)+' -d');

dssave(xok.reshape(1,d),'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/sol32_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='f');
dssave(xok.reshape(1,d),'cmake-build-debug/data/data'+str(d)+'_'+str(nx)+'_'+str(iters)+'/sol64_'+str(d)+'_'+str(nx)+'_'+str(iters)+'.ds2',dtype='d');
