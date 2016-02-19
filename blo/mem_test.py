import numpy as np
import resource
import os

def mem():
	pid = os.getpid()
	with open(os.path.join('/proc', str(pid), 'status')) as f:
		lines = f.readlines()
	_vmsize = [l for l in lines if l.startswith("VmSize")][0]
	vmsize = int(_vmsize.split()[1])

	physical = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss

	print('v: %s | p: %s' % (vmsize, physical))

# print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
mem()
pred = np.arange(225000000).reshape(15000, 15000)
# print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
mem()
truth = np.arange(225000000).reshape(15000, 15000)
# print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
mem()