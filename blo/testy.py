import numpy as np
import os
import resource

def mem(note):
    pid = os.getpid()
    with open(os.path.join('/proc', str(pid), 'status')) as f:
        lines = f.readlines()
    _vt = [l for l in lines if l.startswith("VmSize")][0]
    vt = mem_pretty(int(_vt.split()[1]))
    _vmax = [l for l in lines if l.startswith("VmPeak")][0]
    vmax = mem_pretty(int(_vmax.split()[1]))
    _vram = [l for l in lines if l.startswith("VmRSS")][0]
    vram = mem_pretty(int(_vram.split()[1]))
    _vswap = [l for l in lines if l.startswith("VmSwap")][0]
    vswap = mem_pretty(int(_vswap.split()[1]))

    vrammax = mem_pretty(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)

    print('## MEM -> total: %s (max: %s) | ram: %s (max: %s) | swap: %s @ %s' % (vt, vmax, vram, vrammax, vswap, note))

def mem_pretty(mem):
    denom = 1
    unit = 'kb'
    if (mem > 999999):
        denom = 1000000.0
        unit ='gb'
    elif (mem > 999):
        denom = 1000.0
        unit ='mb'
    return str(mem / denom) + unit

print(os.getpid())
mem('1')

# x = np.arange(225000000, dtype=np.float64).reshape(15000, 15000)
x = np.arange(225000000, dtype=np.int8).reshape(15000, 15000)

mem('2')

x = x.tolist()

mem('4')
# y = np.zeros((15000, 15000))
# mem('4')

# print x

# for i in xrange(2):
# 	x = np.delete(x, range(100), 0)
# 	mem('ok')

# mem('done')

raw_input("Press the <ENTER> key to continue...")