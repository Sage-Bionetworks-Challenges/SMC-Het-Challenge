import numpy as np
import resource

print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
pred = np.arange(225000000).reshape(15000, 15000)
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
truth = np.arange(225000000).reshape(15000, 15000)
print(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
