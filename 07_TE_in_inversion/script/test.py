import numpy as np 

L1 = [1,5,85,9,3,6,9]
L2 = [51,57,89,36,54,74,15,1525,78]

L3 = L1 + L2

print(np.mean(L1))
print(np.mean(L2))
print((np.mean(L1)+np.mean(L2))/2)
print(np.mean(L3))

