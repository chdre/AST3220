import matplotlib.pyplot as plt
import numpy as np

a = np.linspace(0,20,10)
func = lambda av: 1./(np.exp(a)+1)

for i in range(10):
    print np.exp(a)
    plt.plot(a,func(a[i]-10.))
    plt.show()
