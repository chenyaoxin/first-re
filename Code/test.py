import math
import matplotlib.pyplot as plt
import numpy as np

Arr = dict(np.load('/huawei/osv1/chenyaoxin/workspace/Data/redshift-evolution_time10**(-5).npz'))
t_values = Arr['evolution_time']
flag = 0
delta_t_true_values = np.array([])
while flag <= (t_values.shape[0]-2):
    delta_t_true_values = np.append(delta_t_true_values, abs(t_values[flag]-t_values[flag+1]))
    flag += 1

plt.hist(delta_t_true_values, bins = 80)

plt.savefig('/huawei/osv1/chenyaoxin/workspace/Figure/working-test/PDF') 