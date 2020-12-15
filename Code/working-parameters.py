import math
import time
import numpy as np
import multiprocessing
import matplotlib.pyplot as plt
from multiprocessing.pool import Pool

#define function
#define main sequence
def SFR(M_star,z):
    s = A_11*((M_star*1.00/(10**11))**(1+beta))*((1+z)**alpha)
    return s

#define S(z)
def S(z, S_0_para, S_index_para):
    return S_0_para*((1+z)**S_index_para)

#define B(z)
def B(z, B_0_para, B_index_para):
    return B_0_para*((1+z)**B_index_para)

#define M_vir(z)
def M_vir(z, M_vir_0_para):
    return M_vir_0_para*math.exp(0-8.2*a_0/13.0+8.2/(13.0*(1+z)))

#define main function
def Func(para):
    name = multiprocessing.current_process().name
    print(name)
    #get an array of S
    S_values=[S(z, para[0], para[1]) for z in z_values]
    #get an array of B
    B_values=[B(z, para[2], para[3]) for z in z_values]
    #get an array of M_vir
    M_vir_values=[M_vir(z, para[4]) for z in z_values]
    #initial conditions
    SFR_0 = SFR(M_star_0,0)
    M_gas_0 = (SFR_0/para[0])*(10**9)    
    M_tot_0 = M_star_0+M_gas_0 
    #calculate delta_M_vir
    delta_M_vir_values = np.array([])
    flag = 0
    while flag <= num-2:
        delta_M_vir_values = np.append(delta_M_vir_values, M_vir_values[flag+1]-M_vir_values[flag])
        flag += 1
    #calculate delta_M_infall
    delta_M_infall_values = np.array([])
    flag = 0
    while flag <= num-2:
        result = 0.01*delta_M_vir_values[flag]
        delta_M_infall_values = np.append(delta_M_infall_values, result)
        flag +=1
    #calculate M_gas
    M_gas_values = np.array([M_gas_0])
    flag = 0
    while flag <= num-2:
        result1 = 0-(1-R+B_values[flag])*S_values[flag]*M_gas_values[flag]*(t_values[flag+1]-t_values[flag])+delta_M_infall_values[flag]
        result2 = M_gas_values[flag]+result1
        M_gas_values = np.append(M_gas_values, result2)
        flag +=1
    #calculate M_tot
    M_tot_values = np.array([M_tot_0])
    flag = 0
    while flag <= num-2:
        result1 = delta_M_infall_values[flag]-B_values[flag]*S_values[flag]*M_gas_values[flag]*(t_values[flag+1]-t_values[flag])
        result2 = M_tot_values[flag]+result1
        M_tot_values = np.append(M_tot_values,result2)
        flag += 1
    #calculate M_star
    M_star_values = np.array([])
    flag = 0
    while flag <= num-1:
        M_star_values = np.append(M_star_values, M_tot_values[flag]-M_gas_values[flag])
        flag += 1
    #calculate SFR 
    SFR_values = np.array([])
    flag = 0
    while flag <= num-1:
        SFR_values = np.append(SFR_values, S_values[flag]*M_gas_values[flag]/(10**9))
        flag += 1
    #judgement
    flag = 0
    while flag <= num-1:
        if M_star_values[flag]>0:
            flag += 1
        else:
            break
    judge_number = flag-1
    
    #Main sequence
    z_values_obs = np.array([0.0,0.25,0.5,0.75,1.0,1.25,1.5,1.75,2.0,2.25,2.5])

    M_star_values_obs = np.array([10**6,10**6.5,10**7,10**7.5,10**8,10**8.5,10**9,10**9.5,10**10,10**10.5,10**11,10**11.5])

    flag1 = 0
    flag2 = 0
    index = np.array([])
    z_values_ele = np.array([])
    while flag1 <= z_values_obs.shape[0]-1:
        while flag2 <= z_values.shape[0]-1:
            if z_values[flag2] >= z_values_obs[flag1]:
                index = np.append(index, flag2)
                z_values_ele = np.append(z_values_ele, z_values[flag2])
                break
            flag2 += 1
        flag1 += 1
    #plot
    plt.figure(figsize = [15,12])
    plt.subplot(2,2,1)
    plt.plot(t_values[0:judge_number],M_gas_values[0:judge_number],color="black",linestyle=":",label='${M_{gas}}$')
    plt.plot(t_values[0:judge_number],M_tot_values[0:judge_number],color="red",linestyle="-",label='${M_{tot}}$')
    plt.plot(t_values[0:judge_number],M_star_values[0:judge_number],color="blue",linestyle="--",label='${M_{star}}$')
    plt.xlabel("t[Gyr]",fontsize=10)
    plt.ylabel("${M}$/${M_{\odot}}$",fontsize=10)
    plt.legend(loc=1)
    
    plt.subplot(2,2,2)
    plt.plot(t_values[0:judge_number],SFR_values[0:judge_number])
    plt.xlabel("t[Gyr]",fontsize=10)
    plt.ylabel("${SFR}$${[M_{\odot}yr^{-1}]}$",fontsize=10)
    
    plt.subplot(2,2,3)
    flag = 0
    SFR_values_obs = np.array([])
    while flag <= z_values_obs.shape[0]-1:
        SFR_values_obs = np.array([SFR(M_star,z_values_obs[flag]) for M_star in M_star_values_obs])
        plt.loglog(M_star_values_obs,SFR_values_obs,label='z='+str(z_values_obs[flag]))
        flag += 1
        
    flag = 0
    M_star_values_model = np.array([])
    SFR_values_model = np.array([])
    while flag <= z_values_obs.shape[0]-1:
        M_star_values_model = np.append(M_star_values_model, M_star_values[int(index[flag])])
        SFR_values_model = np.append(SFR_values_model, SFR_values[int(index[flag])])
        flag += 1
    plt.loglog(M_star_values_model, SFR_values_model, '-o')
    plt.xlabel("${M_{*}}$$[{M_{\odot}}]$",fontsize=10)
    plt.ylabel("${SFR}$${[M_{\odot}yr^{-1}]}$",fontsize=10)
    
    plt.subplots_adjust(left=0.2, bottom=0.1, right=0.9, top=0.9)
    plt.legend(loc=2,labelspacing=0.02)
    plt.savefig('/huawei/osv1/chenyaoxin/workspace/Figure/working-fit2/pic{}.png'.format(para))
    plt.clf()
    return 0


if __name__ == '__main__':
    Time1 = time.time()
    y_z = 0.07
    f_b = 0.2
    R = 0.46
    Omega_L0 = 0.6911
    Omega_r0 = 0
    Omega_m0 = 0.3089
    h_0 = 0.6774
    A_11 = 3.24 #M_sun/yr
    beta = -0.35
    alpha = 3.45

    M_star_0 = 10**9
    a_0 = 1.0/1.0                                #1.0/(1+z), where z = 0
    S_0 = np.linspace(0.5, 1.5, 4)                #the SFE at z = 0
    S_index = np.linspace(-0.5, -2, 5)          #the parameter in S(z)
    B_0 = np.linspace(0.5, 2, 4)                #the mass loading factor at z = 0
    B_index = np.linspace(0.5, 2, 5)            #the parameter in B(z)
    M_vir_0 = np.linspace(10**10, 10**12, 7)    #the halo mass at z = 0

    #get an array of z and t
    Arr = dict(np.load('/huawei/osv1/chenyaoxin/workspace/Data/redshift-evolution_time10**(-3).npz'))
    #data from the code of upper block
    z_values = Arr['redshift']
    t_values = Arr['evolution_time']
    num = z_values.shape[0]

    #create a list of parameters
    flag1 = 0
    para_space = []
    while flag1 <= S_0.shape[0]-1:
        flag2 = 0
        while flag2 <= S_index.shape[0]-1:
            flag3 = 0
            while flag3 <= B_0.shape[0]-1:
                flag4 = 0
                while flag4 <= B_index.shape[0]-1:
                    flag5 = 0
                    while flag5 <= M_vir_0.shape[0]-1:
                        para_arr = np.array([S_0[flag1], S_index[flag2], B_0[flag3], B_index[flag4], M_vir_0[flag5]])
                        para_space.append(para_arr)
                        flag5 += 1
                    flag4 += 1
                flag3 += 1
            flag2 += 1
        flag1 += 1
    Time2 = time.time()
    print("Sequential execution time:", Time2-Time1)
    pool = Pool(50)
    hh = pool.map(Func, para_space)
    pool.close()
    pool.join()
    print(hh)
    Time3 = time.time()
    print("Parallel execution time:", Time3-Time2)
    print("The total time:", Time3-Time1)