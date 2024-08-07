import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.pyplot import figure
from readTable import readTable


Simu = "./ScatteringResult/ScatteringtoolsREF_He_Corrected.out"
Nesbet = "./data/electron/He/TCS/New_Nesbet.dat"
RPO_McEachran = "./data/electron/He/TCS/New_RPO_McEachran.dat"
BSR_McEachran = "./data/electron/He/TCS/New_BSR_McEachran.dat"
Wan = "./data/electron/He/TCS/New_Wan.dat"

data = readTable(Simu).data[:-1]
data = np.array(data,dtype=np.float64)

data1 = readTable(Nesbet).data[:-1]
data1 = np.array(data1,dtype=np.float64)

data2 = readTable(RPO_McEachran).data[:-1]
data2 = np.array(data2,dtype=np.float64)

data3 = readTable(BSR_McEachran).data[:-1]
data3 = np.array(data3,dtype=np.float64)

data4 = readTable(Wan).data[:-1]
data4 = np.array(data4,dtype=np.float64)

Scontrib = np.zeros(len(data[:,0]))
Pcontrib = np.zeros(len(data[:,0]))
Dcontrib = np.zeros(len(data[:,0]))
for i in range(len(data[:,0])) :
    Scontrib[i] = (2*np.pi/data[i,0])*np.sin(data[i,1])*np.sin(data[i,1])
    Pcontrib[i] = (6*np.pi/data[i,0])*np.sin(data[i,2])*np.sin(data[i,2])
    Dcontrib[i] = (10*np.pi/data[i,0])*np.sin(data[i,3])*np.sin(data[i,3])

    Scontrib[i] /= data[i,4]
    Pcontrib[i] /= data[i,4]
    Dcontrib[i] /= data[i,4]

Nesbet_Scontrib = np.zeros(len(data1[:,0]))
Nesbet_Pcontrib = np.zeros(len(data1[:,0]))
Nesbet_Dcontrib = np.zeros(len(data1[:,0]))
for i in range(len(data1[:,0])) :
    E = data1[i,0]*data1[i,0]*0.5
    Nesbet_Scontrib[i] = (2*np.pi/E)*np.sin(data1[i,1])*np.sin(data1[i,1])
    Nesbet_Pcontrib[i] = (6*np.pi/E)*np.sin(data1[i,2])*np.sin(data1[i,2])
    Nesbet_Dcontrib[i] = (10*np.pi/E)*np.sin(data1[i,3])*np.sin(data1[i,3])

    Nesbet_Scontrib[i] /= data1[i,4]
    Nesbet_Pcontrib[i] /= data1[i,4]
    Nesbet_Dcontrib[i] /= data1[i,4]


RPO_McEachran_Scontrib = np.zeros(len(data2[:,0]))
RPO_McEachran_Pcontrib = np.zeros(len(data2[:,0]))
RPO_McEachran_Dcontrib = np.zeros(len(data2[:,0]))
for i in range(len(data2[:,0])) :
    E = data2[i,0]*data2[i,0]*0.5
    RPO_McEachran_Scontrib[i] = (2*np.pi/E)*np.sin(data2[i,1])*np.sin(data2[i,1])
    RPO_McEachran_Pcontrib[i] = (6*np.pi/E)*np.sin(data2[i,2])*np.sin(data2[i,2])
    RPO_McEachran_Dcontrib[i] = (10*np.pi/E)*np.sin(data2[i,3])*np.sin(data2[i,3])

    RPO_McEachran_Scontrib[i] /= data2[i,4]
    RPO_McEachran_Pcontrib[i] /= data2[i,4]
    RPO_McEachran_Dcontrib[i] /= data2[i,4]


BSR_McEachran_Scontrib = np.zeros(len(data3[:,0]))
BSR_McEachran_Pcontrib = np.zeros(len(data3[:,0]))
BSR_McEachran_Dcontrib = np.zeros(len(data3[:,0]))
for i in range(len(data3[:,0])) :
    E = data3[i,0]*data3[i,0]*0.5
    BSR_McEachran_Scontrib[i] = (2*np.pi/E)*np.sin(data3[i,1])*np.sin(data3[i,1])
    BSR_McEachran_Pcontrib[i] = (6*np.pi/E)*np.sin(data3[i,2])*np.sin(data3[i,2])
    BSR_McEachran_Dcontrib[i] = (10*np.pi/E)*np.sin(data3[i,3])*np.sin(data3[i,3])

    BSR_McEachran_Scontrib[i] /= data3[i,4]
    BSR_McEachran_Pcontrib[i] /= data3[i,4]
    BSR_McEachran_Dcontrib[i] /= data3[i,4]


Wan_Scontrib = np.zeros(len(data4[:,0]))
Wan_Pcontrib = np.zeros(len(data4[:,0]))
Wan_Dcontrib = np.zeros(len(data4[:,0]))
for i in range(len(data4[:,0])) :
    E = data4[i,0]*data4[i,0]*0.5
    Wan_Scontrib[i] = (2*np.pi/E)*np.sin(data4[i,1])*np.sin(data4[i,1])
    Wan_Pcontrib[i] = (6*np.pi/E)*np.sin(data4[i,2])*np.sin(data4[i,2])
    Wan_Dcontrib[i] = (10*np.pi/E)*np.sin(data4[i,3])*np.sin(data4[i,3])

    Wan_Scontrib[i] /= data4[i,4]
    Wan_Pcontrib[i] /= data4[i,4]
    Wan_Dcontrib[i] /= data4[i,4]


fig = plt.figure(figsize = (10,10))
plt.plot(data[:,0]*27.21,Scontrib,'--k', markersize=10, linewidth = 6, label='Onde S')
plt.plot(data[:,0]*27.21,Pcontrib,':k', markersize=10, linewidth = 6, label='Onde P')
plt.plot(data[:,0]*27.21,Dcontrib,'-.k', markersize=10, linewidth = 6, label='Onde D')

plt.plot(data1[:,0]*data1[:,0]*13.605,Nesbet_Scontrib,'*r', markersize=10)
plt.plot(data1[:,0]*data1[:,0]*13.605,Nesbet_Pcontrib,'*r', markersize=10)
plt.plot(data1[:,0]*data1[:,0]*13.605,Nesbet_Dcontrib,'*r', markersize=10)

plt.plot(data4[:,0]*data4[:,0]*13.605,Wan_Scontrib,'+g', markersize=10)
plt.plot(data4[:,0]*data4[:,0]*13.605,Wan_Pcontrib,'+g', markersize=10)
plt.plot(data4[:,0]*data4[:,0]*13.605,Wan_Dcontrib,'+g', markersize=10)

plt.plot(data2[:,0]*data2[:,0]*13.605,RPO_McEachran_Scontrib,'xb', markersize=10)
plt.plot(data2[:,0]*data2[:,0]*13.605,RPO_McEachran_Pcontrib,'xb', markersize=10)
plt.plot(data2[:,0]*data2[:,0]*13.605,RPO_McEachran_Dcontrib,'xb', markersize=10)

plt.plot(data3[:,0]*data3[:,0]*13.605,BSR_McEachran_Scontrib,'oy', markersize=10)
plt.plot(data3[:,0]*data3[:,0]*13.605,BSR_McEachran_Pcontrib,'oy', markersize=10)
plt.plot(data3[:,0]*data3[:,0]*13.605,BSR_McEachran_Dcontrib,'oy', markersize=10)


plt.xticks(fontsize=25)
plt.yticks(fontsize=25)

plt.xlim([0,15])
plt.ylim([-0.05,1.0])

plt.legend(fontsize=25)
plt.savefig("Contrib_He.pdf")
plt.show()

