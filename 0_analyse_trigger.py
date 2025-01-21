import os
import numpy as np
import matplotlib.pyplot as plt

data_path = 'C:/Users/cfo/Documents/Data_Analysis/STL_Acquisition/STL_Data/2024-09-03 #SG19 ALE W UHV LN/STL_16'
dataname = "030924-16-STL-TriggerData.txt"  # File with all the triggers data

# View of the transposed array
trigger_data = np.loadtxt(os.path.join(data_path, dataname), skiprows=1).T

time = trigger_data[0]
STMtoSpec = trigger_data[1]
SpectoSTM = trigger_data[2]

ll = 700  # width ot x axis
tr = 220

for k in range(int(max(time)//ll) + 1):
    # STM to spec triggers
    plt.plot(time, STMtoSpec, 'b')
    plt.xlim(ll*k, ll*(k+1))
    plt.show()

    # spec to STM triggers
    plt.plot(time, SpectoSTM, 'r')
    plt.xlim(ll*k, ll*(k+1))
    plt.show()

trigSTM = []
i = 0
while i < len(time):
    if STMtoSpec[i] > 2:
        trigSTM.append(i)
        i += 5
    else:
        i += 1

trigSpectro = []
i = 0
while i < len(time):
    if SpectoSTM[i] > 2:
        trigSpectro.append(i)
        i += 5
    else:
        i += 1

trigSTM = np.array(trigSTM)
trigSpectro = np.array(trigSpectro)

print(len(trigSTM))
print(len(trigSpectro))

plt.plot(trigSTM)
plt.plot(trigSpectro)
plt.show()

diff = [trigSTM[i+1] - trigSTM[i] for i in range(len(trigSTM) - 1)]
plt.plot(diff)
plt.ylim(0, 400)
plt.xlim(-50, 1050)
plt.show()

l = []
m = []
for i in range(len(trigSTM) - 1):
    if diff[i] < tr:
        l.append(diff[i])
        m.append(i)

m = np.array(m)

print(len(diff))

print(len(l))
print(l)
print(len(m))

# m is the list of missing and perverted spectra
print(m)
np.savetxt(os.path.join(data_path, dataname.strip('.txt') + '_missing_spectra.txt'), m)
