"""Identification of the missing spectra in the STLM acquisition."""

import os

import matplotlib.pyplot as plt
import numpy as np

DATA_PATH = "Demo_Data/STL_8"
DATA_NAME = "280825-8-STL-TriggerData.txt"  # File with all the triggers data

# View of the transposed array
trigger_data = np.loadtxt(os.path.join(DATA_PATH, DATA_NAME), skiprows=1).T

time = trigger_data[0]
STMtoSpec = trigger_data[1]
SpectoSTM = trigger_data[2]

LL = 700  # width ot x axis
TR = 220

for k in range(int(max(time) // LL) + 1):
    # STM to spec TRiggers
    plt.plot(time, STMtoSpec, "b")
    plt.xlim(LL * k, LL * (k + 1))
    plt.show()

    # spec to STM triggers
    plt.plot(time, SpectoSTM, "r")
    plt.xlim(LL * k, LL * (k + 1))
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

diff = [trigSTM[i + 1] - trigSTM[i] for i in range(len(trigSTM) - 1)]
plt.plot(diff)
plt.ylim(0, 400)
plt.xlim(-50, 1050)
plt.show()

l = []
m = []
for i in range(len(trigSTM) - 1):
    if diff[i] < TR:
        l.append(diff[i])
        m.append(i)

m = np.array(m)

print(len(diff))

print(len(l))
print(l)
print(len(m))

# m is the list of missing and perverted spectra
print(m)

if m.size == 0:
    print("No missing spectra")

np.savetxt(os.path.join(DATA_PATH, DATA_NAME[0:-4] + "_missing_spectra.txt"), m)
