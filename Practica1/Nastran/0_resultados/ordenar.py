import os
import numpy as np
from math import exp
import copy

# Para guardar los resultados en un csv
results_dir = './' # Ruta
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
def savedata(name,result):
    np.savetxt(results_dir+str(name)+'.csv', result, delimiter=',', fmt='%s')

def readdata(data):
    with open(data, 'r') as data:
        frecuencias = []
        modulos = []
        argumentos = []
        for line in data:
            if 'G' in line:
                p = line.split()
                frecuencias.append(float(p[1]))
                modulos.append(float(p[4]))
            elif 'T1' not in line:
                p = line.split()
                argumentos.append(float(p[1]))
        return [frecuencias, modulos, argumentos]

name = 'sinusoidal6.txt'
data_read = readdata(name)
savedata('sinusoidal6',data_read)
