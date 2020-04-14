#Como rodar: python3 main.py
from funcoesTermosol import importa, geraSaida
from matrizes import matriz_restricoes, matriz_elemento
import numpy as np
from solver import solver_gauss
from math import *

numero_nos, matriz_nos, numero_elementos, matriz_conexoes, numero_forcas, vetor_forcas, numero_restricoes, vetor_restricoes = importa('grupo4.xlsx')

ke = np.zeros((numero_nos*2, numero_nos*2))

for conexoes in matriz_conexoes:
    no_1 = int(conexoes[0]) - 1
    grau_0 = no_1*2
    grau_1 = no_1*2 + 1
    
    no_2 = int(conexoes[1]) - 1
    grau_2 = no_2*2
    grau_3 = no_2*2 + 1
    
    x_inicio = matriz_nos[0][no_1]
    y_inicio = matriz_nos[1][no_1]
    
    x_fim = matriz_nos[0][no_2]
    y_fim = matriz_nos[1][no_2]
    
    E = conexoes[2]
    A = conexoes[3]

    matriz = matriz_elemento(x_inicio, x_fim, y_inicio, y_fim, E, A)
    ke[grau_0:grau_0+2, grau_0:grau_0+2] += matriz[:2, :2]
    ke[grau_2:grau_2+2, grau_0:grau_0+2] += matriz[:2, 2:]
    ke[grau_0:grau_0+2, grau_2:grau_2+2] += matriz[2:, :2]
    ke[grau_2:grau_2+2, grau_2:grau_2+2] += matriz[2:, 2:]

kg, vetor_forcas_restrito = matriz_restricoes(ke, vetor_forcas, vetor_restricoes)

vetor_deslocamento, _ = solver_gauss(kg, vetor_forcas_restrito)

for i in vetor_restricoes:
    vetor_deslocamento.insert(int(i), 0)

reacoes_apoio = np.array(ke).dot(np.array(vetor_deslocamento))
deformacoes = []
tensoes_internas = []
forcas_internas = []

for conexoes in matriz_conexoes:
    no_1 = int(conexoes[0]) - 1
    grau_0 = no_1*2
    grau_1 = no_1*2 + 1
    
    no_2 = int(conexoes[1]) - 1
    grau_2 = no_2*2
    grau_3 = no_2*2 + 1
    x_inicio = matriz_nos[0][no_1]
    y_inicio = matriz_nos[1][no_1]
    
    x_fim = matriz_nos[0][no_2]
    y_fim = matriz_nos[1][no_2]

    L = ((x_fim - x_inicio)**2 + (y_fim - y_inicio)**2)**0.5
    E = conexoes[2]
    A = conexoes[3]

    sen = (y_fim - y_inicio) / L
    cos = (x_fim - x_inicio) / L
    m = [-cos, -sen, cos, sen]
    graus = [grau_0, grau_1, grau_2, grau_3]
    vetor_desloc_local = [vetor_deslocamento[grau] for grau in graus]
    vetor_desloc_local = np.array(vetor_desloc_local)
    
    desloc = np.matmul(m, vetor_desloc_local)/L
    tensao = E*desloc
    forca = tensao*A

    deformacoes.append(desloc)
    tensoes_internas.append(tensao)
    forcas_internas.append(forca)

AREA = 0.00015
AREA2 = 0.0003

COMPRIMENTO = sqrt(0.08**2+0.04**2)

massa = (0.72*AREA + 0.08*AREA2 + 10*(COMPRIMENTO*AREA))*848

print("Massa:", massa*1000, " (g)")


delete_indexes = []
for i in range(len(reacoes_apoio)):
    if i not in vetor_restricoes:
        delete_indexes.append(i)

for i in range(len(vetor_deslocamento)):
    if abs(vetor_deslocamento[i]) > 0.02:
        print("Nó deslocado:", floor(i/2) + 1)

for i in range(len(tensoes_internas)):
    if abs(int(tensoes_internas[i])) > 18e6:
        print("Elemento acima da tensão:", i + 1)

geraSaida("grupo4", np.array(reacoes_apoio, float), np.array(vetor_deslocamento, float), np.array(deformacoes, float), np.array(forcas_internas, float), np.array(tensoes_internas, float))