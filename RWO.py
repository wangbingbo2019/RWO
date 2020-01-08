# coding: utf-8
import networkx as nx
import numpy as np
import random
import os

# load network
# input file name format:
#     [NetName].txt
# input file format:
#     [gene1]\t[gene2]\r\n
def load_network(netpath,netName,s):
    os.chdir(netpath)
    a = open(netName, "r")
    G1 = nx.Graph()
    for i in a:
        n = i.strip().split("\t")
        G1.add_edge('_'.join([s,n[0]]), '_'.join([s,n[1]]))
    a.close()
    number = list(G1.nodes())
    numNodes = len(number)
    numEdges = G1.number_of_edges()
    return G1, number, numNodes, numEdges

# mapping
def mapping(G1, G2):
    matrix_mappingID = {}
    matrix_mappingName = {}
    matrix_mappingType = {}
    numAllGene = 0
    for x in G1.nodes():
        matrix_mappingID[x] = numAllGene
        matrix_mappingName[numAllGene] = x
        matrix_mappingType[numAllGene] = 'None'
        numAllGene = numAllGene + 1
    for x in G2.nodes():
        matrix_mappingID[x] = numAllGene
        matrix_mappingName[numAllGene] = x
        matrix_mappingType[numAllGene] = 'None'
        numAllGene = numAllGene + 1
    return matrix_mappingID, matrix_mappingName, matrix_mappingType, numAllGene

# load homology information
# input file name format:
#     [homoName].txt
# input file format:
#     [s1Gene]\t[s2Gene]\r\n
def load_homo_info(homoPath, homoName):
    tupleHomoInfo =[]
    os.chdir(homoPath)
    f = open(homoName,'r')
    for line in f.readlines():
        list = line.strip().split('\t')
        tupleHomoInfo.append((list[0],list[1]))
    f.close()
    return(tupleHomoInfo)

# initialization
def intro_initial_RW(G1, G2, count, matrix_mappingID):
    RW = np.zeros(shape=(count, count))
    for (x, y) in G1.edges():
        RW[matrix_mappingID[x]][matrix_mappingID[y]] = 1
        RW[matrix_mappingID[y]][matrix_mappingID[x]] = 1
    for (x, y) in G2.edges():
        RW[matrix_mappingID[x]][matrix_mappingID[y]] = 1
        RW[matrix_mappingID[y]][matrix_mappingID[x]] = 1
    return RW


# inter-initial RW
def inter_initial_RW(RW, tupleHomoInfo,
                     matrix_mappingID, matrix_mappingType):
    count = 0
    for (s1g,s2g) in tupleHomoInfo:
        if s1g in matrix_mappingID and s2g in matrix_mappingID:
            RW[matrix_mappingID[s1g]][matrix_mappingID[s2g]] = 1
            RW[matrix_mappingID[s2g]][matrix_mappingID[s1g]] = 1
            matrix_mappingType[matrix_mappingID[s2g]] = 'homo'
            matrix_mappingType[matrix_mappingID[s1g]] = 'homo'
            count = count + 1
    return RW, matrix_mappingType, count

# create probability transfer matrix
def Pmatrix(a1, a2, RW, matrix_mappingType, numAllGene, numNodeS1):
    Pr = np.zeros(shape=(numAllGene, numAllGene))
    for i in range(numAllGene):
        degree = sum(RW[i, :])
        homo = 0
        normal = 0
        if i < numNodeS1:
            across = sum(RW[i, numNodeS1:numAllGene])
            for x in range(numNodeS1):
                if matrix_mappingType[x] == 'homo' and RW[i][x] > 0:
                    homo = homo + 1
                if matrix_mappingType[x] == 'None' and RW[i][x] > 0:
                    normal = normal + 1

        if i >= numNodeS1:
            across = sum(RW[i, 0:numNodeS1])
            for x in range(numNodeS1, numAllGene):
                if matrix_mappingType[x] == 'homo' and RW[i][x] > 0:
                    homo = homo + 1
                if matrix_mappingType[x] == 'None' and RW[i][x] > 0:
                    normal = normal + 1
        alfa = a1
        beta = a2
        if degree > 0:
            if homo == 0:
                beta = 0
            if across == 0:
                alfa = 0
            if homo == 0 and normal == 0:
                alfa = 1
            for j in range(numAllGene):
                if i < numNodeS1 and j < numNodeS1:
                    if matrix_mappingType[j] == 'homo' and homo > 0:
                        Pr[i][j] = RW[i][j] * (1 - alfa) * beta / homo
                    if matrix_mappingType[j] == 'None' and normal > 0:
                        Pr[i][j] = RW[i][j] * (1 - alfa) * (1 - beta) / normal
                # break
                else:
                    if i >= numNodeS1 and j >= numNodeS1:
                        if matrix_mappingType[j] == 'homo' and homo > 0:
                            Pr[i][j] = RW[i][j] * (1 - alfa) * beta / homo
                        if matrix_mappingType[j] == 'None' and normal > 0:
                            Pr[i][j] = RW[i][j] * (1 - alfa) * (1 - beta) / normal
                       
                    else:
                        if across > 0:
                            Pr[i][j] = RW[i][j] * alfa / across
                           
    return Pr

# RWO
def output_RWO(rho, tau, RW, matrix_mappingType, numAllGene, numNodeS1,
               savePath, matrix_mappingName):
    P = Pmatrix(rho, tau, RW, matrix_mappingType, numAllGene, numNodeS1)
    Pt = np.ones(shape=(1, numAllGene))
    for i in range(200):
        Pt = np.dot(Pt, P)
    os.chdir(savePath)
    f = open('%.2f' % rho + 'result' + '%.2f' % tau + '.txt', 'w')
    f.write('\t'.join(['Gene', 'RWO']));f.write('\n')
    for i in range(numAllGene):
        f.write(matrix_mappingName[i])
        f.write('\t')
        f.write('%.8f' % Pt[0][i])
        f.write('\n')
    f.close()


def RWO(netPath,netNameList,sList,
        homoPath,homoName,
        savePath,rho=0.1, tau=0.4):
    # load network
    # input file name format:
    #     [NetName].txt
    # input file format:
    #     [gene1]\t[gene2]\r\n
    s1NetName = netNameList[0];s1 = sList[0]
    s2NetName = netNameList[1];s2 = sList[1]
    print('Loading network...')
    G1, numberS1, numNodeS1, numEdgeS1 \
        = load_network(netPath, s1NetName, s1)
    G2, numberS2, numNodeS2, numEdgeS2 \
        = load_network(netPath, s2NetName, s2)
    print('    Network information:')
    print('    ' + s1 + ' nodes: ' + str(numNodeS1) + ', edges: ' + str(numEdgeS1))
    print('    ' + s2 + ' nodes: ' + str(numNodeS2) + ', edges: ' + str(numEdgeS2))

    # mapping
    print('mapping...')
    matrix_mappingID, matrix_mappingName, matrix_mappingType, numAllGene \
        = mapping(G1, G2)

    # load homo information
    # input file name format:
    #     [homoName].txt
    # input file format:
    #     [s1Gene]\t[s2Gene]\r\n
    print('load homology information...')
    tupleHomoInfo = load_homo_info(homoPath, homoName)
    print('    Homology information:')
    print('    Number of information: '+str(len(tupleHomoInfo)))

    # initializaition
    print('Initializaition...')
    RW = intro_initial_RW(G1, G2, numAllGene, matrix_mappingID)
    RW, matrix_mappingType, count \
        = inter_initial_RW(RW, tupleHomoInfo, matrix_mappingID, matrix_mappingType)

    # RWO
    print('RWO...')
    # output_RWO(rho, tau, RW, matrix_mappingType, numAllGene, numNodeS1,
    #            savePath, matrix_mappingName)


if __name__ == '__main__':
    netPath = 'C:\Users\Administrator\Desktop\wang_2019_12_31'
    netNameList = ['yeastPPI.txt','HuRI.txt']
    sList = ['yeast', 'human']
    homoPath = 'C:\Users\Administrator\Desktop\wang_2019_12_31'
    homoName = 'yeast_human.txt'
    savePath = 'C:\Users\Administrator\Desktop\wang_2019_12_31'
    RWO(netPath, netNameList, sList,homoPath, homoName,savePath, rho=0.1, tau=0.4)


