## build alternative ancestral sequences

from fasta2seq import *
from itertools import product

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

import json,random

fp = open('spr_all.prob.json','r')
prob = json.load(fp)
numseq = 1000
#mutations = 30

nodes = pd.read_csv('anc.node.txt')
nodes = dict(zip(nodes['node'],nodes['name']))

#nodes = {"471":"Drosophila"}
#nodes = {"471":"Drosophila"}
#471,Drosophila
#435,Diptera
#369,Insecta
#312,Arthropoda
#255,Metazoa
#nodes = {"435":"Diptera"}

for node in nodes:
    sequence = {}
    name = nodes[node]
    node = str(node)
    print(node,name)

    ## all of the alternative residues ##
    alternative_residues = []
    for res_idx in prob[node]:
        alt_res = [res for res in prob[node][res_idx]]
        alternative_residues.append(alt_res)

    ## random sample 100 sequences ##
    N = 0
    index_list = []
    while N<numseq:
        index = []

        residues_with_ALT = []
        for res_idx in prob[node]:
            if len(prob[node][res_idx]) > 1:
                #print(prob[node][res_idx])
                residues_with_ALT.append(res_idx)

        print(len(residues_with_ALT))
        #print(residues_with_ALT)

        mutations = int(len(residues_with_ALT)/2)
        residues_to_mutate = random.choices(residues_with_ALT,k=mutations)
        print(residues_to_mutate)

        #print([idx for idx in prob[node]])
        for res_idx in prob[node]:
            alt_res = [res for res in prob[node][res_idx]]
            alt_idx = [i for i in range(len(alt_res))]
            n = len(alt_res)
            #if n>1:
            #    idx = random.choice(alt_idx)
            #    index.append(idx)
            if res_idx in residues_to_mutate:
                index_pool = [i for i in range(n) if i>0]
                idx = random.choice(index_pool)
                print(res_idx,alt_res,alt_idx,idx)
                index.append(idx)
            else:
                index.append(0)
        if not all([i==0 for i in index]):
            index_list.append(index)
            N += 1

    nres = len(alternative_residues)
    #print(nres,len(index_list))
    for n in range(numseq):
        seq = [alternative_residues[i][index_list[n][i]] for i in range(nres)]
        #seq = [s for s in seq if s!='-']
        seq = ''.join(seq)
        sequence['%s_Alternative_%d'%(name,n)] = seq
#        print(ni,len(seq),seq)
    seq2fasta(sequence,'SPR_Anc_alternative_%s.fasta'%name)
