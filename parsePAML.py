## parse PAML output into a dictionary of ancestral states

from ete3 import Tree
from fasta2seq import *

import json
import numpy as np
import pandas as pd
import subprocess,shutil

def chop_seq(inp,out):
    sequence = fasta2seq(inp)
    for key in sequence:
        seq_ref = sequence[key]
        break
    seq_len = len(seq_ref)
    sequence_new = dict()
    for key in sequence:
        seq = []
        for i in range(seq_len):
            if seq_ref[i]!='-':
                try:
                    seq.append(sequence[key][i])
                except IndexError:
                    seq.append('-')
        sequence_new[key] = ''.join(seq)
    seq2fasta(sequence_new,out)

## get asr
def parseRst(pref,nodelist):
    frst  = "%s.rst"%pref
    fprob = "%s.prob.json"%pref
    fasta = "%s.anc"%pref 
    alt_fasta = "%s.alt.fasta"%pref

    problist     = dict()
    mostprobable = dict() ## possible: most probable and alternatives ANC
    alternative  = dict()
    anc_msa      = dict()

    with open(frst,"r") as p:
        lines  = p.readlines()
        nlines = len(lines)
        ## get the tree with name of ancestral nodes
        idx   = 0
        identifier = "tree with node labels for Rod Page's TreeView"
        if nlines:
            while identifier not in lines[idx]:
                idx += 1
                if idx>=nlines-1:
                    break
        if nlines==0 or idx>=nlines-1:
            treeline = "dmel;"
        else:
            idx  += 1
            treeline = lines[idx][:-1]
        tree  = Tree(treeline,format=1)

        ## name of the root node
        root  = tree.get_tree_root()

        leaflist = [leaf.name.split("_") for leaf in tree]
        leaflist = sorted(leaflist,key=lambda x:int(x[0]))
        leaflist = ["_".join(leaf) for leaf in leaflist]
        #leafindex = dict()
        #for node in nodelist:
        #    leafindex[node] = [leaflist.index(leaf.name) for leaf in tree&node]

        while not lines[idx].startswith("Prob distribution at node"):
            idx += 1
        node = lines[idx].split()[4][:-1]
        problist[node] = dict()
        mostprobable[node] = dict()
        alternative[node]  = dict()
        leafindex = [leaflist.index(leaf.name) for leaf in tree&node]
        idx += 4

        ## get ASR of root node
        while not lines[idx].startswith("Prob of best state at each node"):
            if lines[idx].startswith("Prob distribution at node"):
                node = lines[idx].split()[4][:-1]
                problist[node] = dict()
                mostprobable[node] = dict()
                alternative[node]  = dict()
                leafindex = [leaflist.index(leaf.name) for leaf in tree&node]
                #print(node,len(leafindex),leafindex)
                idx += 4
            elif lines[idx].strip():
                isite,freq,seqd,*info = lines[idx][:-1].split()
                node_seqd = [seqd[i] for i in leafindex]
                gaps_seqd = [s!="-" for s in node_seqd]
                coverage  = sum(gaps_seqd)/len(node_seqd)
                #if node=="171":
                    #print(isite,coverage)
                problist[node][isite] = dict()
                if coverage<0.5:
                    aamax = "-"
                    problist[node][isite][aamax] = 1-coverage
                else:
                    aamax = "A"
                    pmax  = 0
                    for aa in info:
                        ## aa = "A(0.262)" ##
                        try:
                            aa,p = aa[:-1].split("(")
                        except ValueError:
                            print(idx)
                            print(lines[idx])
                        
                        p = float(p)
                        if p>=0.05:
                            problist[node][isite][aa] = p
                        if p>pmax:
                            pmax  = p
                            aamax = aa
                mostprobable[node][isite] = aamax
                if len(problist[node][isite])>1:
                    aalist = sorted(problist[node][isite].items(),
                                    key=lambda x:x[1],reverse=True)
                    problist[node][isite] = dict(aalist)
                    alternative[node][isite] = [aa[0] for aa in aalist]
                    alternative[node][isite].remove(aamax)
                idx += 1
            else:
                idx += 1
                    
            
        ## get most probable and alternative states ##
        alt_msa = dict()
        for node in nodelist:
            name = nodelist[node]
            mpbseq = [mostprobable[node][i] for i in mostprobable[node]] 
            altseq = [s for s in mpbseq]
            for isite in alternative[node]:
                altseq[int(isite)-1] = alternative[node][isite][0]
            anc_msa["%s"%name] = "".join(mpbseq)
            alt_msa["%s|Alternative"%name] = "".join(altseq)


        ## ancestor and probability
        seq2fasta(anc_msa,fasta)
        seq2fasta(alt_msa,alt_fasta)
        #np.savetxt(fprob,problist)
        with open(fprob,"w") as fp:
            json.dump(problist,fp)
        chop_seq(fasta,fasta)
        chop_seq(alt_fasta,alt_fasta)

if __name__ == "__main__":
    df = pd.read_csv('anc.node.txt',dtype=str)
    nodelist = dict(zip(df['node'],df['name']))
    print(nodelist)
    pref = "spr_all"
    parseRst(pref,nodelist)
