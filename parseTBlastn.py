
# Parse tblastn output into Bed format file

import os,sys
import networkx as nx

from tqdm import tqdm

class BlastTarget():
    def __init__(self,query,subject):
        self.query = query
        self.subject = subject
        self.strand = "+"
        self.start = 9999999999999999
        self.end = 0
        self.evalue = 999
        self.identity = 0
        self.positive = 0
        self.bitscore = 0
        self.organism = ""

def boolOverlap(t1,t2,cutoff):
    if t1.subject==t2.subject and t1.strand==t2.strand:
        if t1.start<t2.end+cutoff and t2.start<t1.end+cutoff:
            return True
    return False

def bestOfGroup(grp,targets,cutoff=0):
    ## get the best hits ##
    grp1 = [g for g in grp]

    best = []
    while grp1:
        min_evalue = 1E10
        for i in grp1:
            evalue = targets[i].evalue
            if evalue < min_evalue:
                min_evalue = min_evalue
                min_idx = i
        best.append(min_idx)

        grp2 = [g for g in grp1]
        for i in grp1:
            t1 = targets[min_idx]
            t2 = targets[i]
            if t1.start<t2.end+cutoff and t2.start<t1.end+cutoff:
                grp2.remove(i)
                #print(min_idx,i,grp2)
        grp1 = [g for g in grp2]

    return best

def parseBlastRes(blast="tblastn.out",out="combined.out",cut=1000,ecut=1e-5):
    with open(blast,"r") as p:
        lines = p.readlines()

    nlines = len(lines)
    #print(nlines)

    blastRes = {}
    iline = 0
    target = ""
    #while iline < nlines:
    print("Info:read> reading in tblastn %s"%blast)
    for iline in range(nlines):
        line = lines[iline]
        if line.startswith("Query="):
            query = line[:-1].split()[1]
            blastRes[query] = []

        if line.startswith(">"):
            info = line[1:-1].split()
            subject = info[0]
            organism = " ".join(info[1:3])

        if line.startswith(" Score ="):
            target  = BlastTarget(query,subject) 
            target.organism = organism

            blastRes[query].append(target)

            info = line.split()
            bitscore = float(info[2])
            evalue   = float(info[7][:-1])
            blastRes[query][-1].evalue = evalue
            blastRes[query][-1].bitscore = bitscore
            
        if line.startswith(" Identities = "):
            info = line.split()
            identity = float(info[3][1:-3])
            positive = float(info[7][1:-3])
            blastRes[query][-1].identity = identity
            blastRes[query][-1].positive = positive

        if line.startswith(" Frame ="):
            strand = "+" if "+" in line else "-"
            blastRes[query][-1].strand = strand

        if line.startswith("Sbjct"):
            info = line[:-1].split()
            start = int(info[1])
            end   = int(info[-1])
            if strand == "-":
                if start>end:
                    start,end = end,start
            if start < target.start:
                blastRes[query][-1].start = start
            if end > target.end:
                blastRes[query][-1].end = end

    ## group results by qeury ## 
    ## for each query, connect tiny tblastn hits ##
    blastResBySbjct = {}
    for query in blastRes:
        blastResBySbjct[query] = {}
        for t in blastRes[query]:
            if t.evalue<ecut:
                if t.subject in blastResBySbjct[query]:
                    blastResBySbjct[query][t.subject].append(t)
                else:
                    blastResBySbjct[query][t.subject] = [t]
    print("Info:read> Done reading %s"%blast)
    print("...")
        
    ## combine results ##
    print("Info:connect> Aggregte tiny hits by building connected graphs")
    combinedRes = {}
    for query in blastResBySbjct:
        combinedRes[query] = {}
        for subject in blastResBySbjct[query]:
            n = len(blastResBySbjct[query][subject])

            ## first compute a connected graph ##
            G = nx.Graph()
            for i in range(n):
                G.add_node(i)
            for i in range(n):
                for j in range(i,n):
                    t1 = blastResBySbjct[query][subject][i]
                    t2 = blastResBySbjct[query][subject][j]
                    #if t2.subject==t1.subject and t2.strand==t1.strand:
                    #    if t1.start-cut<t2.end+cut and t2.start-cut<t1.end+cut:
                    #        G.add_edge(i,j)
                    if boolOverlap(t1,t2,cut):
                        G.add_edge(i,j)

            ## get all connected sub graph ##
            groups = []
            for subG in nx.connected_components(G):
                groups += [sorted([s for s in subG])]


            ## combine the subgraphs ##
            combinedRes[query][subject] = []
            subjectGroups = blastResBySbjct[query][subject]
            for grp in groups:
                newtarget  = BlastTarget(query,subject) 
                starts = [subjectGroups[t].start  for t in grp]
                ends   = [subjectGroups[t].end    for t in grp]
                evalues= [subjectGroups[t].evalue for t in grp]
                strand = [subjectGroups[t].strand for t in grp][0]
                newtarget.start = min(starts)
                newtarget.end   = max(ends)
                newtarget.evalue= min(evalues)
                newtarget.strand= strand
                combinedRes[query][subject].append(newtarget)
                #print(query,subject,strand,newtarget.start,newtarget.end,newtarget.evalue)
    print("Info:connect> Done aggregting tiny hits")
    print("...")

    ## for each target, assign it to its closest GPCR family ##
    print("Info:assign> Assign tblastn hit to their closest GPCR family...")
    totaltargets = {}
    for query in combinedRes:
        for subject in combinedRes[query]:
            for t in combinedRes[query][subject]:
                chrom = t.subject+t.strand
                if chrom in totaltargets:
                    totaltargets[chrom].append(t)
                else:
                    totaltargets[chrom] = [t]

    ## start to assign and combine the subgraphs ##
    f = open(out,"w")
    combinedTarget = []
    for chrom in totaltargets:
        G = nx.Graph()
        nt = len(totaltargets[chrom])
        for i in range(nt):
            G.add_node(i)
 
        for i in tqdm(range(nt-1),desc="Info:assign> graph for %s"%chrom):
            for j in range(i+1,nt):
                t1 = totaltargets[chrom][i]
                t2 = totaltargets[chrom][j]
                if boolOverlap(t1,t2,cut):
                    G.add_edge(i,j)

        ## get all connected sub graph ##
        ## each subgraph contains tblastn targets and query GPCR proteins ##
        print("Info:assign> Done building subgraph")
        groups = []
        for subG in nx.connected_components(G):
            groups += [sorted([s for s in subG])]
        
        ## combine the subgraphs ##
        print("Info:assign> Assign hits to closest Pfam by evalue")
        ng = len(groups)
        for k in tqdm(range(ng),desc="Info:assign> Assign"):
            grp = groups[k]
            best = bestOfGroup(grp,totaltargets[chrom])
            for i in best:
                t = totaltargets[chrom][i]
                combinedTarget.append(t)
                #print(t.query,t.subject,t.strand,t.start,t.end,t.evalue)
                bedinfo = [t.subject,str(t.start),str(t.end),
                            t.query,"0",t.strand,"tblastn",str(t.evalue)]
                #info = [t.query,t.subject,t.strand,
                #        str(t.start),str(t.end),str(t.evalue)]
                info = "\t".join(bedinfo)
 
                f.write("%s\n"%info)
    f.close()
    print("Info:assign> All Done!")
    return combinedTarget

def parseMultiBlastRes(flist,cut=1000):
    lines = open(flist,"r")
    for line in lines:
        pref = line.strip()
        inp = "%s_tblastn.out"%pref
        out = "%s_combined.out"%pref
        parseBlastRes(inp,out,cut=cut)


if __name__ == "__main__":
    #parseMultiBlastRes(flist)
    cut = 50
    ecut= 9999
    inp = "sp_euk.out"
    out = "sp_euk.bed"
    #inp = "test.tblastn.out"
    #out = "test.combined.out"
    parseBlastRes(inp,out,cut=cut,ecut=ecut)
