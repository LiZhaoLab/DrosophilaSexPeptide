

from scipy.spatial.distance import cdist
from pdb2seq import *

import os,sys,shutil,json,pickle,matplotlib,argparse
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

mass = {'H': 1.008,
        'C': 12.01,
        'N': 14.01,
        'O': 16.00,
        'S': 32.06,
        'P': 30.97,
        'M': 0.000,
        'ZN':65.0 }

elements = {
         '1H2\'' :'H',
         '1H5\'' :'H',
         '2H2\'' :'H',
         '2H5\'' :'H',
         '1HD1':'H',
         '1HD2':'H',
         '2HD1':'H',
         '2HD2':'H',
         '3HD1':'H',
         '3HD2':'H',
         '1HE2':'H',
         '2HE2':'H',
         '1HG1':'H',
         '1HG2':'H',
         '1HH1':'H',
         '1HH2':'H',
         '2HG1':'H',
         '2HG2':'H',
         '3HG1':'H',
         '3HG2':'H',
         '2HH1':'H',
         '2HH2':'H',
         '3HG2':'H',
         '0C21':'C',
         '1C21':'C',
         '2C21':'C',
         '3C21':'C',
         '4C21':'C',
         '5C21':'C',
         '6C21':'C',
         '7C21':'C',
         '8C21':'C',
         '0C31':'C',
         '1C31':'C',
         '2C31':'C',
         '3C31':'C',
         '4C31':'C',
         '5C31':'C',
         '6C31':'C',
         'ZN':'ZN'
           }

def get_elem(atname):
    if atname[0] in mass.keys():
        return atname[0]
    else:
        return elements[atname]

class atom:
    def __init__(self,atid,atname,resname,resi,coori):
        self.atid = atid
        self.atname = atname
        self.resname = resname
        self.resi = int(resi)
        self.coori = coori

# adopted from Burke et al, NSMB, 2023 #
def sigmoid(X,dcut):
    if dcut <=5:
        #5 SpearmanrResult(correlation=0.765, pvalue=4.75e-280)
        parm=[6.96234405e-01, 2.35483775e+02, 2.25322970e-02, 2.88445245e-02]
        #0.7805034405869632
    elif dcut <=6:
        #6 SpearmanrResult(correlation=0.771, pvalue=2.71e-287)
        parm=[7.02605033e-01, 2.91749822e+02, 2.70621128e-02, 2.25416051e-02]
        #0.7871982094514278
    elif dcut <=7:
        #7 SpearmanrResult(correlation=0.771, pvalue=2.24e-287)
        parm=[7.06385097e-01, 3.32456259e+02, 2.97005237e-02, 2.24488132e-02]
        #0.7859609807320201
    elif dcut <=8:
        #8 SpearmanrResult(correlation=0.763, pvalue=2.34e-278)
        parm=[7.18442739e-01,3.60791204e+02,3.01635944e-02, 2.04076969e-02]
        #0.7764648775754815
    elif dcut <=9:
        #9 SpearmanrResult(correlation=0.750, pvalue=4.54e-263)
        parm=[7.23328534e-01, 3.80036094e+02, 3.06316084e-02, 1.98471192e-02]
        #0.7608417399783565
    elif dcut <=10:
        #10 SpearmanrResult(correlation=0.733, pvalue=7.99e-246)
        parm=[7.20293782e-01, 3.95627723e+02, 3.15235037e-02, 2.37304238e-02]
        #0.7431426093979494
    elif dcut <=11:
        #11 SpearmanrResult(correlation=0.713, pvalue=1.75e-226)
        parm=[7.22015998e-01, 4.09095024e+02, 3.11905555e-02, 2.59467513e-02]
        #0.7219615906164123
    else:
        #12 SpearmanrResult(correlation=0.694, pvalue=9.28e-210)
        parm=[7.20555781e-01, 4.21033584e+02, 3.09024241e-02, 2.88659629e-02]
        #0.7023000652310362
    
    L,x0,k,b = parm
    Q = L / (1 + np.exp(-k*(X-x0)))+b
    return (Q)

def pdockq(fpdb,fpkl,dcut):
    ## load pickle file ##
    data = pickle.load(open(fpkl,"rb"))

    ## get interface ##
    # read pdb #
    lines = open(fpdb,"r")
    pdbatoms = []
    nat = 0
    chains = {}
    for line in lines:
        if line.startswith('ATOM '):
            atid = int(line[6:11])
            atname = line[12:16].strip()
            resname = line[17:21].strip()
            resi = line[22:26].strip()
            chain = line[21]
            resname = "/".join([resname,resi,chain])
            x = line[30:38]
            y = line[38:46]
            z = line[46:54]
            coori = [float(s) for s in [x,y,z]]
            atomi = atom(atid,atname,resname,resi,coori)
            pdbatoms += [atomi]

            ##if atname in ["CA","CB","N"]:
            if get_elem(atname)!="H":
            #if atname:
                if chain not in chains:
                    chains[chain] = [atomi]
                else:
                    chains[chain].append(atomi)

            nat += 1

    chainkeys = [key for key in chains]
    chainA = chains[chainkeys[0]]
    chainB = chains[chainkeys[1]]

    coorA = np.array([ai.coori for ai in chainA])
    coorB = np.array([ai.coori for ai in chainB])

    # compute distance for two list #
    dmat = cdist(coorA,coorB)
    #print(np.shape(dmat))

    # extract pairs with distance < dcut #
    imat = np.where(dmat<dcut)
    indexA = imat[0]
    indexB = imat[1]

    residuesA = np.array([ai.resi for ai in chainA])
    residuesB = np.array([ai.resi for ai in chainB])
    interfaceA = residuesA[indexA]
    interfaceB = residuesB[indexB]
    interfaceA = sorted(np.unique(interfaceA-1))
    interfaceB = sorted(np.unique(interfaceB-1))
    residuesA = sorted(np.unique(residuesA))
    residuesB = sorted(np.unique(residuesB))
    nresA = len(residuesA)
    nresB = len(residuesB)

    ## number of interface residues ##
    num_interface = len(interfaceA) + len(interfaceB)
    ifA_str = numlist_to_str(interfaceA)
    ifB_str = numlist_to_str(interfaceB)
    interacting_res = 'A:%s|B:%s'%(ifA_str,ifB_str)

    plddt = data['plddt']
    iptm  = data['iptm']

    if num_interface:
        plddtA = plddt[interfaceA]
        plddtB = plddt[[i+nresA for i in interfaceB]]
        ave_plddt = np.mean([p for p in plddtA]+[p for p in plddtB])
    else:
        ave_plddt = 0

    ## determine X ##
    if num_interface:
        X = ave_plddt*np.log(num_interface)
    else:
        X = 0

    # get pDockQ #
    #Q = 0.707/(1+np.exp(-0.03148*X+12.2161288)) + 0.03138
    Q = sigmoid(X,dcut)
    
    return(Q,interacting_res)
    

def numlist_to_str(numlist):
    if numlist:
        numstr = []
        idx = 0
        r0 = numlist[0]
        while idx<len(numlist)-1:
            r1 = numlist[idx+1]
            if r1-numlist[idx]>5:
                numstr.append('%d-%d'%(r0,numlist[idx]))
                r0 = r1
            idx += 1
        if idx:
            numstr.append('%d-%d'%(r0,r1))
        else:
            numstr.append('%d-%d'%(r0,r0))
        return ';'.join(numstr)
    else:
        return ''

def human_format(num):
    num = float('{:.3g}'.format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return '{}{}'.format('{:f}'.format(num).rstrip('0').rstrip('.'), ['', 'K', 'M', 'B', 'T'][magnitude])

def plot(modelname,showfigure,dcut):
    ## check if the simulation is finished ##
    rankingf = "%s/ranking_debug.json"%modelname
    if os.path.isfile(rankingf):
        data = json.load(open(rankingf))
        topm = data["order"][0]
        pkl = "%s/result_%s.pkl"%(modelname,topm)
        pdb = "%s/ranked_0.pdb"%modelname
        shutil.copy(pdb,"%s.pdb"%modelname)

        ## determine if it's multimer or monomer ##
        sequence = pdb2seq(pdb)
        if len(sequence)>1:
            model = "multimer"
            prolen = [len(sequence[key]) for key in sequence]
        else:
            model = "monomer"
        #print("Prediction finished, plotting best model")

    ## if the prediction is not finished, check if there are any outputs ##
    else:
        max_pdb = ""
        max_pkl = ""
        max_score = 0
        for f in os.listdir(modelname):
            if f.endswith("pkl"):
                pkl = "%s/%s"%(modelname,f)
                data = pickle.load(open(pkl,"rb"))
                if "iptm" in data:
                    key = "iptm"
                elif "ptm" in data:
                    key = "ptm"
                else:
                    key = ""
                if key:
                    if data[key] > max_score:
                        max_pkl = pkl
                        max_score = data[key]
        if max_pkl == "":
            sys.stderr.write("Prediction with no outputs, quit\n")
            sys.exit(0)

        ## if there are outputs, check the first output model ##
        pkl = max_pkl
        max_name = os.path.basename(pkl)[7:-4]
        max_pdb = "%s/relaxed_%s.pdb"%(modelname,max_name)
        if os.path.isfile(max_pdb):
            pdb = max_pdb
        else:
            pdb = "%s/unrelaxed_%s.pdb"%(modelname,max_name)
        if showfigure=="Yes":
            print("Not finished, plot current best model: %s"%max_name)
        sequence = pdb2seq(pdb)
        if len(sequence)>1:
            model = "multimer"
            prolen = [len(sequence[key]) for key in sequence]
        else:
            model = "monomer"

        ## still extract current best model ##
        shutil.copy(pdb,"%s.pdb"%modelname)

    data = pickle.load(open(pkl,"rb"))
    PAE  = data["predicted_aligned_error"]
    maxp = data["max_predicted_aligned_error"]
    ndim = PAE.shape[0]
    min_PAE = maxp

    if model=="multimer":
        ppi   = "ppi: no?"
        score = data["iptm"] * 100
        diag_PAE1 = PAE[0:prolen[0],prolen[0]:ndim]
        min1 = np.min(diag_PAE1)
        #chainA_set1 = np.where(diag_PAE1<5)[0].tolist()
        #chainB_set1 = np.where(diag_PAE1<5)[1].tolist()
        diag_PAE2 = PAE[prolen[0]:ndim,0:prolen[0]].T
        #chainA_set2 = np.where(diag_PAE2<5)[0].tolist()
        #chainB_set2 = np.where(diag_PAE2<5)[1].tolist()
        #chainA_combined = sorted(list(set(chainA_set1+chainA_set2)))
        #chainB_combined = sorted(list(set(chainB_set1+chainB_set2)))
        #chainA_str = numlist_to_str(chainA_combined)
        #chainB_str = numlist_to_str(chainB_combined)
        #interacting_res = 'A:%s|B:%s'%(chainA_str,chainB_str)
        diag_PAE = np.concatenate((diag_PAE1,diag_PAE2))
        min2 = np.min(diag_PAE2)
        min_PAE = np.max([min1,min2])

        Q,interacting_res = pdockq(pdb,pkl,dcut)
        #min_PAE = np.min(diag_PAE)

        #if score>50 or min_PAE<5:
        #    ppi = "ppi: possible."
        #if score>50 and min_PAE<5:
        #    ppi = "ppi: confident"
        #if score>70 and min_PAE<5:
        #    ppi = "ppi: high conf"
        ppi = "ppi: no?"
        if Q>=0.23:
            ppi = "ppi: confident"    
        elif Q>=0.5:
            ppi = "ppi: highconf"
        out_str = "%s ipTM: %2d PAE: %.1f pDockQ: %.3f %s %s"%(modelname,score,
                                                min_PAE,Q,interacting_res,ppi)
        print(out_str)
    else:
        score = data["plddt"]
        print(modelname,"PLDDT:%2d"%np.mean(score))

    i5   = 1
    l5   = 50
    n5   = ndim//(i5*l5)
    while n5>7:
        if i5>=2:
          i5 += 2
        else:
          i5 += 1
        n5 = ndim//(i5*l5)
        #print(n5,i5,i5*l5)
    n5   = ndim//(i5*l5)
    ticks= [n*(i5*l5) for n in range(n5+1)]
    ticklabels = [human_format(n) for n in ticks]
    
    ## sns heatmap ##
    if showfigure=="Yes":
        fig,ax = plt.subplots(figsize=(4,3.1))
        darkgreen  = "#004b00"
        lightgreen = "#fafffa"
        colors = [darkgreen,lightgreen]
        cmap = matplotlib.colors.LinearSegmentedColormap.from_list("",colors)
        #sns.heatmap(PAE,cmap="Greens_r",vmin=0,vmax=maxp,
        #            ax=ax)
        #plt.imshow(PAE,cmap="Greens_r",vmin=0,vmax=maxp)
        cmap = "bwr"
        sns.heatmap(PAE,cmap=cmap,vmin=0,vmax=maxp,
                    ax=ax,cbar=True,
                    cbar_kws={"shrink":0.85,"label":"Expected position error (Å)"})
        lc = "k"
        #ax.axhline(y=0, color=lc,linewidth=1)
        #ax.axhline(y=ndim, color=lc,linewidth=1)
        #ax.axvline(x=0, color=lc,linewidth=1)
        #ax.axvline(x=ndim, color=lc,linewidth=1)
        if model=="multimer":
            for i in range(len(prolen)-1):
                ax.axvline(x=prolen[i]-.5, color=lc,linewidth=1.5,linestyle="--")
                ax.axhline(y=prolen[i]-.5, color=lc,linewidth=1.5,linestyle="--")
        ax.set_ylabel("Aligned residue") 
        ax.set_xlabel("Scored residue") 
        ax.tick_params(axis='both', which='both', length=1.5)
        plt.xticks(ticks,ticklabels,rotation=0)
        plt.yticks(ticks,ticklabels)
        #plt.axis('on')
        fig.tight_layout()
        #plt.savefig("%s_pae.pdf"%modelname)
        plt.savefig("%s_pae.png"%modelname,dpi=300)
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='PAE, pdockQ, and iPTM of AlphaFold models')
    parser.add_argument('-m','--model',
                        help='Name prefix of the model to plot',
                        required=True)
    parser.add_argument('-s','--show',
                        help='Whether to show and save figure',
                        default="Yes",choices=["Yes", "No"])
    parser.add_argument('-d','--dcut',
                        help='cutoff to determine interactions',
                        default=7)
    args = parser.parse_args()
    plot(args.model,args.show,args.dcut)

if __name__ == "__main__":
    main()
