## extract protein/DNA/RNA sequences from PDB

import os,sys

TRACE = ["CA","O5\'"]
PRIMARY_ALT = [" ","A"]
MAX_LENGTH = 99999

key3to1 = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
    'MSE': 'M',   
            }

def pdb2seq(inp):
    if not os.path.exists(inp):
        print("pdb file: %s not exist"%inp)
        sys.exit(0)

    sequence = {}
    maxresid = {}
    lines = open(inp,"r")
    for line in lines:
        if line.startswith("ATOM "):
            atom = line[12:16].strip()
            if atom in TRACE:
                resid = int(line[22:26])
                resnm = line[17:20].strip()
                alt_r = line[16]
                chain = line[21]
                if alt_r in PRIMARY_ALT:
                    if chain not in sequence:
                        sequence[chain] = ["-" for i in range(MAX_LENGTH)]
                    sequence[chain][resid-1] = key3to1[resnm]
                    maxresid[chain] = resid

    for chain in sequence:
        sequence[chain] = "".join(sequence[chain][:maxresid[chain]])

    return sequence

if __name__ == "__main__":
    inp = sys.argv[1]
    sequence = pdb2seq(inp)
    print(sequence)
