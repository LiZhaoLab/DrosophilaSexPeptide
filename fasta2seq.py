## read fasta format file into python dictionary

import gzip,os,sys

def fasta2seq(inp):

    if not os.path.exists(inp):
        print("WARN: Fasta format input <%s> does not exist"%inp)
        lines = []
        #sys.exit(0)

    else:
        if inp.endswith(".gz"):
            lines = gzip.open(inp,"rb")
        else:        
            lines = open(inp,"r")

    sequence = {}
    for line in lines:
        if line.startswith(">"):
            key = line[1:-1].strip()
            sequence[key] = []
        else:
            sequence[key] += line[:-1].strip()

    for key in sequence:
        sequence[key] = "".join(sequence[key])

    return sequence


def seq2fasta(seq,fasta,linebreak=False,compress=False):
    if compress:
        f = gzip.open("%s.gz"%fasta,'wb')
    else:
        f = open(fasta,'w')

    ## line break for txt and gz format ##
    lb = '\n'
    lb_byte = lb.encode('utf-8')

    for key in seq:
        seqi = seq[key]

        header = '>%s\n'%key   
        if compress:
            header = header.encode('utf-8')
        f.write(header)
        n = 1

        if linebreak:
            for s in seqi:
                if compress:
                    s = s.encode('utf-8')
                f.write(s)
                    
                if n%80 == 0:
                    if compress:
                        lb = lb_byte
                    f.write(lb)
                n += 1
        else:
            if compress:
                seqi = seqi.encode('utf-8')
            f.write(seqi)

        if compress:
            lb = lb_byte
        f.write(lb)

def seq2phylip(seq,outfile,name=''):
    phylines = dict([[key,'%-10s'%"%s%s"%(name,key)] for key in seq])

    for key in seq:
        key = key
        break
    align_length = len(seq[key])
    new_length   = 0
    for j in range(align_length):
        gaps = [seq[key][j] for key in seq]
        if not all([s=="-" for s in gaps]):
            new_length += 1
            for key in seq:
                phylines[key] += seq[key][j]
    for key in seq:
        phylines[key] += '\n'

    with open(outfile,'w') as p:
        p.write("%10d %10d\n"%(len(phylines),new_length))
        for key in phylines:
            p.write(phylines[key])


if __name__ == '__main__':
    fasta = '/ru-auth/local/home/jpeng/scratch/DrosophilaProteins/LiftOver/KaKsTest/dmel_group/fasta/FBgn0003495.fa'
    sequence = fasta2seq(fasta)
    print(sequence)

    #print('This is OK?')
    #seq_aa = translate(sequence)
    #print(seq_aa)

    #from pdb2seq import printscreen
    #printscreen(seq_aa)
