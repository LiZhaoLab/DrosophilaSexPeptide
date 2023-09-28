
from multiprocessing import Pool
from fasta2seq import *
import os,sys,subprocess

#db = "../cactus/MaskedGenomes/worm_genomes"
#db="/ru-auth/local/home/jpeng/scratch/refseqGenomesInsecta/refseq_genomes_insecta"
db="/ru-auth/local/home/jpeng/scratch/refseqGenomes/ref_euk_rep_genomes"
fasta_dir = "./"
alidir = "tmp_spaln_FL"
spaln = "/ru-auth/local/home/jpeng/scratch/softwares/spaln2.4.4.linux64/bin/spaln"
aminoacids = "ACDEFGHIKLMNPQRSTVWYX"

if not os.path.isdir(alidir):
    os.makedirs(alidir)

START_CODON = {
                "ATG":1
                }
STOP_CODON  = {  
                "TGA":1,
                "TAG":1,
                "TAA":1,
                }

class BedLine():
    def __init__(self,line = ""):
        bedinfo = line[0:-1].split()
        ninfo   = len(bedinfo)

        try:
            chrom   = bedinfo[0]
        except IndexError:
            chrom   = -1
        try:
            start   = int(bedinfo[1])
        except IndexError:
            start   = -1
        try:
            end     = int(bedinfo[2])
        except IndexError:
            end     = -1
        try:
            name    = bedinfo[3]
        except IndexError:
            name    = "NA"
        try:
            score   = bedinfo[4]
        except IndexError:
            score   = 0
        try:
            strand  = bedinfo[5]
        except IndexError:
            strand  = "+"

        if ninfo<4:
            valid = False
        else:
            valid = True

        others=[] if ninfo<=6 else bedinfo[6:]

        self.chrom = chrom
        self.start = start
        self.end   = end
        self.name  = name
        self.score = score
        self.strand= strand
        self.others= others
        self.valid = valid

def spalnOnBlast(sp,bedi,trb=100):
    query,subject,strand,start,end = (bedi.name,bedi.chrom,bedi.strand,
                                        bedi.start,bedi.end)
    ## first retrieve tblastn seq #
    name = query.replace("|","_")
    #name = "ele1_"+query
    #print(name)
    #chrom = sp + "|" + subject
    chrom = subject
    ref = os.path.join(fasta_dir,"%s.fasta"%name)
    reflen = len(fasta2seq(ref)[query])
    pref = "%s_%s_%s_%d_%d"%(name,sp,subject,start,end)
    out = os.path.join(alidir,"%s.fasta"%pref)
    alnlen = (end-start+1)/3
    ldiff  = max(0,reflen-alnlen)
    trb = ldiff*3 + trb
    start,end = start-trb,end+trb
    start = 1 if start<=0 else start
    start_end = "%d-%d"%(start,end)
    strand = "minus" if strand=="-" else "plus"
    cmd = "blastdbcmd -db %s -entry \"%s\" "\
          "-range %s -strand %s "\
          "-out %s -outfmt %%f"%(db,chrom,start_end,strand,out)
    print(cmd)
    subprocess.call(cmd,shell=True)

    ## then spaln ##
    cds = out
    #aln = os.path.join(alidir,"%s.spaln"%pref)
    #cmd = "%s -Q0 -O1 -S1 -TInsectDm -t2 -pw "\
    #      "-yL15 -ya0 %s %s > %s"%(spaln,cds,ref,aln)
    #subprocess.call(cmd,shell=True)

    aln = os.path.join(alidir,"%s.spaln2"%pref)
    faa = os.path.join(alidir,"%s.faa"%pref)
    cmd = "%s -Q0 -O7 -S1 -TInsectDm -t2 -pw "\
          "-yL15 -ya0 %s %s > %s"%(spaln,cds,ref,aln)
    print(cmd)
    subprocess.call(cmd,shell=True)
    with open(aln,"r") as p:
        lines = p.readlines()

    frameshift = 0
    cds_region = ""
    score = ""
    
    proseq = []
    for line in lines:
        if line.startswith(">"):
            score = line[:-1].split()[-1]
        elif line.startswith(";M"):
            frameshift += 1
        elif line.startswith(";C"):
            cds_region = cds_region + line[3:-1]
        elif line[0].isupper() :
            prostr = [s if s in aminoacids else "A" for s in line[0:-1]]
            proseq.append("".join(prostr))

    ## output predicted protein ##
    proseq = "".join(proseq)
    proseq = {name:proseq}
    seq2fasta(proseq,faa)

    ## do tmhmm analysis ##
    #cmd = "sh tmhmm.sh %s"%pref
    #subprocess.call(cmd,shell=True,
    #                stdout=subprocess.DEVNULL,
    #                stderr=subprocess.STDOUT)

    #tm_out = os.path.join(alidir,"%s.tm.out"%pref)
    nTM = 0
    #lines = open(tm_out,"r")
    #for line in lines:
    #    if "Number of predicted TMHs" in line:
    #        info = line.split()
    #        num  = info[-1]
    #        nTM  = round(float(num))
    #        break

    frameshift = str(frameshift)
    cds_region = cds_region.strip() ## remove end spaces
    cds_region = cds_region[5:-1]   ## remove join() ###

    dnaseq = fasta2seq(cds)
    key = [key for key in dnaseq][0]
    dnaseq = dnaseq[key]

    cds_out = os.path.join(alidir,"%s.cds"%pref)
    cds_seq = {}

    ## intron len and exon len ##
    exons = cds_region.split(",")
    nexon = len(exons)
    nintron = nexon - 1

    cds_str = []
    boundaries = []
    exon_len = []
    i = 0
    for ex in exons:
        estart,eend = [int(s) for s in ex.split("..")]
        boundaries.append((estart,eend))
        exon_len.append(eend-estart+1)

        ## position in genomic regions
        strand = bedi.strand
        #gstart = estart-1+start-trb if strand=="+" else end+trb-eend+1
        #gend   = eend-1+bedi.start-trb if strand=="+" else end+trb-estart+1
        gstart = estart-1+start if strand=="+" else end-eend+1
        gend   = eend-1+bedi.start if strand=="+" else end-estart+1
        cds_str.append("%d:%d"%(gstart,gend))

        cds_seq[i] = dnaseq[estart-1:eend]
        i += 1
    seq2fasta(cds_seq,cds_out) 

    exon_len = [str(s) for s in exon_len]
    exon_str = ";".join(exon_len)
    cds_str = ";".join(cds_str)
    #cds_str = "join(%s)"%cds_str if strand=="+" else "complement(%s)"%cds_str
    cds_str = "join(%s)"%cds_str

    intron_out = os.path.join(alidir,"%s.intron"%pref)
    intron_seq = {}
    intron_len = []
    for j in range(nintron):
        intron_len.append(boundaries[j+1][0] - boundaries[j][1] - 1)
        intron_seq[j] = dnaseq[boundaries[j][1]-1:boundaries[j+1][0]-1]
    intron_len = [str(s) for s in intron_len]
    intron_str = ";".join(intron_len)
    seq2fasta(intron_seq,intron_out)

    start = boundaries[0][0]-1           ## to determine start codon
    end   = boundaries[-1][1]            ## to determine stop codon

    start_codon = dnaseq[start:start+3]
    stop_codon  = dnaseq[end:end+3]

    #idx = 0
    #nlines = len(lines)
    #while idx < nlines:
    #    if lines[idx].startswith("ALIGNMENT"):
    #        idx += 1
    #        break
    #    idx += 1
    #i = 0
    #ali = "%s.ali"%pref
    #f = open(ali,"w")
    #seq1 = []
    #seq2 = []
    #while i+idx < nlines-1:
    #    #print(i,lines[i+idx]) 
    #    line = lines[i+idx]
    #    if i%2!=0:
    #        if not line.strip():
    #            seq = [line[s] for s in range(8,68) if line[s]!=" " ]
    #        else:
    #            seq = []
    #        if i%4==1:
    #            seq1 = seq1+seq
    #            print("aln:","".join(seq))
    #            f.write("aln: %s\n"%("".join(seq)))
    #        elif i%4==3:
    #            seq2 = seq2+seq
    #            print("ref:","".join(seq))
    #            f.write("ref: %s\n\n"%("".join(seq)))
    #            print("")
    #    i += 1
    #    
    #f.close()
    nexon = str(nexon)
    nintron = str(nintron)
    return bedi,score,frameshift,nexon,exon_str,nintron,intron_str,start_codon,stop_codon,cds_str,nTM


def readAndSpaln(sp,inp,out,ncpu=1,trb=100):
    ## multiprocessing ##
    if ncpu>1:
        pool = Pool(ncpu)

    ## only do spaln on tblastn targets ##
    f = open(out,"w")
    results = []
    lines = open(inp,"r")
    for line in lines:
        bedi = BedLine(line) 
        query = bedi.name
        subject = bedi.chrom
        start = bedi.start
        end = bedi.end
        strand = bedi.strand
        #if bedi.others[0] == "tblastn":
        if ncpu>1: 
            pool.apply_async(spalnOnBlast,args=(sp,bedi,trb),
                                callback=results.append)
        else:
            result = spalnOnBlast(sp,bedi,trb)
            results.append(result)
        #else:
        #    f.write(line)

    pool.close()
    pool.join()

    potential_beds = []
    for result in results:
        bedi,score,frameshift,nexon,exon_str,nintron,intron_str,start_codon,stop_codon,cds_str,nTM = result
        spaln_info = [str(score),str(frameshift),cds_str,str(nexon),exon_str,
                      str(nintron),intron_str,start_codon,stop_codon,str(nTM)]
        spaln_info = ";".join(spaln_info)
        bedinfo = [bedi.chrom,str(bedi.start),str(bedi.end),bedi.name,
                    str(bedi.score),bedi.strand] + bedi.others + [spaln_info]
        bedinfo = "\t".join(bedinfo)

        source,evalue,*_ = bedi.others
        evalue = float(evalue)

        if int(nintron)==1:
            exon_length = [int(s) for s in exon_str.split(";")]
            ratio = exon_length[0]/exon_length[1]
            intron_length = [int(s) for s in intron_str.split(";")][0]
        
        if int(frameshift)==0 and int(nintron)==1 and start_codon in START_CODON and stop_codon in STOP_CODON and abs(ratio-3)<=2.5 and intron_length<90:
            potential_beds.append(bedinfo)
        f.write(bedinfo + "\n")

    f.close()
    print("")
    print("")
    print("")
    f = open("potential_hits.bed","w")
    for bedinfo in potential_beds:
        #print(bedinfo)
        f.write(bedinfo + "\n")
    f.close()
def main(flist,ncpu=1,trb=1000):
    lines = open(flist,"r")
    for line in lines:
        sp = line.strip()
        #if sp!="ele1":
        #    inp = "%s.ele1_GPCR.fasta.bed"%sp
        #    out = "%s.ele1_GPCR.spaln.bed"%sp
        #    readAndSpaln(sp,inp,out,ncpu=ncpu)
        #sp = "nig1"
        inp = "%s.bri1_GPCR.fasta.bed"%sp
        out = "%s.bri1_GPCR.spaln.bed"%sp
        readAndSpaln(sp,inp,out,ncpu=ncpu)
        #break

if __name__ == "__main__":
    flist = "species.pref.txt"
    #flist = "test.txt"
    #flist = "others.txt"
    #main(flist,ncpu=50)
    #main(flist,ncpu=1)
    inp = "sp_euk.bed"
    out = "sp_euk.bed.spaln"
    ncpu = 36
    trb = 100
    readAndSpaln("Dmel",inp,out,ncpu=ncpu,trb=trb)
