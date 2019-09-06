import os
import pysam

def call_snp(bdict,mychr,mypos,myref,myalt):
    snp_counts = []
    rc,ac,nc,hc,fc = 0,0,0,0,0
    for x in sorted(list(bdict)):
        cdict = {'A':0,'T':0,'G':0,'C':0}
        total = 0
        for y in bdict[x].pileup(mychr,mypos-1,mypos):
            if not y.reference_pos == mypos-1: continue
            for z in y.pileups:
                if not z.is_del and not z.is_refskip: 
                    cdict[z.alignment.query_sequence[z.query_position]] += 1
                    total += 1
        if total < 5:
            snp_counts.append('N')
            nc += 1
        elif cdict[myref]/float(total) > .9:
            snp_counts.append(myref)
            rc += 1
        elif cdict[myalt]/float(total) > .9:
            snp_counts.append(myalt)
            ac += 1
        elif (cdict[myalt]+cdict[myref])/float(total) > .9 and cdict[myalt]/float(cdict[myref]+cdict[myalt]) > .2 and cdict[myref]/float(cdict[myalt]+cdict[myref]) > .2:
            snp_counts.append(myref + "/" + myalt)
            hc += 1
        else:
            snp_counts.append('X')
            fc += 1
    snp_counts.extend(map(str,[rc,ac,hc,nc,fc,((ac*2)+hc)/(2*float(ac+rc+hc)+.1),(rc+hc+ac)/float(192),hc/float(hc+rc+ac+.1)]))
    return snp_counts

myfiles = os.listdir("Bam_files_perbc")
bdict = {}
for f in myfiles:
    if f[-3:] != 'bam': continue
    mystub = f.split('.')[0]
    bdict[mystub] = pysam.AlignmentFile("Bam_files_perbc/{0}".format(f),'rb')

import os
import subprocess as sp
mycount = 0
allcount = 0
print ",".join(['blah','blah','chr','pos','ref','alt']+sorted(list(bdict)))
snp_sites = os.listdir("samtools_snps")
for f in snp_sites:
    fh = open("tempfile",'w')
    proc = sp.Popen(['bcftools','view','samtools_snps/{0}'.format(f)],stdout=fh)
    proc.wait()
    fh.close()
    fh = open("tempfile")
    for c in fh:
        if len(c) < 2: continue
        if c[0] == '#': continue
        d = c.split('\t')
        if len(d) < 2: continue
        if float(d[5]) < 50: continue
        if ',' in d[4]: continue
        calls = call_snp(bdict,d[0],int(d[1]),d[3],d[4])
        plist = ['','',d[0],d[1],d[3],d[4]]
        plist.extend(calls)
        allcount += 1
        if float(plist[-1]) > .6: continue
        if float(plist[-2]) < .3: continue
        if float(plist[-3]) < .15: continue
        if float(plist[-3]) > .85: continue
        mycount += 1
        plist[0] = str(mycount)
        plist[1] = str(allcount)
        print ",".join(plist)
