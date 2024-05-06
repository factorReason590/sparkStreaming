
import argparse as arg
import sys
import re

##### Args #####
if __name__ == "__main__":
    parser = arg.ArgumentParser(description='''A simple program to remove duplicate reads from a sam file.\n
                                               Example: deduper.py -u {known_umi_file} -p -f {file} | sort -u -k 1,1 | deduper.py -d > {file}_deduped.sam
                                               ''')
    # Arguments
    parser.add_argument("-f", "--file", help="SAM file, else reads stdin.", type=str)
    parser.add_argument("-p", "--paired", help="Paired reads. Reads need to be 'paired'.", action='store_true')
    parser.add_argument("-u", "--umi", help="Has UMI at end of name.", type=str)
    parser.add_argument("-d", "--decomplex", help="Remove the sorting tag from a processed file.", action='store_true')
    # Parse arguments
    ARGS = parser.parse_args()
    
# Constants
soft_re=re.compile('^(\d)+S')

def yield_sam(file, paired=False):
    '''Make a generator that yields (seq_header, seq) for each entry''' 
    fhs=open(file, 'r')
    l=fhs.readline().strip()
    while l != '':
        if l.startswith('@'):   # If meta data
            sys.stdout.write(l+'\n')
            l=fhs.readline().strip() #skip
        else:
            if paired==False:
                split=l.split('\t')
                yield [l, split] # yield data
                l=fhs.readline().strip()
            else:
                line1=l
                split1=l.split('\t')
                l=fhs.readline().strip()
                split=l.split('\t')
                yield [[line1, split1],[l, split]] # yield data
                l=fhs.readline().strip()

def soft_adj(cigar, pos):
    match=soft_re.match(cigar)
    if match==None:
        return int(pos)
    else:
        return int(pos)+int(match.group(1))
    
if ARGS.umi is not None:
    with open(ARGS.umi, 'r') as ufh:
        known_umis=[x.strip() for x in ufh]
    
if ARGS.decomplex is False:    
    if ARGS.paired is False:    
        for read in sys.stdin if ARGS.file is None else yield_sam(ARGS.file, ARGS.paired):
            if ARGS.umi is None or read[1][0].rpartition(':')[2] in known_umis:
                sys.stdout.write('{chrom}:{pos}:{strand}:{umi}\t{line}\n'.format(
                                 line=read[0],
                                 chrom=read[1][2],
                                 pos=str(soft_adj(read[1][5], read[1][3])),
                                 strand='-' if int(read[1][1]) & 16 is not 0 else '+',
                                 umi=read[1][0].rpartition(':')[2] 
                                 ))
    else:
        for read in sys.stdin if ARGS.file is None else yield_sam(ARGS.file, ARGS.paired):
            if ARGS.umi is None or (read[0][1][0].rpartition(':')[2] in known_umis and read[1][1][0].rpartition(':')[2] in known_umis):
                sys.stdout.write('{chrom0}:{pos0}:{strand0}:{umi0}^{chrom1}:{pos1}:{strand1}:{umi1}\t{line0}\x0b{line1}\n'.format(
                                 line0=read[0][0], line1=read[1][0],
                                 chrom0=read[0][1][2], chrom1=read[1][1][2],
                                 pos0=str(soft_adj(read[0][1][5], read[0][1][3])), pos1=str(soft_adj(read[1][1][5], read[1][1][3])),
                                 strand0='-' if int(read[0][1][1]) & 16 is not 0 else '+', strand1='-' if int(read[1][1][1]) & 16 is not 0 else '+',
                                 umi0=read[0][1][0].rpartition(':')[2] , umi1=read[1][1][0].rpartition(':')[2]
                                 ))
else: 
    for line in sys.stdin if ARGS.file is None else yield_sam(ARGS.file, ARGS.paired):
        sys.stdout.writelines([x+'\n' for x in line.strip().partition('\t')[2].split('\x0b')])
        
sys.stdout.flush()