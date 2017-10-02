import argparse
from contextlib import closing
import toolshed as ts
from multiprocessing import Pool, Manager
from functools import partial 
from collections import defaultdict, namedtuple
import sys
import tabix
import numexpr as ne
import glob
import csv
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--ped', help='ped file')
p.add_argument('--genome', help='genome file with chr lengths')
p.add_argument('-trio', help='three samples to calculate coverage on, comma separated <kid,mom,dad>', type=str)
p.add_argument('-d', '--directory', help='directory containing bgzipped and indexed coverage files', default='./')
p.add_argument('-r', '--region', help='specify region of interest <chr:start-end> OR <chr>', type=str)
p.add_argument('-b', '--bed', help='instead of -r, pass in bed file of regions (i.e., exons, genes, TEs)')
p.add_argument('-s', '--suffix', help="specify the suffix for files containing mosdepth output (default='per-base.bed.gz'", default='per-base.bed.gz')
p.add_argument('-c', '--cutoff', help='depth cutoff (default = 10)', type=int, default=10)
p.add_argument('-p', '--processes', type=int, help='number of processes used to count depth for each individual -- effectively number of chromosomes that can be processed in parallel (default = 1)', default=1)
args = p.parse_args()

def get_chr_sizes(genome_file):
    """
    Generate a dictionary with the bp lengths
    of each chromosome in the provided .genome file.
    """
    from collections import defaultdict
    import csv
    chrom_sizes = defaultdict(int)
    with open(genome_file, 'r') as tsvfile:
        genome_file = csv.reader(tsvfile, delimiter='\t')
        for line in genome_file:
            chrom_sizes[line[0]] = int(line[1])
    return chrom_sizes

def get_regions_from_bed(b):
    region_list = []
    with open(b, 'r') as tsvfile:
        bedfile = csv.reader(tsvfile, delimiter='\t')
        for l in bedfile:
            if int(l[2]) < int(l[1]):
                continue
            region_list.append(l[0] + ':' + l[1] + '-' + l[2])
    return region_list
            
def get_regions(r):
    """
    If applicable, convert the provided args.region
    input into a <chr,start,end> interval for tabix.
    """
    if len(r) in [1,2]:
        c, s, e = r, 1, chrom_sizes[r]
    else:
        f, s = r.split(':')[0], r.split(':')[1]
        c, s, e = f, s.split('-')[0], s.split('-')[1]
    return c, int(s), int(e)

def depth_heap(mdk, cutoff=10):
    """
    m, d, k are iterables returned from a pytabix region query.
    mdk is passed in as a tuple, so could represent a trio,
    tumor-normal pair, or a single sample
    """

    heap = [next(v) for v in mdk]

    while len(heap) > 0:
        chrom = heap[0][0]
        istart = max(int(h[1]) for h in heap)
        iend = min(int(h[2]) for h in heap)

        yn = all(int(h[-1]) >= cutoff for h in heap)
        yield (chrom, istart, iend, yn)

        try:
            heap = [h if int(h[2]) != iend else next(mdk[i]) for i, h in enumerate(heap)]
        except StopIteration:
            break

def generate_trios(pedfile, f1=True):
    """
    Given a PED file, specify whether you want
    to output trios w/r/t the F1 or F2 generation
    (i.e., whether the kid in each trio is an F1 or F2).
    """
    from peddy import Ped
    ped = Ped(pedfile)
    p0s = [k for k in ped.samples() if k.mom is None]
    f1s = [k for k in ped.samples() if k.mom and k.dad and k.mom.mom is None]
    f2s = [k for k in ped.samples() if k.mom and k.mom.mom]
    if f1: trios = f1s
    else: trios = f2s
    for i in trios:
        yield (i.sample_id, i.mom.sample_id, i.dad.sample_id)

def dp_into_array(covfile, region='1:1-100'):
    """
    Read per base depth information into a 1D
    numpy array with length equal to the length
    of the input region/chromosome. 
    """
    import tabix
    tb = tabix.open(covfile)
    c, s, e = get_regions(region)
    records = tb.query(c, s, e)
    length = e - s
    dp_array = np.zeros(length, dtype=np.uint32)
    bed = False
    for r in records:
        # Check if using mosdepth BED output or
        # three-column original output
        if len(r) > 3: bed = True
        if bed: 
            start, end, dp = (int(r[n]) for n in range(1,4))
            if dp == 0: continue
            dp_array[start:end] = dp
            continue
        pos, dp = int(r[1]), int(r[2])
        if dp == 0: continue
        dp_array[s:pos] = dp
        s = pos
    return dp_array

def evaluate_depth_heap(trio, region, return_dict, cutoff):
    """
    Method to evaluate proportion of bases in region
    that have depth above cutoff. Smaller memory footprint
    than `evaluate_dp_thresh`, since we're not creating three
    depth arrays per function call.
    """
    import tabix
    k_tx, m_tx, d_tx = (tabix.open(covfile) for covfile in trio)
    c, s, e = get_regions(region)
    kr, mr, dr = (f.query(c, s, e) for f in (k_tx, m_tx, d_tx))

    shared_bases = 0
    for r in depth_heap((kr, mr, dr), cutoff=cutoff):
        if r[-1] == True: 
            shared_bases += (r[2] - r[1])
    return_dict[str(region)] = str(shared_bases)

def evaluate_dp_thresh(region, trio, return_dict, cutoff):
    ka, ma, da = (dp_into_array(c, region=region) for c in trio)
    shared_bases = numexpr("ka >= cutoff & ma >= cutoff & da >= cutoff")
    shared_bases = shared_bases.sum()
    return_dict[str(region)] = str(shared_bases)

chrom_sizes = get_chr_sizes(args.genome)
def main(args):
    chroms = chrom_sizes.keys()
    wd = args.directory
    for trio in generate_trios(args.ped):
        k, m, d = trio
        if args.trio and sorted((k, m, d)) != sorted(args.trio.split(',')): 
            continue
        if args.region and args.bed: 
            raise NameError("can't specify region and bed")
        elif args.region: chroms = [args.region]
        elif args.bed: chroms = get_regions_from_bed(args.bed)
        else: pass
        k, m, d = (wd + samp + '.' + args.suffix for samp in (k, m, d))
        # create iterable that pool.starmap needs for parallelization
        # across provided chromosomes/regions
        dict_proxy = Manager().dict()
        iterable = [((k, m, d), c, dict_proxy, args.cutoff) for c in chroms]
        # confirm that all samples in the trio have a corresponding depth file.
        if all([x in glob.glob(wd + '*.gz') for x in (k, m, d)]):
            outfile = open(k + '.trio', 'w')
            outfile.write('\t'.join(['chrom','bases_at_depth_cutoff','\n']))
        else: continue
        with Pool(processes=args.processes) as pool:
            pool.starmap(evaluate_depth_heap, iterable)
        # convert dict proxy to 'real' dict
        dict_ = dict(dict_proxy)
        combined_bp = sum([int(x) for x in dict_.values()])
        outfile.write('\t'.join(['all', str(combined_bp), '\n']))
        for chrom in dict_:
            outfile.write('\t'.join([chrom, dict_[chrom], '\n']))

if __name__ == '__main__':
    main(args)
