import argparse
from contextlib import closing
import toolshed as ts
import multiprocessing
from multiprocessing import Pool, Process
from functools import partial 
from collections import defaultdict, namedtuple
import numexpr as ne
import pickle
import sys
import csv
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--ped', help='ped file')
p.add_argument('--g', help='genome file with chr lengths')
p.add_argument('-trio', help='three samples to calculate coverage on, comma separated <kid,mom,dad>')
p.add_argument('-v', '--vcf', help='vcf file')
p.add_argument('-r', '--region', help='specify region of interest <chr:start-end> OR <chr>', type=str)
p.add_argument('-c', '--cutoff', help='depth cutoff (default=10)', type=int, default=10)
p.add_argument('-s', '--shared', help='only export shared bases in trio', action='store_true')
p.add_argument('-p', '--processes', type=int, help='number of processes used to count depth for each individual -- effectively number of chromosomes that can be processed in parallel (default=1)', default=1)
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
    f2s = [k for k in ped.samples() if k.mom and k.mom.mom is not None]
    if f1: trios = f1s
    else: trios = f2s
    for i in trios:
        yield ','.join([i.sample_id, i.mom.sample_id, i.dad.sample_id])

def dp_into_array(covfile, region="1:1-100"):
    """
    Read per base depth information into a 1D
    numpy array with length equal to the length
    of the input region/chromosome. Currently this
    operates on original mosdepth output (not BED).
    """
    import tabix
    tb = tabix.open(covfile)
    c, s, e = get_regions(region)
    records = tb.query(c, s, e)
    length = e - s
    dp_array = np.zeros(length, dtype=np.uint32)
    s = 0
    for r in records:
        pos = int(r[1])
        dp = int(r[2])
        dp_array[s:pos] = dp
        s = pos
        if s % 1000000 == 0:
            print ("on chromosome %s processed %.2f percent of bases" % (c, 100 * (s / length)), file=sys.stderr)
    return dp_array

def variant_into_array(vcffile, region, trio, cutoff=10):
    """
    Iterate over a VCF, and catalog each putative DNM,
    as well as the depth at that site, and index those
    observations into two separate numpy arrays.
    """
    from cyvcf2 import VCF
    vcf = VCF(vcffile)
    smp2idx = dict(zip(vcf.samples, range(len(vcf.samples))))
    k, m, d = trio
    ki, mi, di = smp2idx[k], smp2idx[m], smp2idx[d]
    c, s, e = get_regions(region)
    length = e - s
    v_array = np.zeros(length, dtype=bool) 
    d_array = np.zeros(length, dtype=np.uint32) 
    for v in vcf(region): 
        v_len = v.end - v.start
        if v_len > 1: continue
        gts = v.gt_types
        ref, alt = v.gt_ref_depths, v.gt_alt_depths
        d_array[v.start:v.end] = ref[ki] + alt[ki]
        if (gts[mi], gts[di], gts[ki]) == (0, 0, 1):
            v_array[v.start:v.end] = True 
        if v.end % 100000 == 0:
            print ("chr %s, processed %.2f percent" % (v.CHROM, 100 * (v.end / length)), file=sys.stderr)
    return d_array, v_array

def evaluate_variant_detect_power(region, trio, vcffile, return_dict):
    d_array, v_array = variant_into_array(vcffile, region, trio)
    x, y = [], []
    for c in range(1, 100):
        v_at_cutoff = ne.evaluate("(v_array == True) & (d_array == c)")
        x.append(c)
        y.append(v_at_cutoff.sum())
    return_dict[str(region)] = (x, y)

def evaluate_dp_thresh(region, trio, return_dict, cutoff=10):
    """
    Create a 1D array of depths across a given chromosom
    in each sample from a trio. Evaluate the number of sites
    at which depth is >= some cutoff.
    """
    kd, md, dd = (dp_into_array(f+'.cov.gz', region=region) for f in trio) 
    shared_bases = ne.evaluate("(kd >= cutoff) & (md >= cutoff) & (dd >= cutoff)")
    shared_bases = shared_bases.sum()
    return_dict[str(region)] = str(shared_bases)

def plot_variant_power(v_dict):
    """
    Plot the number of putative DNMs that
    were identified at a specific depth.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    x_ax_len = max([len(v_dict[chrom][0]) for chrom in v_dict])
    xs, ys = np.zeros(x_ax_len), np.zeros(x_ax_len)

    # NOTE: will probably need special logic for X.
    wrote_x = False
    for chrom in v_dict:
        if chrom == 'X': continue
        x, y = v_dict[chrom]
        x, y = np.array(x), np.array(y)
        if not wrote_x:
            xs += x
            wrote_x = True
        ys += y
    plt.plot(xs, ys)
    power = ys[np.where(xs == 10)[0][0]]
    #plt.plot([10, 10], [0, 1], color='red', label="power = %.2f" % (power), linestyle='--')
    plt.ylabel('Number of de novo variants recovered')
    plt.xlabel('Depth (exact) required in kid')
    plt.legend()
    plt.savefig('power-plot.png')

chrom_sizes = get_chr_sizes(args.g)
def main(args):
    for trio in generate_trios(args.ped):
        k, m, d = trio.split(',')
        chroms = chrom_sizes.keys()
        if args.region: chroms = args.region
        # NOTE: this assumes gzipped files are sitting in the working directory
        if args.trio and sorted(trio.split(',')) != sorted(args.trio.split(',')): continue
        iter_trios = [(k, m, d) for i in range(len(chroms))]
        # CRITICAL to initialize dictionary this way.
        manager = multiprocessing.Manager() 
        return_dict = manager.dict()
        # Create iterables that pool.starmap needs for parallelization
        iter_dicts = [return_dict for i in range(len(chroms))]
        iter_vcf = [args.vcf for i in range(len(chroms))]
        iter_ = [i for i in zip(chroms, iter_trios, iter_dicts)]
        iter_v = [i for i in zip(chroms, iter_trios, iter_vcf, iter_dicts)]
        # If we only want to output # bp above a depth cutoff in a trio...
        if args.shared: 
            import glob
            # Make sure we haven't already created a trio-shared depth file,
            # and confirm all samples in the trio have a corresponding mosdepth file
            already_written = (k + '.trio' in glob.glob('*.trio')
            if all([t + '.cov.gz' in glob.glob('*.gz') for t in (k, m, d)]) and not already_written:
                outfile = open(k + '.trio', 'a')
                outfile.write('\t'.join(['chrom','bases_at_depth_cutoff','\n']))
            else: continue

            with Pool(processes=args.processes) as pool:
                pool.starmap(evaluate_dp_thresh, iter_)
        # ...or if we just want to plot # DNM vs depth.
        else:
            with Pool(processes=args.processes) as pool:
                pool.starmap(evaluate_variant_detect_power, iter_v)
        # convert dict proxy to 'real' dict
        d_return_dict = dict(return_dict)
        if args.shared and not already_written:
            outfile.write('\t'.join(['all', str(sum([int(x) for x in d_return_dict.values()])), '\n']))
            for chrom in d_return_dict:
                outfile.write('\t'.join([chrom, d_return_dict[chrom], '\n']))
        else: plot_variant_power(d_return_dict)

if __name__ == '__main__':
    main(args)
