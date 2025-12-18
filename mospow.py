import argparse
import pysam
import tqdm
import pandas as pd


def depth_heap(mdk, cutoff=10):
    """
    m, d, k are iterables returned from a pytabix region query.
    mdk is passed in as a tuple, so could represent a trio,
    tumor-normal pair, or a single sample
    """
    heap = [next(v) for v in mdk]
    while len(heap) > 0:
        chrom = heap[0].contig
        istart = max([h.start for h in heap])
        iend = min([h.end for h in heap])

        yn = all([float(h.name) >= cutoff for h in heap])
        yield (chrom, istart, iend, yn)
        try:
            heap = [h if h.end != iend else next(mdk[i]) for i, h in enumerate(heap)]
        except StopIteration:
            break


def evaluate_depth_heap(trio, contig, cutoff: int = 12):
    """
    Method to evaluate proportion of bases in region
    that have depth above cutoff. Smaller memory footprint
    than `evaluate_dp_thresh`, since we're not creating three
    depth arrays per function call.
    """
    k_tx, m_tx, d_tx = (
        pysam.TabixFile(covfile, parser=pysam.asBed()) for covfile in trio
    )
    kr, mr, dr = (f.fetch(contig) for f in (k_tx, m_tx, d_tx))
    shared_bases = 0
    total_bases = 0
    for c, s, e, yn in tqdm.tqdm(depth_heap((kr, mr, dr), cutoff=cutoff)):
        total_bases += (e - s)
        if yn:
            shared_bases += (e - s)
    return shared_bases, total_bases


def main(args):
    d, t = evaluate_depth_heap(args.fhs, args.contig)
    res = pd.DataFrame({"thresh_bp": [d], "total_bp": [t]})
    res["contig"] = args.contig if args.contig is not None else "UNK"
    res.to_csv(args.out, index=False)


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--fhs', nargs="+", help='mosdepth files')
    p.add_argument('--out')
    p.add_argument('-c', help='depth cutoff (default = 10)', type=int, default=10)
    p.add_argument('-contig', help='contig', default=None)
    args = p.parse_args()
    main(args)
