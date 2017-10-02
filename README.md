## mospow

`mospow` is a tool for quickly determining your power to detect mutation in WGS alignments

specifically, `mospow` calculates the fraction of the genome covered by at least X reads in a single individual, across a trio of samples, or across a tumor-normal pair

```
python mospow.py \
    --ped       PED file [REQUIRED]
    --genome    genome file, formatted as <chr:length> [REQUIRED]
    -trio       trio of specific sample names, comma-separated <8410,8402,8403>
    -d          directory containing gzipped and index mosdepth output (default = '.')
    -r          specify region of interest <chr:start-end> or <chr>
    -s          file suffix for mosdepth output (default = 'per-base.bed.gz')
    -c          depth cutoff (default = 10)
    -p          number of processes (default = 1)
```



