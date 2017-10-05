## mospow

`mospow` is a tool for quickly determining your power to detect mutation in WGS alignments

specifically, `mospow` calculates the fraction of the genome covered by at least X reads in a single individual, across a trio of samples, or across a tumor-normal pair

```
mospow [options] BED \
    -d --depth_cutoff           depth cutoff (default = 10)
    -b --bed_regions            BED file of regions to restrict power calculation 
    -r --region                 instead of -b, single region <chr:start-end>
```



