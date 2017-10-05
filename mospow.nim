import hts
import strutils
import kexpr
import os

type
  interval = tuple
    chrom: string
    start: int
    stop: int
    depth: int

  Gen = iterator (): interval

# https://forum.nim-lang.org/t/435
proc fastSubStr(dest: var string; src: string, a, b: int) {.inline.} =
  # once the stdlib uses TR macros, these things should not be necessary
  template `+!`(src, a): expr = cast[pointer](cast[int](cstring(src)) + a)
  setLen(dest, b-a)
  copyMem(addr dest[0], src+!a, b-a)


proc gen_from_bed(per_base_bed:string, chrom:string, rstart:int, rstop:int, depth_col:int=3): Gen =
  result = iterator(): interval {.inline.} =
    var col = depth_col
    var b = ropen_bgzi(per_base_bed)
    var pos0, pos1: int
    var tmp = new_string_of_cap(20)

    var start, stop, depth: int
    # a lot of extra stuff in here to minimize garbage creation and
    # copying
    for line in b.query(chrom, rstart, rstop):
      pos0 = line.find('\t')
      pos1 = line.find('\t', pos0+1)
      fastSubstr(tmp, line, pos0+1, pos1)
      start = parseInt(tmp)

      pos0 = line.find('\t', pos1+1)
      fastSubStr(tmp, line, pos1+1, pos0)
      stop = parseInt(tmp)

      pos1 = line.find('\t', pos0+1)
      if pos1 == -1:
        pos1 = line.len
      fastSubStr(tmp, line, pos0+1, pos1)
      depth = parseInt(tmp)
      yield (chrom, start, stop, depth)

iterator meets(depth_cutoff:int, depths: seq[Gen]): interval =
  var heap = new_seq[interval](depths.len)
  var maxi = new_seq[int](depths.len)
  for i in 0..<depths.len:
    heap[i] = depths[i]()

  var done = false
  var chrom: string
  while not done:
    maxi.set_len(0)
    chrom = heap[0].chrom
    var istart = heap[0].start
    var istop = heap[0].stop
    var allok = true
    for i, iv in heap:
      if iv.start > istart:
        istart = iv.start
      if iv.stop < istop:
        # this is the new min stop so we clear
        maxi.set_len(0)
        istop = iv.stop
        # and save
        maxi.add(i)
      elif iv.stop == istop:
        # or we add to the existing set
        maxi.add(i)
      if iv.depth < depth_cutoff:
        allok = false

    if allok:
      yield (chrom, istart, istop, 1)
    # any site that ended at the max end
    # we replace with the next value.
    for i in maxi:
      heap[i] = depths[i]()
      if finished(depths[i]):
        done = true

proc main(cutoff:int, chrom:string="1", start:int=0, stop:int=249250621, beds : seq[string]) =
  var gens = new_seq[Gen](len(beds))
  for i, b in beds:
    gens[i] = gen_from_bed(b, chrom, start, stop)
  var sum: int
  for value in meets(cutoff, gens):
    sum += value.stop - value.start
  echo sum

when isMainModule:
  var beds = new_seq[string]()
  for s in commandLineParams():
    beds.add(s)
  main(10, beds=beds)
