import hts
import strutils
import sequtils
import kexpr
import os
import docopt
import parsecsv

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
  template `+!`(src, a): untyped = cast[pointer](cast[int](cstring(src)) + a)
  setLen(dest, b-a)
  copyMem(addr dest[0], src+!a, b-a)

proc gen_from_bgzi(b:BGZI, iv:interval, depth_col:int=3): Gen =
  result = iterator(): interval {.inline.} =
    var col = depth_col
    var pos0, pos1: int
    var tmp = new_string_of_cap(20)
    var chrom = iv.chrom
    var start, stop, depth: int
    # a lot of extra stuff in here to minimize garbage creation and
    # copying
    for line in b.query(iv.chrom, iv.start, iv.stop):
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

iterator meets(depth_cutoff:int, expr:Expr, depths: seq[Gen], names: seq[string]): interval =
  var heap = new_seq[interval](depths.len)
  var maxi = new_seq[int](depths.len)
  for i in 0..<depths.len:
    heap[i] = depths[i]()
    if expr != nil:
      discard ke_set_int(expr.ke, names[i], cint(heap[i].depth))

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
    if expr != nil:
      if bool(expr):
        yield (chrom, istart, istop, 1)
    elif allok:
      yield (chrom, istart, istop, 1)
    # any site that ended at the max end
    # we replace with the next value.
    for i in maxi:
      heap[i] = depths[i]()
      if finished(depths[i]):
        done = true
      if expr != nil:
        discard ke_set_int(expr.ke, names[i], cint(heap[i].depth))

proc bed_line_to_region(line: string): interval =
  var
    cse = sequtils.to_seq(line.strip().split("\t"))

  if len(cse) < 3:
    stderr.write_line("[mosdepth] skipping bad bed line:", line.strip())
    return ("", 0, 0, 0)
  var
    s = parse_int(cse[1])
    e = parse_int(cse[2])
  return (cse[0], s, e, 0)

proc bed_to_table(bed: string): OrderedTableRef[string, seq[interval]] =
  ## this is used to read in a bed file of, e.g. capture regions into a table.
  var bed_regions = newOrderedTable[string, seq[interval]]()
  var hf = hts.hts_open(cstring(bed), "r")
  var kstr = hts.kstring_t(l:0, m:0, s:nil)
  while hts_getline(hf, cint(10), addr kstr) > 0:
    var l = $kstr.s
    if l.startswith("track "):
      continue
    if l.startswith("#"):
      continue
    var v = bed_line_to_region(l)
    if v.stop == 0: continue
    discard bed_regions.hasKeyOrPut(v.chrom, new_seq[interval]())
    bed_regions[v.chrom].add(v)

  hts.free(kstr.s)
  return bed_regions

proc region_main(cutoff:int, expr:string, bed_regions:string, depths: seq[BGZI], names: seq[string]) =
  var e : Expr
  if expr != "":
    e = expression(expr)
    if e.error() != 0:
      stderr.write_line("[mospow] error parsing expression:", expr)
      quit(e.error())

  var regions = bed_to_table(bed_regions)
  var gens = new_seq[Gen](len(depths))
  var sum: int
  for chrom, regs in regions:
    for reg in regs:
      for i, d in depths:
        gens[i] = d.gen_from_bgzi(reg)
      for iv in meets(cutoff, e, gens, names):
        sum += iv.stop - iv.start
      echo sum

proc isdigit(s:string): bool =
  for c in s:
    if not c.isdigit: return false
  return true


proc get_names(paths: seq[string], expr:string): seq[string] =
  if expr == "": return
  result = new_seq[string](len(paths))
  for i, p in paths:
    var tmp = p.split("/")
    var name = tmp[tmp.high]
    for suffix in @[".gz", ".bed", ".per-base"]:
      if name.ends_with(suffix):
        name = name[0..<(name.len - len(suffix))]
    if name.isdigit():
      name = "s" & name
    stderr.write_line("[mospow] using name:", name, " for:", p)
    result[i] = name

proc main(cutoff:int, expr:string, bed_regions:string, chrom: string, depth_beds: seq[string]) =
  var bgzis = new_seq[BGZI](len(depth_beds))
  for i, b in depth_beds:
    bgzis[i] = ropen_bgzi(b)

  if bed_regions != "":
    region_main(cutoff, expr, bed_regions, bgzis, get_names(depth_beds, expr))
    quit(0)

  # TODO:

when isMainModule:
  let doc = """
  mospow

  Usage:
    mospow [options] <DEPTH_BED>...

    -d --cutoff <cutoff>              depth required at every base in each input BED [default: 10]
    -e --expr <expr>                  expression on which to filter (use instead of cutoff)
    -r --region <region>              region to restrict power calculation (use instead of bed_regions)
    -b --bed-regions <bed_regions>    BED file of regions to restrict power calculation (TE, exons, etc.)

  """
  let args = docopt(doc, version = "mospow v0.0.2")
  var bed_regions: string = ""
  var expr = ""
  var region = ""
  # if a chrom other than 21 is specified, mospow won't know the length
  if $args["--region"] != "nil":
    region = $args["--region"]
  if $args["--expr"] != "nil":
    expr = $args["--expr"]

  var cutoff = parseInt($args["--cutoff"])

  if $args["--bed-regions"] != "nil":
    bed_regions = $args["--bed-regions"]

  main(cutoff, expr, bed_regions, region, @(args["<DEPTH_BED>"]))
