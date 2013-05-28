#!/usr/bin/python
"""
python ~/source_code/cel_scanupc_norm/scratch.py fdir=$HOME/brca/GSE31448/raw/cel outdir=$HOME/brca/GSE31448/normed dry=True 
"""
from lab_util import *
import qsub
import sys, os, shutil, re

PLATFORMS = set(["hgu133plus2hsentrezg"])
RX_BATCHDIR = re.compile("batch_\d+")
R_SCRIPT_TMP = """library(SCAN.UPC)
library(%(name)sprobe)
library(%(name)s.db)

celFilePath <- "%(path)s/*%(ptn)s"
norm.cust = %(cmd)s(celFilePath, probeSummaryPackage=hgu133plus2hsentrezgprobe)
IDS <- rownames(exprs(norm.cust)) # note: IDS are entrez_ids appended by _at
EIDS <- mget(IDS, %(name)sENTREZID, ifnotfound=NA)
rm.noann <- is.na(EIDS)
E.%(cmd)s.%(batch)s <- norm.cust[!rm.noann,]
rownames(E.%(cmd)s.%(batch)s) <- EIDS[!rm.noann]
# every row in EID is now uniquely identified by an entrez ID
save(E.%(cmd)s.%(batch)s, file="%(outdir)s/%(cmd)s.%(batch)s.RData")
"""

def split_cels(fdir, n=50, ptn=".CEL.gz", dry=False):
  """Given a directory of .CEL.gz files, split into subdirectories of <=n files each."""
  c = 1; i = 0
  dirpath = os.path.join(fdir,"batch."+str(c))
  members = {dirpath:set()}; 
  for fname in os.listdir(fdir):
    if ptn not in fname: continue
    fpath = os.path.join(fdir, fname)
    if len(members[dirpath]) < n:
      members[dirpath].add(fpath)
      i += 1
    else:
      c += 1; i += 1
      dirpath = os.path.join(fdir,"batch."+str(c))
      members[dirpath] = set([fpath])
  print "Created %d directories of %d files, <=%d each" % (len(members), i, n)
  
  for cdir, mems in members.items():
    print "Making directory %s..." %cdir
    make_dir(cdir)
    for s in mems:
      if not dry:
        shutil.move(s,cdir)
  return members


def read_split(fdir, ptn):
  i = 0
  for fname in os.listdir(fdir):
    if RX_BATCHDIR.match(fname):
      fpath = os.path.join(fdir,fname)
      members[fpath] = set()
      for s in os.listdir(fpath):
        if ptn in s:
          i += 1
          members[fpath].add(os.path.join(fpath,s))
  print "Read %d directories of %d files." % (len(members), i)
  

      
def main(fdir=None, n=50, ptn=".CEL.gz", outdir=None, dosplit=True, platform="hgu133plus2hsentrezg", dry=False):
  assert fdir
  assert n
  assert ptn
  assert outdir
  n = int(n)
  if isinstance(dosplit, basestring) and dosplit.lower() in ('f','false','none'): dosplit = False
  if isinstance(dry, basestring) and dry.lower() in ('f','false','none'): dry = False
    
  if dosplit:
    members = split_cels(fdir, n, ptn, dry)
  else:
    members = read_split(fdir, ptn)

  if not os.path.exists(outdir):
    print "outdir %s does not exist. creating..." % outdir
    make_dir(outdir)
  
  for mdir, mems in members.items():
    batchname = os.path.basename(mdir)
    script_SCAN = R_SCRIPT_TMP % {'name':platform, "path":mdir, "ptn":ptn, "cmd":"SCAN", "batch":batchname,"outdir":outdir}
    script_UPC = R_SCRIPT_TMP % {'name':platform, "path":mdir, "ptn":ptn, "cmd":"UPC", "batch":batchname,"outdir":outdir}
    print script_SCAN
    print script_UPC
    scan_fpath = os.path.join(outdir,"%s.%s.R"%(batchname,"SCAN"))
    upc_fpath = os.path.join(outdir,"%s.%s.R"%(batchname,"UPC"))
    print "Writing %s, %s..." % (scan_fpath,upc_fpath)
    open(scan_fpath,"w").write(script_SCAN)
    open(upc_fpath,"w").write(script_UPC)
    # submit scripts to qsub
    if not dry:
      Q1 = qsub.Qsub(n_nodes=1, n_ppn=8, hours=99, email=True, jobname=scan_fpath)
      Q1.add("R CMD BATCH %s %s.Rout" % (scan_fpath, scan_fpath))
      Q1.submit(dry)
      Q2 = qsub.Qsub(n_nodes=1, n_ppn=8, hours=99, email=True, jobname=upc_fpath)
      Q2.add("R CMD BATCH %s %s.Rout" % (upc_fpath, upc_fpath))
      Q2.submit(dry)
      
    

if __name__ == "__main__":
  args = dict((s.split('=') for s in sys.argv[1:]))
  print args
  main(**args)
