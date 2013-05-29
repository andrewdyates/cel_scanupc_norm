#!/usr/bin/python
"""
python ~/source_code/cel_scanupc_norm/script.py fdir=$HOME/brca/GSE31448/raw/cel outdir=$HOME/brca/GSE31448/normed dry=True
python ~/source_code/cel_scanupc_norm/script.py fdir=$HOME/brca/GSE31448/raw/cel.bkp outdir=$HOME/brca/GSE31448/normed gse=GSE31448 dry=True 
"""
from lab_util import *
import qsub
import sys, os, shutil, re
LOCALDIR = os.path.abspath(os.path.dirname(__file__))

PLATFORMS = set(["hgu133plus2hsentrezg"])
RX_BATCHDIR = re.compile("batch_\d+")
R_NORM_TMP = open(os.path.join(LOCALDIR,"R_norm.R.tmp")).read()
R_COMPILE_TMP = open(os.path.join(LOCALDIR,"R_compile.R.tmp")).read()


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
  if i==0:
    print "WARNING: no CEL files found. Maybe the path %s or the file pattern %s are wrong?" % (fdir, ptn)
    return None
  print "Sorted %d files into %d directories of <=%d files each" % (i, len(members), n)
  
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
  

      
def main(fdir=None, n=50, ptn=".CEL.gz", outdir=None, dosplit=True, platform="hgu133plus2hsentrezg", dry=False, gse=None):
  assert fdir
  assert n
  assert ptn
  assert outdir
  assert gse
  n = int(n)
  if isinstance(dosplit, basestring) and dosplit.lower() in ('f','false','none'): dosplit = False
  if isinstance(dry, basestring) and dry.lower() in ('f','false','none'): dry = False
    
  if dosplit:
    members = split_cels(fdir, n, ptn, dry)
  else:
    members = read_split(fdir, ptn)
  if members is None:
    print "Exiting..."
    sys.exit(1)

  if not os.path.exists(outdir):
    print "outdir %s does not exist. creating..." % outdir
    make_dir(outdir)
    
  scan_outfiles = set(); upc_outfiles = set()
  pids = set()
  for mdir, mems in members.items():
    batchname = os.path.basename(mdir)
    scan_outfile = "%s/%s.%s.RData" % (outdir,"SCAN",batchname)
    upc_outfile = "%s/%s.%s.RData" % (outdir,"UPC",batchname)
    script_SCAN = R_NORM_TMP % {'name':platform, "path":mdir, "ptn":ptn, "cmd":"SCAN", "batch":batchname, "outfile":scan_outfile}
    script_UPC = R_NORM_TMP % {'name':platform, "path":mdir, "ptn":ptn, "cmd":"UPC", "batch":batchname, "outfile":upc_outfile}
    scan_outfiles.add(scan_outfile)
    upc_outfiles.add(upc_outfile)
    
    scan_fpath = os.path.join(outdir,"%s.%s.R"%(batchname,"SCAN"))
    upc_fpath = os.path.join(outdir,"%s.%s.R"%(batchname,"UPC"))
    print "Writing %s, %s..." % (scan_fpath,upc_fpath)
    open(scan_fpath,"w").write(script_SCAN)
    open(upc_fpath,"w").write(script_UPC)
    # submit scripts to qsub
    Q1 = qsub.Qsub(n_nodes=1, n_ppn=8, hours=6, email=True, jobname=gse+scan_fpath)
    Q1.add("R CMD BATCH %s %s.Rout" % (scan_fpath, scan_fpath))
    pids.add(Q1.submit(dry))
    Q2 = qsub.Qsub(n_nodes=1, n_ppn=8, hours=6, email=True, jobname=gse+upc_fpath)
    Q2.add("R CMD BATCH %s %s.Rout" % (upc_fpath, upc_fpath))
    pids.add(Q2.submit(dry))
      
  # create compile script
  load_cmds = make_load_cmds(scan_outfiles, upc_outfiles)
  scan_expr_list = make_expr_list(scan_outfiles)
  upc_expr_list = make_expr_list(upc_outfiles)
  script_compile = R_COMPILE_TMP % {'load_cmds': load_cmds, 'gse':gse, 'scan_expr_list':scan_expr_list, 'upc_expr_list':upc_expr_list, 'outdir':outdir}
  compile_fname = os.path.join(outdir,"compile.R")
  print "Writing compile script %s..." % compile_fname
  open(compile_fname,"w").write(script_compile)
  Q3 = qsub.Qsub(hours=1, email=True, jobname="COMPILE", after_jobids=pids)
  Q3.add("R CMD BATCH %s %s.Rout" % (compile_fname, compile_fname))
  print Q3.script()
  print Q3.submit(dry)
    

def make_load_cmds(scan_outfiles, upc_outfiles):
  z = []
  for s in scan_outfiles:
    z.append('load("%s")'%os.path.abspath(s))
  for s in upc_outfiles:
    z.append('load("%s")'%os.path.abspath(s))
  return "\n".join(z)

def make_expr_list(outfiles):
  z = []
  for s in outfiles:
    z.append("exprs(E."+os.path.basename(s).replace('.RData','')+")")
  return ", ".join(z)

if __name__ == "__main__":
  args = dict((s.split('=') for s in sys.argv[1:]))
  print args
  main(**args)
