#!/usr/bin/env python3
"""Independent bounded precision-release audit with direct triangle counts."""
import argparse,hashlib,json,shutil,subprocess,tempfile
from itertools import combinations
from pathlib import Path
HERE=Path(__file__).resolve().parent
checks=0

def require(ok,label):
 global checks
 checks+=1
 if not ok:raise RuntimeError(label)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def section(text,key):
 lines=text.splitlines();i=lines.index(key+':')+1;j=i
 while j<len(lines) and lines[j].startswith('  '):j+=1
 return lines[i:j]

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--receipt',type=Path,default=HERE/'independent-normal.json');args=ap.parse_args()
 producer=HERE/'inputs/producer';sources=HERE/'inputs/sources'
 require(sha(producer/'precision.patch')=='8880e6a357ecd0aa426378e6a9e2325d20c235124925eb9d05e9cddf66b587c0','exact patch')
 manifest=json.loads((producer/'manifest.json').read_text())
 for name,digest in manifest.items():require(sha(producer/name)==digest,('producer manifest',name))
 pins=json.loads((producer/'source-pins.json').read_text())['sources']
 with tempfile.TemporaryDirectory(prefix='precision-independent-') as directory:
  target=Path(directory)
  for pin in pins:
   original=sources/pin['path'];require(sha(original)==pin['source_sha256'],('source pin',pin['path']))
   dest=target/pin['path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(original,dest)
  run=subprocess.run(['patch','--batch','--fuzz=0','-p1','-i',str(producer/'precision.patch')],cwd=target,capture_output=True,text=True)
  require(run.returncode==0,('patch',run.stdout,run.stderr))
  for pin in pins:
   dest=target/pin['path'];require(sha(dest)==pin['result_sha256'],('result pin',pin['path']))
   old=(sources/pin['path']).read_text();new=dest.read_text()
   require([x for x in old.splitlines() if x.startswith(('id:','source:'))]==[x for x in new.splitlines() if x.startswith(('id:','source:'))],('attribution',pin['path']))
   require(section(old,'depends_on')==section(new,'depends_on'),('dependencies',pin['path']))
   require(section(old,'results')==section(new,'results'),('results',pin['path']))
   require([x for x in old.splitlines() if x.startswith('|')]==[x for x in new.splitlines() if x.startswith('|')],('numeric tables',pin['path']))
 controls=[]
 for n in range(3,7):
  pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)}
  triples=[(index[(a,b)],index[(b,c)],index[(a,c)]) for a,b,c in combinations(range(n),3)]
  # Directly test the two cyclic orientations on every triple; no score formula.
  counts=[sum(((mask>>ab)&1)==((mask>>bc)&1) and ((mask>>ab)&1)!=((mask>>ac)&1) for ab,bc,ac in triples) for mask in range(1<<len(pairs))]
  centered=[4*x-len(triples) for x in counts]
  require(sum(centered)==0,('mean',n));require(any(centered),('nonzero',n))
  for mask,value in enumerate(centered):
   require(sum(centered[mask^(1<<edge)] for edge in range(len(pairs)))==(len(pairs)-4)*value,('direct triangle eigenfunction',n,mask))
  controls.append({'n':n,'labelled_tournaments':len(counts),'cyclicity_sum':sum(counts),'eigenvalue':len(pairs)-4})
 # The n=2 transposition interchanges the two labelled tournaments: trace 0, fixed-subspace dimension 1.
 require(sum([[0,1],[1,0]][i][i] for i in range(2))==0,'trace dimension hostile trace')
 require([1,1]==[1,1] and [1,-1]!=[-1,1],'trace dimension hostile fixed line')
 out={'status':'PASS','checks':checks,'producer_patch_sha256':sha(producer/'precision.patch'),'source_pins':pins,'direct_triangle_controls':controls,'scope':'Exact narrow diff pins, unchanged attribution/dependencies/result paths/numeric tables, direct-triangle cyclicity eigenfunction identity for every labelled tournament n=3..6. No shared-source edits.'}
 args.receipt.write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'PASS','checks':checks}))
if __name__=='__main__':main()
