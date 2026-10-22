"""Exact polynomial comparison of independently generated kernel tensors."""
import json,itertools,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent
DEFAULT=HERE.parent/"session30-decorated-engine"/"symbolic-kernels.json"
from collections import defaultdict
from fractions import Fraction
from decorated_engine import symbolic_kernel
raw=json.load(open(sys.argv[1] if len(sys.argv)>1 else DEFAULT));out=[]
for r in raw:
 ed=tuple(tuple(x)for x in r['edges']);aut=r['automorphisms'];sy=symbolic_kernel(ed,aut)
 ours=defaultdict(Fraction)
 for key,c in sy.items():
  opts=[[(0,-1),(1,-1),(2,-1)] if k==3 else [(3 if k==4 else k,1)] for k in key]
  for vals in itertools.product(*opts):
   kk=tuple(t[0]for t in vals);cc=c
   for _,s in vals:cc*=s
   ours[kk]+=cc
 ours={k:v for k,v in ours.items()if v}
 theirs={tuple(t['indices']):Fraction(t['numerator'],r['denominator']) for t in r['terms'] if t['numerator']}
 if ours!=theirs:
  diff={k:ours.get(k,0)-theirs.get(k,0)for k in ours.keys()|theirs.keys() if ours.get(k,0)!=theirs.get(k,0)}
  print('MISMATCH',ed,len(diff),list(diff.items())[:4],flush=True)
  raise RuntimeError('exact tensor mismatch')
 rec={'edges':ed,'automorphisms':aut,'weight':r['weight'],'nonzero_tensor_entries':len(ours),'exact_equal':True};out.append(rec);print('EXACT',r['weight'],ed,len(ours),flush=True)
json.dump(out,open(HERE/'symbolic_comparison.json','w'),indent=2)
print('ALL31 EXACTLY EQUAL')
