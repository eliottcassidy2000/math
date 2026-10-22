"""Exact parity-basis evaluator with per-edge palette compression.
Source metadata and denominators are retained; no frozen source is modified.
"""
import pathlib,sys,json,copy,functools
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'session30-decorated-engine'))
from decorated import compiled_reduced
CORES=json.loads((ROOT/'parity-kernels.json').read_text())
def basis(M,z=1.):
 X=z*M;I=np.eye(len(M));Y=X@X;A=np.linalg.inv(I-Y);C=np.linalg.inv(I+Y)
 return np.array([-Y@C,Y@A,X@A,X@C])
@functools.lru_cache(None)
def compiled(weight):
 out=[]
 for c,(_,old) in zip([c for c in CORES if c['weight']==weight],compiled_reduced(weight)):
  ob=copy.copy(old);W=np.zeros((4,)*len(c['edges']))
  for t in c['terms']:W[tuple(t['indices'])]=t['numerator']
  ob.Wb=W;ob.pathb=None
  palettes=[sorted({t['indices'][i] for t in c['terms']}) for i in range(len(c['edges']))]
  out.append((c,ob,palettes,W[np.ix_(*palettes)]))
 return out
PATHS={}
def core_value(c,ob,palettes,W,B):
 # Preserve the old optimized K4/doubled-square contractions when all edges use four channels.
 if all(len(p)==4 for p in palettes):return ob.value_b(B)
 ops=[np.diagonal(B[p],axis1=1,axis2=2) if a==b else B[p] for (a,b),p in zip(c['edges'],palettes)]
 key=(c['weight'],tuple(map(tuple,c['edges'])),len(B[0]))
 if key not in PATHS:PATHS[key]=np.einsum_path(ob.spec,W,*ops,optimize=('greedy', 2**24))[0]
 return np.einsum(ob.spec,W,*ops,optimize=PATHS[key])/c['denominator']
def evaluate(M,weight,z=1.,detail=False):
 B=basis(M,z);out=[(c,core_value(c,ob,p,W,B)) for c,ob,p,W in compiled(weight)]
 return out if detail else sum(v for c,v in out)
