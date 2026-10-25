import itertools,json,sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from edge_actions import relabel,pairs,require

def rank(a):
 a=[[F(x) for x in row] for row in a];r=0
 for j in range(len(a[0])):
  z=next((i for i in range(r,len(a)) if a[i][j]),None)
  if z is None:continue
  a[r],a[z]=a[z],a[r];v=a[r][j];a[r]=[x/v for x in a[r]]
  for i in range(len(a)):
   if i!=r:
    v=a[i][j];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
  r+=1
 return r

def main():
 n=6;N=1<<15;full=N-1;perms=list(itertools.permutations(range(n)));gc={};tc={};reps=[];sizes={}
 for x in range(N):
  if x not in gc:
   orb={relabel(x,n,p) for p in perms};r=min(orb);gc.update(dict.fromkeys(orb,r))
  if x not in tc:
   orb={relabel(x,n,p,True) for p in perms};orb|={y^full for y in orb};r=min(orb);tc.update(dict.fromkeys(orb,r));reps.append(r);sizes[r]=len(orb)
 hist={r:Counter(min(gc[y],gc[y^full]) for y in [relabel(r,n,p,True) for p in perms]) for r in reps}
 cols=sorted(set().union(*hist.values()));A=[[hist[r].get(c,0) for c in cols] for r in reps]
 require(len(reps)==34 and len({tuple(row) for row in A})==34,'direct classes/profiles')
 require(rank(A)==32,'direct profile matrix rank')
 # S_X for zero-sum graph mask255, independently evaluated on representatives.
 X=255;terms={}
 for p in perms:
  y=relabel(X,n,p);sign=(-1)**sum(p[i]>p[j] for k,(i,j) in enumerate(pairs(n)) if X>>k&1)
  require(y not in terms or terms[y]==sign,'coherent signed coefficients');terms[y]=sign
 f={r:sum(s*(-1)**((r&y).bit_count()%2) for y,s in terms.items()) for r in reps}
 w={r:sizes[r]*f[r] for r in reps};require(sum(w.values())==0,'kernel total mass')
 for c in cols:require(sum(w[r]*hist[r].get(c,0) for r in reps)==0,'mixture equality')
 mass=sum(v for v in w.values() if v>0);require(mass>0,'nonzero signed kernel')
 out={'n':n,'graph_pair_columns':len(cols),'converse_classes':len(reps),'distinct_profiles':len({tuple(row) for row in A}),'matrix_rank':rank(A),'columns':cols,'representatives':reps,'matrix':A,'signed_kernel':w,'positive_mass':mass,'mixture_positive':{r:str(F(v,mass)) for r,v in w.items() if v>0},'mixture_negative':{r:str(F(-v,mass)) for r,v in w.items() if v<0}}
 Path(__file__).with_name('direct_certificate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({k:v for k,v in out.items() if k not in ['matrix','columns']},indent=2))
if __name__=='__main__':main()
