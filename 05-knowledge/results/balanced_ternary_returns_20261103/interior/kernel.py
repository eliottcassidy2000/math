from pathlib import Path
import importlib.util,itertools,json
ROOT=Path(__file__).resolve().parents[4]
p=ROOT/'05-knowledge/results/orbit_spectrum_session_20261103/coverage/check.py'
s=importlib.util.spec_from_file_location('old',p);c=importlib.util.module_from_spec(s);s.loader.exec_module(c)

def cut_kernel(a):
 n=len(a);full=(1<<n)-1;_,ends=c.tables(a)
 rev=tuple(sum(1<<j for j in range(n) if a[j]>>i&1) for i in range(n));_,starts=c.tables(rev)
 W=[[0]*n for _ in range(n)]
 for mask in range(1,full):
  comp=full^mask
  for u in range(n):
   if ends[mask][u]:
    for v in range(n):
     W[u][v]+=ends[mask][u]*starts[comp][v]
 return ends[full],starts[full],W

def cut_h(kernel,imask):
 ends,starts,W=kernel;n=len(ends)
 I=[i for i in range(n) if imask>>i&1];O=[i for i in range(n) if not imask>>i&1]
 return sum(ends[i] for i in I)+sum(starts[j] for j in O)+sum(W[i][j] for i in I for j in O)

def add(a,imask):
 n=len(a);b=[a[i]|((1<<n) if imask>>i&1 else 0) for i in range(n)]
 return tuple(b+[((1<<n)-1)^imask])

if __name__=='__main__':
 for w in ['202032','2020032','20200032','20203332','202030332']:
  a=c.graph(w);h=c.word.hamilton(w);kernel=cut_kernel(a);hits={b:[] for b in (2,4)}
  for imask in range(1<<len(a)):
   x=cut_h(kernel,imask)
   for b in (2,4):
    if x==3*h+b:hits[b].append(imask)
  print(w,h,hits,flush=True)
