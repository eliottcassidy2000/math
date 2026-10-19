"""Exact lone3 risk via proved binary-tail shield, parallel bitplane CA."""
import numpy as np
from fractions import Fraction
import sys,json,time

def scan(F):
 t0=time.time();n=1<<(F-1);mask=(1<<F)-1;q=1<<(F-1).bit_length()
 b=np.arange(n,dtype=np.uint64)|(1<<(F-1));h=np.full(n,1<<(F-1),dtype=np.uint64)
 fail=np.full(n,-1,dtype=np.int32)
 checkpoints={};repeated=False
 for t in range(4*q+1):
  hit=(fail<0)&((h&1)!=0);fail[hit]=t
  if t%q==0:
   if t>=q:
    oldb,oldh=checkpoints['last'];live=(fail<0)
    repeated=bool(np.all(b[live]==oldb[live])and np.all(h[live]==oldh[live]))
    if repeated:break
   checkpoints['last']=(b.copy(),h.copy())
  c=b>>1;k=h>>1
  h=((h&~k&(b|~c)) | (k&~h&(c|~b)))&mask
  b=(b^c)&mask
 if not repeated:raise RuntimeError(f'No complete cycle certificate F={F}')
 # Time t first hits firstcolumn. It is apex of initialprefix of lengtht+1.
 # If t<F, failure already inside w3, else needs t-F+1 zero sea bits.
 exponents=np.maximum(0,fail[fail>=0]-(F-1));bins=np.bincount(exponents)
 val=sum((int(ct)*Fraction(1,1<<i)for i,ct in enumerate(bins)),Fraction(0))/n
 excess=val*2**F-2*F
 print(json.dumps({'F':F,'initial_words':n,'zero_input_rows_checked':t,'live_states_period_divides':q,'bad_prefix_count':int(bins[0]),'threshold_counts':{i:int(c)for i,c in enumerate(bins)if c},'max_threshold':len(bins)-1,'p':str(val),'scaled_excess':str(excess),'float_excess':float(excess),'seconds':round(time.time()-t0,3)}),flush=True)
 return val
if __name__=='__main__':
 for F in range(int(sys.argv[1]),int(sys.argv[2])+1):scan(F)
