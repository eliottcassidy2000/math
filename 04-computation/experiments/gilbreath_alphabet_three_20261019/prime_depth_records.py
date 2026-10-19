"""Exact integer CA scan; NumPy accelerates but does not approximate."""
import argparse, json, time
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=200_000_000);a=p.parse_args();t=time.time()
sieve=np.ones((a.limit+1)//2,dtype=np.bool_);sieve[0]=False
for q in range(3,int(a.limit**.5)+1,2):
 if sieve[q//2]:sieve[q*q//2::q]=False
primes=np.concatenate((np.array([2],dtype=np.int32),2*np.flatnonzero(sieve).astype(np.int32)+1));primes=primes[primes<a.limit];del sieve
N=len(primes);row=np.diff(primes);gmax=np.maximum.accumulate(row);r=1;bounds=[1];details=[];leadok=True
while True:
 leadok=leadok and int(row[0])==1
 defect=row>2
 if not bool(np.any(defect)):
  bounds.append(N+1);depth=r;break
 i=int(np.argmax(defect));b=i+r+1;bounds.append(b)
 details.append({'r':r,'first_defect_index_0':i,'prefix_threshold':b,'value':int(row[i]),'source_start_prime':int(primes[i]),'source_end_prime':int(primes[i+r]),'source_max_gap':int(np.max(np.diff(primes[i:i+r+1]))),'prefix_max_gap':int(gmax[b-2])})
 row=np.abs(np.diff(row));r+=1
 if r>2000:raise RuntimeError('too deep')
records=[]
for k in range(1,depth+1):
 lo=max(3,bounds[k-1]);hi=min(N,bounds[k]-1)
 if lo>hi:continue
 records.append({'lo':lo,'hi':hi,'depth':k,'p_lo':int(primes[lo-1]),'p_hi':int(primes[hi-1]),'G_lo':int(gmax[lo-2]),'G_hi':int(gmax[hi-2]),'ratio_max':k/int(gmax[lo-2]),'ratio_min':k/int(gmax[hi-2])})
for cutoff in [3,100,1000,10000,100000,1000000,10000000]:
 cand=[dict(d,lo=max(d['lo'],cutoff))for d in records if d['hi']>=cutoff]
 if cand:
  for d in cand:d['ratio_max']=d['depth']/int(gmax[d['lo']-2])
  print('CUTOFF',cutoff,'MAX',max(cand,key=lambda d:d['ratio_max']),'MIN',min(cand,key=lambda d:d['ratio_min']),flush=True)
print('N',N,'pN',int(primes[-1]),'depth',depth,'maxgap',int(gmax[-1]),'leading_ones_checked',leadok,'seconds',round(time.time()-t,3),flush=True)
print('DEPTH_RECORDS',json.dumps(records),flush=True)
print('LAST_WITNESS',json.dumps(details[-1]),flush=True)
out=Path(__file__).parent/f'prime_depth_{a.limit}.json';out.write_text(json.dumps({'limit':a.limit,'N':N,'pN':int(primes[-1]),'depth':depth,'maxgap':int(gmax[-1]),'leading_ones_checked':leadok,'records':records,'first_defects':details},indent=2)+'\n')
