"""Independent mod-eight profile and ambient audit."""
from itertools import combinations,permutations
from math import comb,factorial
from functools import lru_cache
from pathlib import Path
import json,random,hashlib

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def rows(n,code):
 a=[0]*n
 for bit,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(j,i) if code>>bit&1 else(i,j);a[u]|=1<<v
 return tuple(a)

def induced(a,mask):
 names=[i for i in range(len(a)) if mask>>i&1]
 return tuple(sum(1<<j for j,v in enumerate(names) if a[u]>>v&1) for u in names)

@lru_cache(None)
def q_words(a):
 n=len(a)
 if not n:return (1,)
 q=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(not(a[u]>>v&1) for u,v in zip(w,w[1:]));good=n-1-bad
  for k in range(bad+1,n+1):q[k]+=comb(good,k-1-bad)
 return tuple(q)

@lru_cache(None)
def hp(a):
 n=len(a);full=(1<<n)-1
 @lru_cache(None)
 def finish(mask,last):
  if mask==full:return 1
  out=0;rem=a[last]&~mask
  while rem:
   b=rem&-rem;rem-=b;out+=finish(mask|b,b.bit_length()-1)
  return out
 return sum(finish(1<<i,i) for i in range(n)) if n else 1

def blow(Q,blocks):
 offsets=[0]
 for b in blocks:offsets.append(offsets[-1]+len(b))
 a=[0]*offsets[-1]
 for i,b in enumerate(blocks):
  for u,r in enumerate(b):
   a[offsets[i]+u]=r<<offsets[i]
   for j,c in enumerate(blocks):
    if Q[i]>>j&1:a[offsets[i]+u]|=((1<<len(c))-1)<<offsets[j]
 return tuple(a)

def quotient(a,blocks):
 reps=[(M&-M).bit_length()-1 for M in blocks]
 return tuple(sum(1<<j for j,v in enumerate(reps) if a[u]>>v&1) for u in reps)

def extend(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for v in range(n):
  if x>>v&1:b[v]|=1<<n
  else:b[n]|=1<<v
  if c>>v&1:b[v]|=1<<(n+1)
  else:b[n+1]|=1<<v
 return tuple(b)


def cycles(a):
 n=len(a)
 if n<3:return 0
 return sum(all(a[u]>>v&1 for u,v in zip(w,w[1:]+w[:1]))
            for tail in permutations(range(1,n)) for w in [(0,)+tail])

def stirling(n,k):
 if n==k==0:return 1
 if not n or not k:return 0
 return stirling(n-1,k-1)+k*stirling(n-1,k)

@lru_cache(None)
def q_partitions(a):
 n=len(a);H={S:hp(induced(a,S)) for S in range(1,1<<n)}
 @lru_cache(None)
 def count(S,k):
  if not S:return int(k==0)
  if not k:return 0
  bit=S&-S;rest=S^bit;part=rest;total=0
  while True:
   B=part|bit;total+=H[B]*count(S^B,k-1)
   if not part:break
   part=(part-1)&rest
  return total
 return (0,)+tuple(factorial(k)*count((1<<n)-1,k) for k in range(1,n+1))

def verify_profile(a,literal):
 n=len(a);q=q_words(a) if literal else q_partitions(a);H=q[1];z=cycles(a) if n%2 else 0
 at=lambda k:q[k] if k<len(q) else 0
 need((at(2)-2*H-(1<<n)+4+4*z)%8==0,'q2 mod8')
 need((at(3)-6*stirling(n,3)-(2*H-2-4*z if n%2 else 0))%8==0,'q3 mod8')
 need(all(v%8==0 for v in q[4:]),'higher cover divisibility')
 if n<=5:need(q==q_partitions(a),'literal word/canonical partition agreement')
 return q,z

def put(Q,r,A):
 blocks=tuple(A if i==r else(0,) for i in range(len(Q)))
 return blow(Q,blocks)

def main():
 profiles=0;cycles_checked=0;rng=random.Random(188127)
 for n in range(1,9):
  count=(1<<(n*(n-1)//2)) if n<=5 else {6:32,7:8,8:4}[n]
  codes=range(count) if n<=5 else [rng.randrange(1<<(n*(n-1)//2)) for _ in range(count)]
  for code in codes:
   a=rows(n,code);verify_profile(a,True);profiles+=1;cycles_checked+=n%2
 ambient=0;residue_counts={str(k):0 for k in(0,2,4,6)}
 A=rows(4,0)
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   Q=rows(n,code)
   for r in range(n):
    factor=hp(put(Q,r,rows(2,0)));need(factor%2==1,'ambient multiplier not odd')
    for b in range(64):
     B=rows(4,b);delta=hp(B)-hp(A)
     actual=hp(put(Q,r,B))-hp(put(Q,r,A))
     need((actual-factor*delta)%8==0,'literal ambient odd-unit law')
     need((actual%8==0)==(delta%8==0),'mod8 zero not preserved')
     ambient+=1;residue_counts[str(delta%8)]+=1
 # Use exact physical rows independently retained before the producer packet.
 controls=json.loads((Path(__file__).resolve().parent/'physical_controls.json').read_text())
 checked=[]
 for r in controls:
  X,Y=tuple(r['port_core_T']),tuple(r['port_core_U']);qx,zx=verify_profile(X,False);qy,zy=verify_profile(Y,False)
  need(list(qx[1:])==r['q_T'] and list(qy[1:])==r['q_U'],'stored physical core profiles')
  D=[a-b for a,b in zip(qx,qy)]
  gap=hp(put((2,4,1),0,X))-hp(put((2,4,1),0,Y))
  need(gap==3*D[1]+D[2],'triangle context coefficient identity')
  if len(X)%2:need(gap%8==0 and gap!=0 and (zx-zy)%2==1,'odd boundary')
  else:need(gap%8==4,'even physical lift')
  checked.append({'core_order':len(X),'q_difference':D,'cycle_difference':zx-zy,'ambient_gap':gap})
 root=Path(__file__).resolve().parent.parent/'mod8_transfer'
 pins={f:hashlib.sha256((root/f).read_bytes()).hexdigest() for f in ['PROOF.md','verify.py','receipt.json','receipt.optimized.json']}
 return {'status':'PASS','independent_profile_cases':profiles,'literal_cycle_cases':cycles_checked,
         'literal_ambient_comparisons':ambient,'core_difference_residue_counts':residue_counts,'physical_controls':checked,'producer_pins':pins}

if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
