from itertools import product
from functools import lru_cache
from pathlib import Path
import json,hashlib,importlib.util
P={'N': ((6, 4, 10, 9), (4, 3, 7, 6), (10, 7, 15, 14), (9, 6, 14, 13)), 'M': ((11, 8, 20, 18), (8, 6, 14, 13), (20, 14, 34, 31), (18, 13, 31, 29)), 'A': ((14, 10, 21, 20), (10, 7, 16, 15), (21, 16, 33, 31), (20, 15, 31, 29)), 'B': ((14, 10, 26, 23), (11, 8, 20, 18), (20, 15, 35, 32), (19, 14, 34, 31)), 'C': ((23, 18, 44, 40), (18, 14, 34, 31), (44, 34, 80, 73), (40, 31, 73, 67)), 'T': ((14, 11, 20, 19), (10, 8, 15, 14), (26, 20, 35, 34), (23, 18, 32, 31)), 'J0': ((12, 12, 24, 24), (12, 12, 21, 21), (15, 12, 30, 27), (15, 12, 30, 27)), 'J1': ((12, 12, 21, 21), (12, 12, 18, 18), (18, 15, 32, 30), (18, 15, 31, 29)), 'J2': ((12, 12, 18, 18), (12, 12, 15, 15), (21, 18, 32, 31), (21, 18, 30, 29)), 'J3': ((12, 12, 15, 15), (12, 12, 12, 12), (24, 21, 30, 30), (24, 21, 27, 27))}
ROOT=Path(__file__).resolve().parents[1]

def need(c,m):
 if not c:raise RuntimeError(m)
def step(v,l):
 D,E,U,d=v
 return (D,0,2*U,d) if l=='O' else ((D-E+U,d,2*U-d,d) if l=='Z' else (D-E+U,D-E,D+U,D))
def run(v,w):
 for l in w:v=step(v,l)
 return v
def y(v):
 D,E,U,d=v;return(D-E-d,U-D+E-d,E,d-E)
def dot(v,w):return sum(a*b for a,b in zip(v,w))
def qrow(P,v):
 z=y(v);r=[dot(row,z) for row in P]
 return(r[0]-r[1],-r[0]+r[1]+r[2]-r[3],r[1],-r[0]-r[1]+r[3])
def F(v):return v[0]-v[1]+v[2]
@lru_cache(None)
def g(n):
 if n==1:return 1
 if n%2==0:return g(n//2)
 h=(n+1)//2;return g(h)+h-h%2
seed=(1,0,2,1)
states={run(seed,w) for k in range(4) for w in product('OZJ',repeat=k)}
partners={run(seed,w) for k in range(3) for w in product('OZJ',repeat=k)}
rows={(1,-1,1,0),(5,-2,16,-4),(3,-3,9,-6)}
for A in P.values():
 for B in [A,tuple(zip(*A))]:
  for p in partners:rows.add(qrow(B,p))
exchanges=0
G=[((1,0,1,0),(1,2,1,2),(0,0,0,0),(0,0,1,1)),((2,1,1,1),(0,1,0,0),(0,0,1,1),(0,0,0,0)),((0,1,0,1),(0,0,0,0),(1,0,1,1),(0,0,1,0))]
cone_checks=0
for alpha,beta,gamma,eta in rows:
 r=(alpha+gamma,gamma,2*alpha+beta+2*gamma+eta,alpha+2*gamma+eta)
 current={r}
 for depth in range(6):
  for v in current:
   a,b,c,d=v
   need(c>a>b>0 and d>0 and d<a+c,'strict arbitrary suffix cone')
   cone_checks+=1
  current={tuple(sum(v[i]*A[i][j] for i in range(4)) for j in range(4)) for v in current for A in G}
for r in rows:
 al,be,ga,et=r;need(al>0 and ga>0 and 2*al+be+ga>0,'readout positivity')
 for s in states:
  for k in range(4):
   for w in product('OZJ',repeat=k):
    large=dot(r,run(s,('O','Z')+w));small=dot(r,run(s,('Z','O')+w))
    need(large>small,'strict suffix exchange');exchanges+=1
minmax={};cutrows={}
for k in range(1,13):
 counts=[];cutmin=None;cutargs=[]
 for word in product('OZ',repeat=k):
  counts.append(F(run(seed,word)))
  bits=[int(x=='Z') for x in word];A=sum(e*2**i for i,e in enumerate(bits));B=2**k
  energy=sum((2**(j-i-1)+int(j==i+1))*abs(bits[i]-bits[j]) for i in range(k) for j in range(i+1,k))
  need(energy==g(A+1)+g(B-A)-2,'independent cut identity')
  if 0<A<B-1:
   if cutmin is None or energy<cutmin:cutmin=energy;cutargs=[A]
   elif energy==cutmin:cutargs.append(A)
 if k>=2:
  need(min(counts)==(2**(k//2)+1)*(2**((k+1)//2)+1),'minimum formula')
  need(max(counts)==(2**(k//2+1)-1)*(2**((k+1)//2+1)-1)+2,'maximum formula')
  delta=2 if k==2 else 2**((k-1)//2)+2**(k//2)
  need(cutmin==delta,'cut minimum')
 if k>=4:
  c={1<<((k-1)//2),1<<(k//2)};c|={((1<<k)-1)^x for x in c};need(set(cutargs)==c,'cut equality')
 minmax[k]=[min(counts),max(counts)];cutrows[k]=cutmin
report={'status':'PASS','proof_sha256':hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),'script_sha256':hashlib.sha256((ROOT/'verify.py').read_bytes()).hexdigest(),'independent_fixed_partner_readout_rows':len(rows),'physical_initial_profiles':len(states),'strict_exchange_checks':exchanges,'arbitrary_OZJ_suffix_cone_checks':cone_checks,'block_and_cut_exhaustion_through':12,'minmax':minmax,'cutmin':cutrows,'manual_audit':'Cut minimum proof, equality boundary, signed versus positive interaction, arbitrary O/Z/J suffix exchange, sorted extrema, and both plus/Mersenne product envelopes checked. No mathematical correction required.','scope':'Extrema are not exact-level coverage.'}
Path(__file__).with_name('INDEPENDENT.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
