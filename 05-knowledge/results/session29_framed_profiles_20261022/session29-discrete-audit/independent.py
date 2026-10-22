from itertools import combinations,permutations,product
from collections import Counter
from math import prod
import json

def ck(x,s):
 if not x: raise RuntimeError(s)
def trajectory(n,h):
 out=[]; nodes=[n]
 for i in range(h):
  n=3*n+1;a=0
  while not n%2:a+=1;n//=2
  out.append(a);nodes.append(n)
 return tuple(out),nodes

def comps(a):
 if a==0:yield ();return
 for b in range(1,a+1):
  for tail in comps(a-b):yield (b,)+tail
checks=0
for A in range(1,10):
 M=2**(A+1)
 # Direct modular exponent enumeration, independent of the binary lift.
 exp_to_u={k:((pow(3,k,4*M)+1)//4)%M for k in range(1,M,2)}
 ck(len(set(exp_to_u.values()))==M//2,'oddpart bijection')
 for w in comps(A):
  actual=[n for n in range(1,M,2) if trajectory(n,len(w))[0]==w]
  ck(len(actual)==1,'one cylinder by exhaustive residues')
  r=actual[0];C=0;B=0
  for v in w:C=3*C+2**B;B+=v
  ck((3**len(w)*r+C)%M==2**A,'closed cylinder')
  ks=[k for k,u in exp_to_u.items() if u==r]
  ck(len(ks)==1,'one odd exponent class')
  # Exact word and all prefix clocks for a root-free lift.
  n=(3**(ks[0]+3*M)+1)//4
  vw,nodes=trajectory(n,len(w))
  ck(vw==w and 1 not in nodes,'root-free finite lift')
  registers={((3**(ks[0]+j*M)+1)//4)%105 for j in range(3)}
  ck(len(registers)==3,'exactly three sidecars')
  checks+=1
# Cardinality check over many integer annuli; proof still supplies all N.
u=[(3**k+1)//4 for k in range(1,31,2)]
ck(max(sum(4*N<=s<256*N for s in u) for N in range(1,10001))==2,'annulus cardinality')

def matrix(mask,n):
 a=[[0]*n for _ in range(n)]
 for b,(i,j) in enumerate(combinations(range(n),2)):
  a[i][j]=(mask>>b)&1;a[j][i]=1-a[i][j]
 return a

def encode(a,p):return sum(a[p[i]][p[j]]<<b for b,(i,j) in enumerate(combinations(range(len(p)),2)))
def canon(a,vs):return min(encode(a,p) for p in permutations(vs))
seeds=[]
for mask in [40,41]:
 a=matrix(mask,5);ps=list(permutations(range(5)))
 paths=[p for p in ps if all(a[x][y] for x,y in zip(p,p[1:]))]
 aut=[p for p in ps if encode(a,p)==mask]
 anti=[p for p in ps if encode(a,p)==1023-mask]
 canonical_path=lambda p:min(tuple(g[x] for x in p) for g in aut)
 tiles={canonical_path(p) for p in paths}
 mirrors=[{p:canonical_path(tuple(g[x] for x in p[::-1])) for p in tiles} for g in anti]
 ck(all(m==mirrors[0] for m in mirrors),'anti choice full map equality')
 ck(all(mirrors[0][mirrors[0][p]]==p for p in tiles),'mirror involution')
 fixed=sum(mirrors[0][p]==p for p in tiles)
 decks=[dict(Counter(canon(a,vs) for vs in combinations(range(5),k))) for k in range(1,5)]
 seeds.append({'mask':mask,'H':len(paths),'aut':len(aut),'anti':len(anti),'tiles':len(tiles),'fixed':fixed,'decks':decks})
ck(seeds[0]['decks']==seeds[1]['decks'],'all proper decks')
ck([(s['H'],s['aut'],s['tiles'],s['fixed']) for s in seeds]==[(15,3,5,1),(13,1,13,3)],'seed totals')
coeffs=[]
for r in range(11):
 words=list(product(range(2),repeat=r));orbits={min(w,w[::-1]) for w in words}
 row={'r':r,'nodes':len(orbits)}
 ck(row['nodes']==(2**r+2**((r+1)//2))//2,'Burnside nodes')
 for ws in [(15,13),(5,13)]:
  mass=lambda w:prod(ws[i] for i in w)
  one=sum(mass(w) for w in orbits)
  formula=(sum(ws)**r+sum(t*t for t in ws)**(r//2)*(sum(ws) if r%2 else 1))//2
  ck(one==formula,'weighted coefficient')
 actualmass=sum(prod((5,13)[i] for i in w)*(1 if w==w[::-1] else 2) for w in orbits)
 blue=sum(prod((5,13)[i] for i in w[:r//2])*((1,3)[w[r//2]] if r%2 else 1) for w in orbits if w==w[::-1])
 ck(actualmass==18**r and blue==18**(r//2)*(4 if r%2 else 1),'merged mass blue')
 row.update(mass=actualmass,blue=blue);coeffs.append(row)
# Strict uncrossing local binary inequalities, independent of four printed formulas.
local_cases=0
for zinA,targetB in product([False,True],repeat=2):
 for p,r,q,s in product([0,1],repeat=4):
  if zinA and r!=0 or not zinA and p!=1:continue
  if targetB and q!=0 or not targetB and s!=1:continue
  ck(abs(p-s)+abs(r-q)<=abs(p-q)+abs(r-s),'uncrossing pair contribution')
  local_cases+=1
print(json.dumps({'status':'PASS','valuation_words':checks,'annulus_N_max':10000,'seeds':seeds,'coefficients':coeffs,'uncrossing_local_cases':local_cases},indent=2))
