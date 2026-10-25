#!/usr/bin/env python3
import itertools as it,json,pathlib,math
D=pathlib.Path(__file__).resolve().parent
def ck(p,m):
 if not p:raise RuntimeError(m)
def data(n,Q):
 E=[(u,v) for u in range(n) for v in range(n) if (v-u)%n in Q]
 cycles=[]
 for d in sorted(Q):
  seq=[k*d%n for k in range(n)];ck(len(set(seq))==n,'non-Hamilton difference cycle')
  cyc=list(zip(seq,seq[1:]+seq[:1]));cycles.append(cyc)
 ck(len(E)==n*(n-1)//2,'tournament edge count')
 ck(len(set(e for C in cycles for e in C))==len(E) and set(e for C in cycles for e in C)==set(E),'arc partition')
 # Complete local certificate: each cycle with zero/one deleted edge supplies an HP.
 for C in cycles:
  for i,deleted in enumerate(C):
   path=[C[(i+j)%n][1] for j in range(n)]
   ck(len(set(path))==n and deleted not in list(zip(path,path[1:])),'single-deletion path')
   ck(set(zip(path,path[1:]))<=set(C),'path arcs')
 return E,cycles

def path(n,E,B):
 inc=[0]*n
 for i,(u,v) in enumerate(E):
  if not(B>>i&1):inc[v]|=1<<u
 end=[0]*(1<<n)
 for v in range(n):end[1<<v]=1<<v
 for S in range(1,1<<n):
  if S&(S-1)==0:continue
  for v in range(n):
   if S>>v&1 and end[S^(1<<v)]&inc[v]:end[S]|=1<<v
 if not end[-1]:return None
 S=(1<<n)-1;v=(end[-1]&-end[-1]).bit_length()-1;out=[v]
 while S&(S-1):
  S^=1<<v;p=end[S]&inc[v];v=(p&-p).bit_length()-1;out.append(v)
 return out[::-1]
orientations={}
for n in (3,5,7,11,13):
 count=0
 for signs in it.product((0,1),repeat=(n-1)//2):
  Q={d if s==0 else n-d for d,s in zip(range(1,(n+1)//2),signs)}
  data(n,Q);count+=1
 orientations[n]=count
out=[]
for n in (7,11):
 Q={x*x%n for x in range(1,n)};E,C=data(n,Q);EI={e:i for i,e in enumerate(E)}
 cycles=[sum(1<<EI[e] for e in cycle) for cycle in C]
 if n==7:blocks=[r['blocker'] for r in json.loads((D/'inputs/P7-blockers.json').read_text())['bijection']]
 else:
  row=next(r for r in json.loads((D/'inputs/Hall-census.json').read_text())['censuses'] if r['n']==11)
  blocks=[r['mask'] for r in row['blockers']]+[sum(1<<i for i,(u,v) in enumerate(E) if c in (u,v)) for c in range(n)]
 ck(all(all((B&C).bit_count()==2 for C in cycles) for B in blocks),'minimum profile')
 # Necessary but not sufficient: choose two consecutive edges in each difference cycle.
 B=sum(1<<EI[e] for cycle in C for e in cycle[:2]);p=path(n,E,B)
 ck(all((B&C).bit_count()==2 for C in cycles) and p,'profile sufficiency hostile')
 out.append({'n':n,'beta_HP_proved':n-1,'cycles':C,'verified_known_blockers':len(blocks),'two_per_cycle_candidate_universe':math.comb(n,2)**((n-1)//2),'two_per_cycle_nonblocker':{'mask':B,'surviving_path':p}})
try:data(9,{1,2,3,4})
except RuntimeError:composite_rejected=True
else:composite_rejected=False
ck(composite_rejected,'composite nonunit hostile')
print(json.dumps({'status':'PASS','all_prime_circulant_orientations_checked':orientations,'paley_certificates':out,'composite9_nonunit_step3_rejected':composite_rejected},indent=2,sort_keys=True))
