#!/usr/bin/env python3
"""Independent exact arithmetic and graph checks; assertions not used."""
from pathlib import Path
from itertools import permutations, combinations
from math import isqrt
from fractions import Fraction
import hashlib,json,subprocess,sys
D=Path(__file__).resolve().parent
checks=[]
def require(ok,label):
 if not ok: raise RuntimeError(label)
 checks.append(label)
def lines(g): return list(g)+list(zip(*g))+[[g[i][i] for i in range(3)],[g[i][2-i] for i in range(3)]]
def edges(p): return {tuple(sorted((a,b))) for a,b in zip(p,p[1:])}
g=[[15,35,12],[57,1,4],[49,26,46]]; entries={v for row in g for v in row}; squares={i*i%59 for i in range(59)}
require(len(entries)==9 and entries <= squares-{0},'Nine distinct nonzero square residues')
require(all(sum(row)%59==3 for row in lines(g)),'All eight modular magic totals equal 3')
require(all((57+v)%59 not in squares for v in entries-{57}),'57 isolated for both square-sum conventions')
A=[48,0,11];B=[46,1,15];L=[[(a+b)%59 for b in B] for a in A]
require(L==[[35,49,4],[46,1,15],[57,12,26]],'Cross-label matrix exact')
require({v for row in L for v in row}==entries,'Edge-label multiset equals original grid entries')
require(all((A[i+1]-A[i])%59==11 for i in range(2)) and all((B[i+1]-B[i])%59==14 for i in range(2)),'Both AP guards')
require(len(set(A+B))==6,'Six distinct carrier potentials')
E={(i,j) for i in range(3) for j in range(3,6)}
P=(0,3,1,4,2,5);Q=(4,0,5,1,3,2)
require(set(P)==set(Q)==set(range(6)) and edges(P)<=E and edges(Q)<=E,'Both displayed orders are Hamilton paths')
require(edges(P)|edges(Q)==E and edges(P)&edges(Q)=={(1,3)},'Nine-edge union and unique shared edge A1 B0, label 46')
HP=[p for p in permutations(range(6)) if p<p[::-1] and edges(p)<=E]
require(len(HP)==36,'36 unoriented K33 Hamilton paths')
require(sum(edges(p)|edges(q)==E for p,q in combinations(HP,2))==36,'36 unordered covering pairs')
rows=[]
for i in (1,2):
 for j in (1,2):
  r=[0]*9
  for k,c in [(3*i+j,1),(3*i,-1),(j,-1),(0,1)]:r[k]+=c
  rows.append(r)
rows.extend([[1,0,0,-2,0,0,1,0,0],[1,-2,1,0,0,0,0,0,0]])
m=[[Fraction(v) for v in r] for r in rows];rank=0
for c in range(9):
 p=next((i for i in range(rank,len(m)) if m[i][c]),None)
 if p is None:continue
 m[rank],m[p]=m[p],m[rank];x=m[rank][c];m[rank]=[v/x for v in m[rank]]
 for i in range(len(m)):
  if i!=rank:
   x=m[i][c];m[i]=[a-x*b for a,b in zip(m[i],m[rank])]
 rank+=1
require(rank==6,'Six independent linear label constraints over Q')
require(all(sum(a*b for a,b in zip(r,sum(L,[])))%59==0 for r in rows),'Displayed modular labels satisfy all six constraints')
roots=[[2655,2976,320],[2624,2520,1665],[1440,895,3624]]; G=[[v*v for v in r] for r in roots]; totals=[sum(r) for r in lines(G)]
require(len({v for r in G for v in r})==9 and min(sum(G,[]))>0,'Nine distinct positive integer square entries')
require(totals==[16008001]*6+[26532801,8526400],'All eight integer line totals independently recomputed')
require([isqrt(t) for t in totals]==[4001]*6+[5151,2920] and all(isqrt(t)**2==t for t in totals),'Every integer line total is square')
require(len(set(totals))==3,'Integer variant is not magic')
for pin in json.loads((D/'SOURCE_PINS.json').read_text()):
 require(hashlib.sha256((D/pin['local']).read_bytes()).hexdigest()==pin['sha256'],'Source SHA256 '+pin['local'])
names=['edge-label-repair.svg','square-totals-not-magic.svg'];before={n:hashlib.sha256((D/n).read_bytes()).hexdigest() for n in names}
subprocess.run([sys.executable,str(D/'generate.py')],check=True,stdout=subprocess.DEVNULL)
require(before=={n:hashlib.sha256((D/n).read_bytes()).hexdigest() for n in names},'Deterministic SVG replay is byte-identical')
print(json.dumps({'status':'PASS','checks':checks,'check_count':len(checks),'svg_sha256':before},indent=2))
