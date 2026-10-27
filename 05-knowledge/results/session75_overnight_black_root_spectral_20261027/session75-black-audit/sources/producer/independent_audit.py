#!/usr/bin/env python3
"""Independent arc-set + score-vector enumeration audit. No producer imports."""
import json,itertools,functools,math,hashlib,sys
from pathlib import Path

def need(p,s):
 if not p:raise RuntimeError(s)
def arcs(rows):return {(i,j) for i,r in enumerate(rows) for j in range(len(rows)) if r>>j&1}
def pathblue(E,p):
 r=dict(zip(p,p[::-1]));return all((r[v],r[u]) in E for u,v in E)
def coretarget(E,n):
 for v in range(n):
  d=sum((v,u) in E for u in range(n))
  if d not in (0,n-1):continue
  S=set(range(n))-{v}
  good=True
  for root in S:
   seen={root}
   while True:
    new=seen|{y for x,y in E if x in seen and y in S}
    if new==seen:break
    seen=new
   if seen!=S:good=False;break
  if good:return True
 return False
@functools.lru_cache(None)
def scorecount(d):
 n=len(d)
 if n==0:return 1
 if min(d)<0 or max(d)>=n or sum(d)!=n*(n-1)//2:return 0
 total=0
 for outs in itertools.combinations(range(1,n),d[0]):
  out=set(outs)
  residual=tuple(sorted(d[i]-(i not in out) for i in range(1,n)))
  total+=scorecount(residual)
 return total

def check_record(rec):
 n=len(rec['initial_rows']);E=arcs(rec['initial_rows'])
 requirecuts=0
 for move in rec['moves']:
  p=move['K_path'];need(sorted(p)==list(range(n)),'path labels')
  need(E==arcs(move['K_source_rows']),'source carrier')
  need(all((u,v) in E for u,v in zip(p,p[1:])),'K word not directed')
  need(not pathblue(E,p),'K word blue')
  pe={frozenset((u,v)) for u,v in zip(p,p[1:])}
  E={(u,v) if frozenset((u,v)) in pe else (v,u) for u,v in E}
  need(E==arcs(move['K_target_rows']),'actual K edge labels')
  need(all((u,v) in E for u,v in zip(p,p[1:])),'K path directed after move')
  need(not pathblue(E,p),'K blackness not stable')
 need(coretarget(E,n),'core target lost')
 d=[sum((u,v) in arcs(rec['initial_rows']) for v in range(n)) for u in range(n)]
 need(len(rec['moves'])==min(min(x,n-1-x) for x in d),'score lower bound')
 return len(rec['moves'])
def main():
 base=Path(__file__).resolve().parent;raw=(base/'normal.json').read_bytes();r=json.loads(raw)
 checked=sum(check_record(x) for x in r['seeds']+r['examples']+[r['nonregular_target_separation']])
 orbits={};counts={}
 for n in (4,6,7):
  seeds=[x['initial_rows'] for x in r['seeds'] if len(x['initial_rows'])==n];U=set()
  for a in seeds:
   E=arcs(a)
   for perm in itertools.permutations(range(n)):
    F={(perm[u],perm[v]) for u,v in E}
    U.add(tuple(sum(1<<v for x,v in F if x==u) for u in range(n)))
  d=tuple([n//2-1]*(n//2)+[n//2]*(n//2)) if n%2==0 else ((n-1)//2,)*n
  exact=scorecount(d)*(math.comb(n,n//2) if n%2==0 else 1)
  need(len(U)==exact==int(r['complete_balanced_counts'][str(n)]),'orbit coverage vs independent score DP')
  need(all(max(x.bit_count() for x in a)-min(x.bit_count() for x in a)<=1 for a in U),'orbit balance')
  orbits[n]=len(U);counts[n]=exact
 hostile=[]
 for x in r['pure_blue_controls']:
  E=arcs(x['rows']);n=len(x['rows']);hp=[]
  for p in itertools.permutations(range(n)):
   if all((u,v) in E for u,v in zip(p,p[1:])):hp.append(p)
  need(len(hp)==x['hamilton_paths'] and all(pathblue(E,p) for p in hp),'pure blue permutation control')
  hostile.append([x['name'],len(hp)])
 print(json.dumps({'status':'PASS','producer_receipt_sha256':hashlib.sha256(raw).hexdigest(),'literal_K_moves_checked':checked,'complete_balanced_orbit_coverage':orbits,'independent_score_vector_counts':counts,'pure_blue_all_permutations':hostile,'independence':'no producer imports; arc sets and permutation anti-map; orbit coverage vs memoized score-count recurrence'},sort_keys=True,indent=2))
if __name__=='__main__':main()
