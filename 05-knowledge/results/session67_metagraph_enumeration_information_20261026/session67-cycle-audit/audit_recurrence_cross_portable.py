#!/usr/bin/env python3
"""Portable derivative of the frozen recurrence cross-audit.

Only peer-lane path selection is changed; numerical checks are identical.
Use --recursion-dir to select the frozen producer lane at any location.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb,factorial
from pathlib import Path
import argparse,hashlib,importlib.util,json
from verify import adjacency,hpairs

B=Path(__file__).resolve().parent
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--recursion-dir',type=Path,default=B.parent/'recursion',
                help='Peer lane containing REPORT.md and verify_recursion.py')
args=ap.parse_args();P=args.recursion_dir.resolve();checks=Counter()
def need(ok,what):
 checks[what]+=1
 if not ok:raise RuntimeError(what)
def triangle(w):
 own=[i for i,k in enumerate(w) for _ in range(k)];rows=[0]*sum(w)
 for a,b in combinations(range(len(own)),2):
  u,v=(a,b) if own[a]==own[b] or (own[b]-own[a])%3==1 else (b,a)
  rows[u]|=1<<v
 return tuple(rows)
def ext(rows,mask):
 n=len(rows);return tuple(r|((1-(mask>>v&1))<<n) for v,r in enumerate(rows))+(mask,)
def chains(rows):
 par=list(range(len(rows)))
 def root(i):
  while par[i]!=i:i=par[i]
  return i
 for i,j in hpairs(rows):par[root(i)]=root(j)
 return sorted(Counter(root(i) for i in par).values())
def auts(rows):
 # Exact arc-consistent bijection backtracking, not peer's degree block min.
 n=len(rows);deg=[r.bit_count() for r in rows]
 candidates=[[j for j in range(n) if deg[j]==deg[i]] for i in range(n)]
 order=sorted(range(n),key=lambda i:(len(candidates[i]),i));mapping=[-1]*n;out=[]
 def rec(t,used):
  if t==n:out.append(tuple(mapping));return
  v=order[t]
  for w in candidates[v]:
   if used>>w&1:continue
   if all(((rows[v]>>u)&1)==((rows[w]>>mapping[u])&1) for u in order[:t]):
    mapping[v]=w;rec(t+1,used|1<<w);mapping[v]=-1
 rec(0,0);return out
def brute_code(rows,order):
 return sum(((rows[order[a]]>>order[b])&1)<<k for k,(a,b) in enumerate(combinations(range(len(rows)),2)))
def rooted_code(rows,v):
 return min(brute_code(rows,(v,)+p) for p in permutations([u for u in range(len(rows)) if u!=v]))

families=[]
for k in range(2,13):
 rows_out=[]
 for name,weights in [('A',(3*k-3,1,1)),('B',(k-1,k,k))]:
  rows=triangle(weights);n=len(rows);group=auts(rows)
  need(len(group)==1,'both parents rigid')
  need(chains(rows)==sorted(weights),'parents exact canonical weights')
  masks={r|bit<<v for v,r in enumerate(rows) for bit in (0,1)}
  need(len(masks)==n+3,'literal clone mask count')
  aut_hist=Counter();total=0;uniform=0
  for mask in masks:
   child=ext(rows,mask);sizes=chains(child);g=len(auts(child));aut_hist[g]+=1
   need(len(sizes)==3,'each literal clone mask keeps core3')
   total+=g
   if sizes==[k]*3:uniform+=g
  mass=Q(total,n+1);umass=Q(uniform,n+1)
  expected=Q(3*k+2+(2*k if name=='B' else 0),3*k)
  need(mass==expected,'infinite formula exact at sampled k')
  need(umass==(1 if name=='B' else 0),'uniform child exact 0 versus1 transition')
  row=dict(k=k,name=name,n=n,clone_aut_histogram=dict(aut_hist),core3_mass=str(mass),uniform_child_mass=str(umass))
  if k<=4:
   raw=Counter();weighted=Counter();core3_masks=set();child_aut=Counter()
   for mask in range(1<<n):
    child=ext(rows,mask);sizes=chains(child);g=len(auts(child));raw[len(sizes)]+=1;weighted[len(sizes)]+=g;child_aut[g]+=1
    if len(sizes)==3:core3_masks.add(mask)
   need(core3_masks==masks,'all masks core3 iff literal clones')
   need(Q(weighted[3],n+1)==mass,'full attachment and clone-only mass agree')
   row.update(raw_histogram=dict(raw),weighted_histogram={b:str(Q(v,n+1)) for b,v in weighted.items()},all_child_aut_histogram=dict(child_aut))
  rows_out.append(row)
 need(Q(rows_out[1]['core3_mass'])-Q(rows_out[0]['core3_mass'])==Q(2,3),'infinite hostile gap two thirds')
 families.extend(rows_out)
# Inspect producer's rooted-key implementation against complete brute rooted codes.
s=importlib.util.spec_from_file_location('peer_recursion',P/'verify_recursion.py');peer=importlib.util.module_from_spec(s);s.loader.exec_module(peer)
rootrows=[]
for n in range(1,7):
 seen=set();nclass=nroots=0
 for code in range(1<<comb(n,2)):
  if code in seen:continue
  rows=tuple(adjacency(n,code));orb={brute_code(rows,p) for p in permutations(range(n))};seen.update(orb);nclass+=1
  g=auts(rows);need(len(g)*len(orb)==factorial(n),'backtracking automorphisms full orbit mass')
  exact=[rooted_code(rows,v) for v in range(n)];pruned=[peer.root_key(rows,v) for v in range(n)]
  selected=[v for v,c in enumerate(pruned) if c==min(pruned)]
  need(set(selected)=={p[selected[0]] for p in g},'selected roots exactly one actual automorphism orbit')
  for u in range(n):
   for v in range(n):need((exact[u]==exact[v])==(pruned[u]==pruned[v]),'rooted code equality iff exact rooted isomorphism')
  nroots+=len(set(exact))
 rootrows.append(dict(n=n,classes=nclass,rooted_classes=nroots))
need(rootrows[1]['rooted_classes']==2 and rootrows[2]['rooted_classes']==4,'minimal card and divide-order boundaries')
out=dict(status='PASS',
 frozen_parent_script_sha256='84ee2a8c0eabc9655fae98230a7967240415a23750a99694ce8327138c2ea613',
 checks=dict(sorted(checks.items())),total_checks=sum(checks.values()),sampled_k=families,exact_rooted_census=rootrows,
 peer_report_sha256=hashlib.sha256((P/'REPORT.md').read_bytes()).hexdigest(),
 peer_script_sha256=hashlib.sha256((P/'verify_recursion.py').read_bytes()).hexdigest())
print(json.dumps(out,indent=2,sort_keys=True))
