#!/usr/bin/env python3
"""Complete native line matrices and full-slot partition closures through n=7."""
from array import array
from collections import Counter,defaultdict
from fractions import Fraction as F
from itertools import combinations,permutations
from math import comb,factorial
from pathlib import Path
import argparse,hashlib,json,subprocess

HERE=Path(__file__).resolve().parent
CHECKS=0
def need(ok,msg):
 global CHECKS
 CHECKS+=1
 if not ok:raise RuntimeError(str(msg))
def pin(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def orbit_table(n):
 scratch=HERE/'_scratch';scratch.mkdir(exist_ok=True);exe=scratch/'orbits'
 if not exe.exists() or exe.stat().st_mtime<(HERE/'orbits.cpp').stat().st_mtime:
  subprocess.run(['c++','-O2','-std=c++17',str(HERE/'orbits.cpp'),'-o',str(exe)],check=True)
 path=scratch/f'owner{n}.bin'
 subprocess.run([str(exe),str(n),str(path)],check=True,stdout=subprocess.DEVNULL)
 out=array('I');out.frombytes(path.read_bytes());need(len(out)==1<<(n*(n-1)//2),'orbit universe')
 return out

def partition(keys):
 groups=defaultdict(list)
 for i,k in enumerate(keys):groups[k].append(i)
 return sorted(groups.values(),key=lambda g:g[0])
def refine(cells,channels,d,normalized):
 trace=[len(cells)];levels=[cells]
 while True:
  labels={v:i for i,g in enumerate(cells) for v in g};sig=[]
  for i in range(len(d)):
   rows=tuple(tuple(F(sum(M[i].get(j,0) for j in g),d[i] if normalized else 1) for g in cells) for M in channels)
   sig.append((labels[i],rows))
  new=partition(sig);trace.append(len(new))
  if new==cells:break
  cells=new;levels.append(cells)
 return {'counts':trace,'strict_rounds':len(levels)-1,'cells':cells,'levels':levels}

def profiles(n,code,E,pos,allp):
 def arc(i,j):return bool((code>>pos[min(i,j),max(i,j)])&1)==(i<j)
 literal=[0]*(n+1)
 for p in allp:
  bad=sum(not arc(a,b) for a,b in zip(p,p[1:]));good=n-1-bad
  for extra in range(good+1):literal[bad+extra+1]+=comb(good,extra)
 hp=[0]*(1<<n);hp[0]=1;ends=[[0]*n for _ in hp]
 for mask in range(1,1<<n):
  for j in range(n):
   if mask>>j&1:
    old=mask^(1<<j);ends[mask][j]=1 if not old else sum(ends[old][i] for i in range(n) if old>>i&1 and arc(i,j))
  hp[mask]=sum(ends[mask])
 covers=[[0]*(1<<n) for _ in range(n+1)];covers[0][0]=1
 for j in range(1,n+1):
  for mask in range(1,1<<n):
   sub=mask
   while sub:
    covers[j][mask]+=hp[sub]*covers[j-1][mask^sub];sub=(sub-1)&mask
 need(literal==[row[-1] for row in covers],('two full cover paths',n,code))
 Q=[F(0)]*(n+1);falling=[F(1)]
 for j in range(1,n+1):
  new=[F(0)]*(j+1)
  for k,a in enumerate(falling):new[k]-=a*F(j-1,j);new[k+1]+=a/F(j)
  falling=new
  for k,a in enumerate(falling):Q[k]+=literal[j]*a
 need(all(a.denominator==1 for a in Q),('integral Q',n,code))
 need(all(Q[k]==0 for k in range(n+1) if k%2!=n%2),('slot parity',n,code))
 return literal[1:],list(map(int,Q))

def census(n):
 E=list(combinations(range(n),2));pos={e:k for k,e in enumerate(E)};top=(1<<len(E))-1
 owner=orbit_table(n);orbits=Counter(owner);ordinary=sorted(orbits);oi={c:i for i,c in enumerate(ordinary)}
 dual={c:int(owner[top^c]) for c in ordinary};merged=sorted({min(c,dual[c]) for c in ordinary});mi={c:i for i,c in enumerate(merged)}
 aut={c:factorial(n)//orbits[c] for c in ordinary};O=[Counter() for _ in ordinary];OB=[Counter() for _ in ordinary]
 free=[e for e in E if e[1]>e[0]+1];spine=sum(1<<pos[i,i+1] for i in range(n-1));freebits=top^spine
 def arc(t,i,j):return bool((t>>pos[min(i,j),max(i,j)])&1)==(i<j)
 def blue(t,p):return all(arc(t,p[i],p[j])==arc(t,p[n-1-j],p[n-1-i]) for i,j in E)
 for w in range(1<<len(free)):
  t=spine|sum(((w>>k)&1)<<pos[e] for k,e in enumerate(free));u=int(owner[t]);v=int(owner[t^freebits])
  O[oi[u]][oi[v]]+=1
  if blue(t,tuple(range(n))):OB[oi[u]][oi[v]]+=1
 allp=list(permutations(range(n)));H={};path_entries=0
 for c in ordinary:
  counts=Counter();bcounts=Counter();h=0
  for p in allp:
   if not all(arc(c,a,b) for a,b in zip(p,p[1:])):continue
   h+=1;pathmask=sum(1<<pos[min(a,b),max(a,b)] for a,b in zip(p,p[1:]));u=int(owner[c^top^pathmask]);counts[oi[u]]+=1
   if blue(c,p):bcounts[oi[u]]+=1
  H[c]=h
  need(all(v%aut[c]==0 for v in counts.values()),('ordinary orbit divisors',n,c))
  need(O[oi[c]]==Counter({j:v//aut[c] for j,v in counts.items()}),('native two incidence paths',n,c))
  need(OB[oi[c]]==Counter({j:v//aut[c] for j,v in bcounts.items()}),('blue two incidence paths',n,c))
  path_entries+=h
 # Unrestricted minimum is independently checked on every ordinary representative.
 for c in ordinary:
  need(min(sum(arc(c,p[i],p[j])<<k for k,(i,j) in enumerate(E)) for p in allp)==c,('unrestricted code',n,c))
 A=[Counter() for _ in merged];B=[Counter() for _ in merged]
 for i,c in enumerate(ordinary):
  u=mi[min(c,dual[c])]
  for j,v in O[i].items():A[u][mi[min(ordinary[j],dual[ordinary[j]])]]+=v
  for j,v in OB[i].items():B[u][mi[min(ordinary[j],dual[ordinary[j]])]]+=v
 d=list(map(lambda row:sum(row.values()),A));K=[A[i]-B[i] for i in range(len(A))]
 need(len(ordinary)=={3:2,4:4,5:12,6:56,7:456}[n],('ordinary class control',n))
 need(len(merged)=={3:2,4:3,5:10,6:34,7:272}[n],('merged class control',n))
 need(sum(d)==1<<len(free),'dart total')
 need(sum(sum(row.values()) for row in B)==1<<((n-1)**2//4),'blue total')
 nodes=[];qcache={}
 for i,c in enumerate(merged):
  covers,Q=profiles(n,c,E,pos,allp);qcache[c]=(covers,Q)
  score=sorted(sum(arc(c,a,b) for b in range(n) if a!=b) for a in range(n));fold=min(score,sorted(n-1-x for x in score))
  need(covers[0]==H[c],('H vs covers',n,c));need(d[i]*aut[c]==(1+(dual[c]!=c))*H[c],('source mass',n,c))
  nodes.append({'code':c,'ordinary_codes':[c] if c==dual[c] else sorted([c,dual[c]]),'H':H[c],'aut':aut[c],
   'self_converse':c==dual[c],'degree':d[i],'blue_degree':sum(B[i].values()),'score':score,'folded_score':fold,
   'ordered_covers':covers,'Q':Q})
 for M in (A,B,K):
  for i in range(len(merged)):
   need(M[i].get(i,0)%2==0,'merged loops doubled')
   for j,v in M[i].items():need(M[j].get(i,0)==v,'merged symmetric')
 for i in range(len(merged)):
  for v in K[i].values():need(v%2==0,'black twins')
 for i,c in enumerate(ordinary):
  need(sum(O[i].values())*aut[c]==H[c],'ordinary source mass')
  for j,v in O[i].items():need(O[j].get(i,0)==v,'ordinary native symmetry')
 # Explicit observer keys; H is already Q(1), retained only as a diagnostic.
 keys={
  'Q':[tuple(x['Q']) for x in nodes],
  'Q_augmented':[(tuple(x['Q']),tuple(x['folded_score']),x['aut'],x['self_converse'],x['degree']) for x in nodes],
  'Q_degree':[(tuple(x['Q']),x['degree']) for x in nodes],
  'Q_blue_degree':[(tuple(x['Q']),x['blue_degree']) for x in nodes],
  'Q_blue_parity':[(tuple(x['Q']),x['blue_degree']%2) for x in nodes],
  'Q_TT_incidence':[(tuple(x['Q']),bool(A[i].get(mi[0],0))) for i,x in enumerate(nodes)],
  'Q_SC_TT_incidence':[(tuple(x['Q']),x['self_converse'],bool(A[i].get(mi[0],0))) for i,x in enumerate(nodes)]}
 closure={}
 for label,key in keys.items():
  initial=partition(key);closure[label]={}
  for normalized in (False,True):
   for colored in (False,True):
    name=('P' if normalized else 'A')+('_colored' if colored else '')
    result=refine(initial,[B,K] if colored else [A],d,normalized)
    result['cells']=[[merged[i] for i in cell] for cell in result['cells']]
    result['levels']=[[[merged[i] for i in cell] for cell in level] for level in result['levels']]
    closure[label][name]=result
 need(len(partition(keys['Q']))=={3:2,4:3,5:8,6:24,7:156}[n],('inherited full Q class count',n))
 def sparse(M,codes):return [[codes[i],codes[j],v] for i,row in enumerate(M) for j,v in sorted(row.items())]
 ordinary_nodes=[{'code':c,'converse':dual[c],'merged_code':min(c,dual[c]),'H':H[c],'aut':aut[c],
  'self_converse':c==dual[c],'degree':sum(O[oi[c]].values()),'blue_degree':sum(OB[oi[c]].values()),
  'Q':qcache[min(c,dual[c])][1],'ordered_covers':qcache[min(c,dual[c])][0]} for c in ordinary]
 return {'n':n,'labelled_universe':len(owner),'fixed_spines':1<<len(free),'physical_paths_checked':path_entries,
  'ordinary_nodes':ordinary_nodes,'merged_nodes':nodes,
  'ordinary_A':sparse(O,ordinary),'ordinary_blue':sparse(OB,ordinary),
  'merged_A':sparse(A,merged),'merged_blue':sparse(B,merged),'closure':closure}

def main():
 parser=argparse.ArgumentParser();parser.add_argument('--max-n',type=int,default=7);parser.add_argument('--output',default='data.json');args=parser.parse_args()
 rows=[]
 for n in range(3,args.max_n+1):
  rows.append(census(n));print('completed',n,'nodes',len(rows[-1]['merged_nodes']),flush=True)
  data={'status':'FINITE-EXACT; proof candidate RESERVED until audit','base_commit':'a8b6d8dbf0c6b601813785efb0b5904b46c23d43',
   'gauge':'bit1 means increasing arc; lexicographic pairs; unrestricted minimum over all vertex permutations; native nonpath complement before converse merge',
   'matrix_schema':'sparse triples [source canonical code,target canonical code,integer dart count]; loops doubled; black=A-blue',
   'Q_schema':'Q monomial coefficients in ascending exponent; ordered_covers[j-1]=F_j with j nonempty labelled path slots',
   'checks':CHECKS,'source_sha256':pin(Path(__file__)),'orbits_source_sha256':pin(HERE/'orbits.cpp'),'rows':rows}
  (HERE/args.output).write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
 print(json.dumps({'status':'PASS','checks':CHECKS,'max_n':args.max_n,'data_sha256':pin(HERE/args.output),
                   'source_sha256':pin(Path(__file__)),'orbits_source_sha256':pin(HERE/'orbits.cpp')},sort_keys=True))
if __name__=='__main__':main()
