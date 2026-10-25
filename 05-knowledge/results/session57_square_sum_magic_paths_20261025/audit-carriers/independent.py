#!/usr/bin/env python3
"""Independent bounded audit; no imports from either research verifier."""
from pathlib import Path
from math import isqrt,gcd
from itertools import combinations
import json,hashlib
ROOT=Path(__file__).resolve().parent
checks=0
def need(ok,msg):
 global checks
 checks+=1
 if not ok:raise RuntimeError(msg)
def rejects(f):
 try:f()
 except (RuntimeError,ValueError):return True
 raise RuntimeError('negative control was accepted')
def sq(n):return n>=0 and isqrt(n)**2==n
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def validate_manifest(root):
 for line in (root/'SHA256SUMS').read_text().splitlines():
  h,n=line.split(maxsplit=1);need(digest(root/n)==h,'manifest '+n)
for packet in ('carriers','inheritance'):validate_manifest(ROOT/'input'/packet)
I=ROOT/'input/inheritance'; C=ROOT/'input/carriers'
files=set(json.loads((I/'FILES.json').read_text()))|{'FILES.json','SHA256SUMS'}
need({str(p.relative_to(I)) for p in I.rglob('*') if p.is_file()}==files,'publication closure')
for p in json.loads((I/'SOURCE_PINS.json').read_text()):need(digest(I/p['snapshot'])==p['sha256'],'pin '+p['id'])
for i,p in enumerate(json.loads((C/'SOURCE_PINS.json').read_text())['internal'],1):need(digest(ROOT/f'input/carrier-dependencies/S{i:02}.md')==p['sha256'],'carrier source snapshot pin')
# Generate G_n by square shells rather than evaluating all ordered pairs.
hcounts=[]
for n in range(15,26):
 adj=[set() for _ in range(n+1)]
 for k in range(2,isqrt(2*n)+1):
  for x in range(1,n+1):
   y=k*k-x
   if x<y<=n:adj[x].add(y);adj[y].add(x)
 found=set()
 def walk(path,unused):
  if not unused:
   p=tuple(path);found.add(min(p,p[::-1]));return
  for y in adj[path[-1]]&unused:walk(path+[y],unused-{y})
 for x in range(1,n+1):walk([x],set(range(1,n+1))-{x})
 hcounts.append(len(found))
 need(sum(map(len,adj))//2==json.loads((C/'certificate.json').read_text())['square_sum_graphs'][n-15]['edges'],'graph edges')
need(hcounts==[1,1,1,0,0,0,0,0,3,0,10],'Hamilton census')
# Independent Euclid generation of all opposite-square gaps up to e=1105.
# Nonprimitive triples are represented by their unique positive gcd scale.
limit=1105; euclid={e:set() for e in range(1,limit+1)}
for m in range(1,isqrt(limit)+1):
 for n in range(1,m):
  c=m*m+n*n
  if c>limit or gcd(m,n)!=1 or (m-n)%2==0:continue
  A,B=sorted((m*m-n*n,2*m*n),reverse=True)
  for g in range(1,limit//c+1):euclid[g*c].add(2*g*g*A*B)
petals=0; allgaps=0; examples={}
for e in range(1,limit+1):
 direct={e*e-a*a for a in range(1,e) if sq(2*e*e-a*a)}
 need(direct==euclid[e],'Euclid completeness '+str(e));allgaps+=len(direct)
 images=set()
 for d in direct:
  a=isqrt(e*e-d);b=isqrt(e*e+d)
  need((a+b)%2==0 and (b-a)%2==0,'inverse parity')
  P,Q=(b+a)//2,(b-a)//2;g=gcd(P,Q);A,B,C0=P//g,Q//g,e//g
  need(A>B>0 and gcd(A,B)==1 and A*A+B*B==C0*C0,'inverse primitive')
  need(e==g*C0 and d==2*g*g*A*B,'inverse exact')
  t=(2*e*e-4*d,2*e*e,2*e*e+4*d)
  need((t[0]>0)==(C0*C0>4*A*B),'positivity')
  need(all(gcd(abs(t[i]),abs(t[j]))==2*g*g for i,j in combinations(range(3),2)),'vertex gcd')
  roots=(2*g*(A-B),2*e,2*g*(A+B))
  need(all(gcd(roots[i],roots[j])==2*g for i,j in combinations(range(3),2)),'root gcd')
  if t[0]>0:images.add(t)
 # Alternate inverse enumerates roots of the smaller pair sum, not gaps.
 triangles=set()
 for p in range(1,2*e):
  r2=8*e*e-p*p
  if sq(r2) and p*p>2*e*e:triangles.add((p*p-2*e*e,2*e*e,6*e*e-p*p))
 need(images==triangles,'root-pair enumeration '+str(e))
 if e<=500:petals+=len(images)
 if e in (5,41,65,1105):
  diamonds=[(x,y) for x in direct for y in direct if x>y and x!=2*y and x-y in direct and x+y in direct]
  examples[str(e)]={'total_gaps':len(direct),'short_gaps':len(images),'diamonds':diamonds}
need(petals==68 and examples['1105']=={'total_gaps':13,'short_gaps':5,'diamonds':[]},'petal counts')
# Derive grid lines from geometry, not research constants.
lines=[tuple(3*r+c for c in range(3)) for r in range(3)]+[tuple(3*r+c for r in range(3)) for c in range(3)]+[(0,4,8),(2,4,6)]
need(sum(bool(set(a)&set(b)) for a,b in combinations(lines,2))==22,'line geometry')
# Universal parametrization is checked in all Z/nZ, 2<=n<=19.
for mod in range(2,20):
 for c in range(mod):
  for a in range(mod):
   for b in range(mod):
    v=[a,3*c-a-b,b,c-a+b,c,c+a-b,2*c-b,a+b-c,2*c-a]
    need(all(sum(v[i] for i in l)%mod==3*c%mod for l in lines),'ring parametrization')
# Original arithmetic witnesses and shared-label path recovery.
for mod,roots in [(29,[9,11,1,6,0,14,12,16,8]),(59,[20,12,7,2,1,23,22,25,29])]:
 v=[r*r%mod for r in roots];q={r*r%mod for r in range(mod)};T=3*v[4]%mod
 need(len(set(v))==9 and set(v)<=q,'field injectivity')
 paths=[]
 for line in lines:
  prefix=[0]
  for i in line:prefix.append((prefix[-1]+v[i])%mod)
  need(prefix[-1]==T,'field sum');paths.append(prefix)
 def decode(paths):
  cell={}
  for line,path in zip(lines,paths):
   need(len(path)==4 and path[0]==0 and path[-1]==T,'path endpoints')
   for j,i in enumerate(line):
    value=(path[j+1]-path[j])%mod
    need(value in q,'square increment')
    need(i not in cell or cell[i]==value,'shared cell agreement');cell[i]=value
  need(len(set(cell.values()))==9,'injective cell assignment')
  return [cell[i] for i in range(9)]
 need(decode(paths)==v,'inverse path carrier')
 bad=[p[:] for p in paths];bad[0][1]=(bad[0][1]+1)%mod
 rejects(lambda:decode(bad))
q={x*x%29 for x in range(29)}
need(all(((-a-b)%29 in q) for a,b in combinations((0,1,4),2)) and sum((0,1,4))%29!=0,'shadow false triangle')
# Negative controls that remain active under -O.
rejects(lambda:need(13==5,'capacity loses positivity'))
rejects(lambda:need(sq(2*1105**2),'opposite-shell misidentification'))
rejects(lambda:need(sq(3*1105**2),'square-total misidentification'))
rejects(lambda:need(hcounts==[1]*11,'fake Hamilton census'))
print(json.dumps({'status':'PASS','checks':checks,'Hamiltonian_counts':hcounts,'center_root_limit':limit,'total_gaps_checked':allgaps,'positive_petals_e_le_500':petals,'examples':examples,'negative_controls':6,'source_scope':'14 inheritance snapshots hashed; 3 original carrier pins hashed; no promotion of source conjectures'},indent=2,sort_keys=True))
