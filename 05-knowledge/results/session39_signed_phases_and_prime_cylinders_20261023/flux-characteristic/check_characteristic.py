#!/usr/bin/env python3
import csv,gzip,json,hashlib,random,sys
from pathlib import Path
from fractions import Fraction as F
base=Path(sys.argv[1]) if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'session38-wiring'
def require(c,m):
 if not c: raise RuntimeError(m)
def rows(name):
 with gzip.open(base/name,'rt') as f: return [{k:int(v) for k,v in r.items()} for r in csv.DictReader(f,delimiter='\t')]
nodes=rows('nodes.tsv.gz'); edges=rows('edges.tsv.gz'); N=len(nodes)
deg=[0]*N; diag=[0]*N
for e in edges:
 u,v,w=e['source'],e['target'],e['blue']+e['black']
 deg[u]+=w; deg[v]+=w
 if u!=v: diag[u]+=w; diag[v]+=w
for r in nodes:
 i=r['node_id']; h=r['H']//r['Aut']; require(r['H']%r['Aut']==0 and h%2==1,'odd fibre')
 require(deg[i]==h*(2-r['SC']),'literal degree identity')
 require(diag[i]%2==r['SC'],'characteristic indicator')
root=0; rng=random.Random(39)
for _ in range(128):
 z=[rng.randrange(-3,4) for _ in range(N)]; z[root]=0
 energy=sum((e['blue']+e['black'])*(z[e['source']]-z[e['target']])**2 for e in edges)
 require(energy%2==sum(r['SC']*z[r['node_id']] for r in nodes)%2,'native local parity')
# Independent odd-lattice characteristic refinement descent on a 2x2 example.
K=((3,-1),(-1,2)); inv=((F(2,5),F(1,5)),(F(1,5),F(3,5))); w=(1,0)
def dotinv(x,y): return sum(x[i]*inv[i][j]*y[j] for i in range(2) for j in range(2))
def q(x):return (dotinv(x,x)-dotinv(x,w))/2
count=0
for a in range(-3,4):
 for b in range(-3,4):
  x=(a,b)
  for c in range(-2,3):
   for d in range(-2,3):
    n=(c,d); y=tuple(x[i]+sum(K[i][j]*n[j] for j in range(2)) for i in range(2))
    require((q(y)-q(x)).denominator==1,'refinement descent');count+=1
print(json.dumps({'status':'PASS','native_order':8,'nodes':N,'edges':len(edges),'self_converse_nodes':sum(r['SC'] for r in nodes),'root':root,'retained_characteristic_ones':sum(r['SC'] for r in nodes if r['node_id']!=root),'degree_and_diagonal_checks':N,'native_integer_potential_checks':128,'small_refinement_descent_checks':count,'input_sha256':{n:hashlib.sha256((base/n).read_bytes()).hexdigest() for n in ['nodes.tsv.gz','edges.tsv.gz']}},indent=2))
