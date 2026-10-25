#!/usr/bin/env python3
"""Independent reconstruction; imports no producer Python."""
import pathlib,json,itertools,collections,hashlib
D=pathlib.Path(__file__).resolve().parent
P=D/'input'
def ck(x,m):
 if not x:raise RuntimeError(m)
ck(hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest()=='9b52b7d7495d418de0f561b8ae7e95c08736f316b2d0c7926a01a54d08248439','pinned manifest')
M=json.loads((P/'MANIFEST.json').read_text())['files']
ck({str(f.relative_to(P)) for f in P.rglob('*') if f.is_file()}==set(M)|{'MANIFEST.json'},'closed input')
for f,h in M.items():ck(hashlib.sha256((P/f).read_bytes()).hexdigest()==h,'digest '+f)
V=set(range(7));E=[(u,v) for u in range(7) for v in range(7) if (v-u)%7 in {1,2,4}];idx={e:i for i,e in enumerate(E)}
def mask(es):return sum(1<<idx[e] for e in es)
paths=[mask(zip(p,p[1:])) for p in itertools.permutations(range(7)) if all((a,b) in idx for a,b in zip(p,p[1:]))]
ck(len(paths)==189,'path census')
B=[]
for k in (5,6):
 for comb in itertools.combinations(range(21),k):
  b=sum(1<<i for i in comb)
  if all(b&p for p in paths):
   ck(k==6,'minimum size');B.append(b)
ck(len(B)==63,'independent blocker census')
rows=json.loads((P/'result.json').read_text())['map'];ck(set(B)=={r['blocker'] for r in rows},'producer universe')
Nout={c:{v for u,v in E if u==c} for c in V};Nin={c:{u for u,v in E if v==c} for c in V}
def boundary(b):
 parity=collections.Counter(v for i,e in enumerate(E) if b>>i&1 for v in e)
 return frozenset(v for v,n in parity.items() if n%2)
def code(b,phase=10):
 remain=set(e for i,e in enumerate(E) if not b>>i&1)
 sources=V-{v for u,v in remain};sinks=V-{u for u,v in remain}
 if len(sources)==2 or len(sinks)==2:
  S=sources if len(sources)==2 else sinks
  u,v=next((u,v) for u,v in E if {u,v}==S);k={1:0,2:1,4:2}[(v-u)%7]
  return ((9*u+phase*16**k) if len(sources)==2 else (9*v-phase*16**k))%63
 if sources&sinks:return 9*next(iter(sources&sinks))
 for c in V:
  if sum({v for u,v in remain if u==x}=={c} for x in V)==3:return (9*c+21)%63
  if sum({u for u,v in remain if v==x}=={c} for x in V)==3:return (9*c+42)%63
 raise RuntimeError('classification')
ex={b:code(b) for b in B};ck(set(ex.values())==set(range(63)),'bijection')
ck(all(ex[r['blocker']]==r['exponent'] and sum(1<<x for x in boundary(r['blocker']))==r['boundary'] for r in rows),'full saved map')
rho={}
for c in V:
 rho[9*c%21]=frozenset(V-{c});rho[(9*c+7)%21]=frozenset({c}|Nin[c]);rho[(9*c-7)%21]=frozenset({c}|Nout[c])
ck(len(set(rho.values()))==21,'decoder bijection')
for b in B:ck(boundary(b)==rho[ex[b]%21],'decoder')
def act(b,a,t,dual):
 es=[((a*u+t)%7,(a*v+t)%7) for i,(u,v) in enumerate(E) if b>>i&1]
 if dual:es=[(-v%7,-u%7) for u,v in es]
 return mask(es)
for b,a,t,d in itertools.product(B,(1,2,4),range(7),(0,1)):
 target=act(b,a,t,d);k={1:0,2:1,4:2}[a]
 ck(ex[target]==((-1)**d*(16**k*ex[b]+9*t))%63,'G42 action')
# Alternative field implementation: carryless convolution followed by long division.
def mul(a,b):
 z=0
 for i in range(6):
  for j in range(6):
   if (a>>i&1) and (b>>j&1):z^=1<<(i+j)
 for i in range(10,5,-1):
  if z>>i&1:z^=67<<(i-6)
 return z
def power(a,n):
 z=1
 for _ in range(n):z=mul(z,a)
 return z
g=[power(2,e) for e in range(63)];ck(len(set(g))==63 and mul(g[-1],2)==1,'primitive')
log={x:i for i,x in enumerate(g)};F={x for x in range(64) if power(x,4)==x};ck(len(F)==4,'subfield')
inv={0:0,**{x:g[-log[x]%63] for x in g}}
for x in g:ck(mul(x,inv[x])==1,'inverse')
def trace(x):return x^power(x,4)^power(x,16)
for a in g:
 counts=collections.Counter(inv[x^a]^inv[x] for x in range(64))
 ck(collections.Counter(counts.get(y,0) for y in range(64))=={0:33,2:30,4:1},'derivative histogram')
 ck({x for x in range(64) if inv[x^a]^inv[x]==inv[a]}=={mul(a,f) for f in F},'exceptional quartet')
for b in B:ck({g[ex[c]] for c in B if boundary(c)==boundary(b)}=={mul(g[ex[b]],f) for f in F if f},'native scalar fibers')
lines={frozenset(log[x]%21 for x in g if trace(mul(a,x))==0) for a in g}
ck(len(lines)==21 and all(len(L)==5 for L in lines),'trace lines')
# Independently construct projective lines from two-generator vector spans.
spanlines=set()
for a,b in itertools.combinations(g,2):
 if log[a]%21==log[b]%21:continue
 span={mul(a,s)^mul(b,t) for s in F for t in F};ck(len(span)==16,'span dimension')
 spanlines.add(frozenset(log[x]%21 for x in span if x))
ck(lines==spanlines,'trace lines equal geometric lines')
ck(not any(frozenset(-e%21 for e in L) in lines for L in lines),'bare inverse hostile')
def incident(p,l):return trace(g[(p-l)%63])==0
flags={(p,l) for p,l in itertools.product(range(21),repeat=2) if incident(p,l)}
ck(len(flags)==105,'base flags')
for q in range(21):ck(sum(p==q for p,l in flags)==5 and sum(l==q for p,l in flags)==5,'base degrees')
for p,l,a,t,d in itertools.product(range(21),range(21),(1,2,4),range(7),(0,1)):
 k={1:0,2:1,4:2}[a];pp=(16**k*p+9*t)%21;ll=(16**k*l+9*t)%21
 if d:pp,ll=-ll%21,-pp%21
 ck(incident(p,l)==incident(pp,ll),'typed G42 flags')
lift={(b,c) for b,c in itertools.product(B,repeat=2) if incident(ex[b],ex[c])}
ck(len(lift)==945,'lift edges')
for sheet in (0,1):
 neigh={b:frozenset(c for c in B if ((b,c) if sheet==0 else (c,b)) in lift) for b in B}
 ck(all(len(n)==15 for n in neigh.values()),'lift degree')
 for b,c in itertools.product(B,repeat=2):ck((neigh[b]==neigh[c])==(boundary(b)==boundary(c)),'both-sheet twins')
for b,c in itertools.product(B,repeat=2):ck((boundary(b)==boundary(c))==(ex[b]%21==ex[c]%21),'all pair quotient')
for b,c in lift:
 for a,t,d in itertools.product((1,2,4),range(7),(0,1)):
  bb,cc=act(b,a,t,d),act(c,a,t,d)
  ck(((cc,bb) if d else (bb,cc)) in lift,'native typed G42')
ck((63,28679) in lift and (63,455) not in lift and (63&28679).bit_count()==(63&455).bit_count()==3,'intersection hostile')
ck(any(boundary(b)==boundary(c) and code(b,1)%21!=code(c,1)%21 for b,c in itertools.product(B,repeat=2)),'old phase hostile')
absolute=[p for p in range(21) if incident(p,-p%21)];ck(absolute==[3,6,7,12,14],'absolute labels')
print(json.dumps({'status':'PASS','producer_manifest':hashlib.sha256((P/'MANIFEST.json').read_bytes()).hexdigest(),'independent_minimum_blockers':63,'all_native_group_checks':2646,'all_quotient_pairs':3969,'all_typed_group_flag_checks':18522,'base_objects':42,'base_flags':105,'lift_objects':126,'lift_edges':945,'lift_degree':15,'native_flag_symmetry_checks':39690,'absolute_points':absolute,'field_subfield':sorted(F),'inverse_preserved_lines':0,'derivative_directions':63,'both_sheet_twins':'exactly native boundary fibers','hostiles':['old phase 1','bare point inversion','equal native intersections with unequal transported incidence'],'scope':'Paley-7 minimum blockers only; no Hamiltonian-count or all-size transfer'},sort_keys=True,indent=2))
