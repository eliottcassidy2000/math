from itertools import combinations,permutations,product
from collections import defaultdict,Counter,deque
from pathlib import Path
import json
import sympy as s
HERE=Path(__file__).resolve().parent
SRC=HERE.parent/'session29-target-wiring'
def req(b,m):
 if not b:raise RuntimeError(m)
def pairs(n):return list(combinations(range(n),2))
def matrix(mask,n):
 A=[[0]*n for _ in range(n)]
 for k,(i,j) in enumerate(pairs(n)):A[i][j]=(mask>>k)&1;A[j][i]=1-A[i][j]
 return A
def mask(A):return sum(A[i][j]<<k for k,(i,j) in enumerate(pairs(len(A))))
def relabel(A,p):return [[A[i][j] for j in p] for i in p]
def dual(A):return list(map(list,zip(*A)))
def canonical(A,roots):
 n=len(A); fixed=set(roots); free=[i for i in range(n) if i not in fixed];vals=[]
 for q in permutations(free):
  p=list(range(n))
  for j,i in enumerate(free):p[i]=q[j]
  vals.append((mask(relabel(A,p)),p))
 v=min(x[0] for x in vals)
 return v,sum(x[0]==v for x in vals)
def hp(A):
 n=len(A);out=[sum(A[i][j]<<j for j in range(n)) for i in range(n)]
 D=[[0]*n for _ in range(1<<n)]
 for i in range(n):D[1<<i][i]=1
 for q in range(1,1<<n):
  for i in range(n):
   if not D[q][i]:continue
   rest=out[i]&~q
   while rest:
    bit=rest&-rest;rest-=bit;j=bit.bit_length()-1;D[q|bit][j]+=D[q][i]
 return sum(D[-1])
def source_tiles(seed):
 A=matrix(seed,5)
 return sorted({mask(relabel(A,p)) for p in permutations(range(5)) if all(A[p[i]][p[i+1]] for i in range(4))})
alpha={t:typ for typ,seed in [('A',40),('B',41)] for t in source_tiles(seed)}
req(len(alpha)==18,'alphabet size')
flip=1023^sum(1<<pairs(5).index((i,i+1)) for i in range(4))
blocks={t:matrix(t^flip,5) for t in alpha}
mirror={t:mask(relabel(dual(matrix(t,5)),range(4,-1,-1))) for t in alpha}
cap={};entry={};bi={};capgroups=defaultdict(list);bigroups=defaultdict(list)
for t,A in blocks.items():
 cap[t]=canonical(A,[4]);entry[t]=canonical(relabel(dual(A),[4,1,2,3,0]),[4]);bi[t]=canonical(A,[0,4])
 capgroups[cap[t][0]].append(t);bigroups[bi[t][0]].append(t)
req(len(capgroups)==13 and len(bigroups)==17,'root class dimensions')
mirbi={b:bi[mirror[ts[0]]][0] for b,ts in bigroups.items()}
req(sum(k==v for k,v in mirbi.items())==5,'fixed bire-root classes')
req([sorted(ts) for b,ts in bigroups.items() if mirbi[b]==b and not any(mirror[t]==t for t in ts)]==[[755,949]],'invisible blue center')
def seam(t,u):
 X,Y=blocks[t],blocks[u];A=[[0]*10 for _ in range(10)]
 for i in range(10):
  for j in range(10):
   A[i][j]=(X if i<5 else Y)[i%5][j%5] if i//5==j//5 else int(i>j)
 A[4][5]=1;A[5][4]=0
 return A
def onecuts(A):
 n=len(A);out=[]
 for q in range(1,(1<<n)-1):
  count=sum(A[i][j] for i in range(n) if q>>i&1 for j in range(n) if not q>>j&1)
  if count==1:out.append(q)
 return out
W=json.loads((SRC/'wiring-results.json').read_text());R=W['witness_rows'];observed=defaultdict(lambda:Counter());H={};edges=[];knownids={tuple(k):i+3 for i,k in enumerate(sorted({tuple(r['merged_key']) for r in R}))};sources={'AA':0,'AB/BA':1,'BB':2}
ordinary=set();autos=Counter();strong_bridge_hist=Counter()
for rec in R:
 t,u=rec['tiles'];A=seam(t,u);q=cap[t][0];r=entry[u][0];key=tuple(sorted((q,r)));ordinary.add((q,r));source=alpha[t]+alpha[u];source=source if source in ['AA','BB'] else 'AB/BA';blue=u==mirror[t]
 req(mask(A)==rec['target'],'literal target')
 req(key==tuple(rec['merged_key']) and (q,r)==tuple(rec['ordinary_key']),'root keys')
 req(source==rec['source'] and blue==rec['blue'],'source/color')
 val=hp(A);req(val==rec['H'],'full target Hamiltonian count');H[key]=val
 aut=cap[t][1]*entry[u][1];req(aut==rec['target_aut_order'],'factored automorphism order');autos[aut]+=1
 cs=onecuts(A);mid=[x for x in cs if x.bit_count()==5];req(mid==[31],'intrinsic balanced cut')
 interior=[x for x in cs if 2<=x.bit_count()<=8]
 req(all((x&y)==x or (x&y)==y for x,y in combinations(interior,2)),'nested actual cuts')
 C=matrix(next(x['canonical_target_mask'] for x in W['incidence_table'] if tuple(x['root_pair'])==key),10);p=rec['old_to_canonical_vertex'];req(sorted(p)==list(range(10)),'witness bijection');AA=dual(A) if rec['converse_used'] else A
 req(all(AA[i][j]==C[p[i]][p[j]] for i in range(10) for j in range(10)),'literal witness')
 observed[key][source,blue]+=1;edges.append((sources[source],knownids[key],blue))
req(len(R)==324 and len(ordinary)==169 and len(observed)==91,'target census')
for row in W['incidence_table']:
 key=tuple(row['root_pair']);q,r=key
 typq=Counter(alpha[t] for t in capgroups[q]);typr=Counter(alpha[t] for t in capgroups[r]);factor=1 if q==r else 2
 exp={'AA':factor*typq['A']*typr['A'],'BB':factor*typq['B']*typr['B'],'AB/BA':factor*(typq['A']*typr['B']+typq['B']*typr['A'])}
 for src in sources:
  req(sum(observed[key][src,b] for b in [False,True])==exp[src],'weighted symmetric-square incidence')
  for blue in [False,True]:req(observed[key][src,blue]==row['incidences'][src]['blue' if blue else 'black'],'published incidence')
# Independently recompute every multigraph bridge.
def reachable(start,drop=None,black=False):
 seen={start};todo=[start]
 while todo:
  a=todo.pop()
  for i,(u,v,c) in enumerate(edges):
   if i==drop or (black and c):continue
   if a in (u,v):
    b=v if a==u else u
    if b not in seen:seen.add(b);todo.append(b)
 return seen
bridges=[i for i,(a,b,c) in enumerate(edges) if b not in reachable(a,i)]
req(len(bridges)==8 and all(edges[i][2] for i in bridges),'restricted bridge census')
req(len(reachable(0))==94,'boundary connected')
blackverts={x for a,b,c in edges if not c for x in [a,b]};req(reachable(0,black=True)==blackverts,'black boundary connected')
cycledata=json.loads((SRC/'cycle-results.json').read_text())
for cert in cycledata['cycle_certificates']:
 i=cert['edge'];a,b,_=edges[i];cur=a
 for j in cert['alternative_path']:
  req(i!=j and cur in edges[j][:2],'alternative path edge');u,v,_=edges[j];cur=v if cur==u else u
 req(cur==b,'alternative path terminates')
req(len(cycledata['cycle_certificates'])==316,'cycle certificate count')
# Scalar matrices are already independently certified in the earlier audit;
# here check that they descend to the newly reconstructed classes.
M={int(k):s.Matrix(v) for k,v in json.loads((HERE.parent/'session29-profiles/compressed-matrices.json').read_text()).items()};J=s.diag(1,-1,1);e=s.Matrix([1,0,0])
for ts in capgroups.values():req(all(M[t]*e==M[ts[0]]*e for t in ts),'cap observer descends')
for ts in bigroups.values():req(all(M[t]==M[ts[0]] for t in ts),'interior observer descends')
B=sum((M[ts[0]] for ts in bigroups.values()),s.zeros(3));v=sum((M[ts[0]]*e for ts in capgroups.values()),s.zeros(3,1));z=s.symbols('z');gf=s.cancel(z*z*(v.T*J*(s.eye(3)-z*B).inv()*v)[0]);expected=z*z*(93313+3359832*z+4057812*z*z)/(1-823*z-37140*z*z-45180*z**3)
req(s.cancel(gf-expected)==0,'ordinary distinct-target H GF')
req(sum(H.values())==int(json.loads((SRC/'all-length-results.json').read_text())['merged_H_r2']),'merged H mass')
req(sum(H[k]*(1 if k[0]==k[1] else 2) for k in H)==93313,'ordinary H mass')
# Independent nine-coordinate second-moment transfer, rather than six monomials.
L=sum((s.kronecker_product(M[ts[0]],M[ts[0]]) for ts in bigroups.values()),s.zeros(9));v0=sum((s.kronecker_product(M[ts[0]]*e,M[ts[0]]*e) for ts in capgroups.values()),s.zeros(9,1));fixedM=sum((M[ts[0]] for b,ts in bigroups.items() if mirbi[b]==b),s.zeros(3));ce=s.Matrix([[J[i,j] for i in range(3) for j in range(3)]]);JO=J*fixedM;co=s.Matrix([[JO[i,j] for i in range(3) for j in range(3)]]);ALL=json.loads((SRC/'all-length-results.json').read_text());w=s.symbols('w');E=s.sympify(ALL['fixed_even_H_GF_in_w']);O=s.sympify(ALL['fixed_odd_H_GF_in_w'])
# Verify rational identities by polynomial matrix annihilation, without costly 9x9 inversion.
for cost,fn in [(ce,E),(co,O)]:
 num,den=s.fraction(s.cancel(fn));coeffs=s.Poly(den,w);seq=[];state=v0
 for k in range(18):seq.append((cost*state)[0]);state=L*state
 for k in range(18):
  lhs=sum(coeffs.coeff_monomial(w**j)*seq[k-j] for j in range(min(k,coeffs.degree())+1));rhs=s.Poly(num,w).coeff_monomial(w**k)
  req(lhs==rhs,'nine-state fixed-target rational moments')
 # Denominator annihilates all moments after the numerator transient, via CH-sized state witness.
 d=coeffs.degree();operator=sum((coeffs.coeff_monomial(w**j)*L**(d-j) for j in range(d+1)),s.zeros(9))
 req(cost*operator==s.zeros(1,9),'exact nine-state denominator annihilator')
result={'status':'PASS','alphabet':{k:source_tiles(v) for k,v in [('A',40),('B',41)]},'a':13,'b':17,'f':5,'full_target_DP_checks':324,'literal_isomorphism_witnesses':324,'ordinary_targets_r2':169,'merged_targets_r2':91,'r2_blue':sum(c for a,b,c in edges),'restricted_bridges':len(bridges),'globally_certified_nonbridges':316,'black_boundary_vertices':len(blackverts),'black_boundary_edges':sum(not c for a,b,c in edges),'root_automorphism_products':dict(autos),'B':B.tolist(),'v':list(v),'ordinary_H_GF':str(gf),'ordinary_H_r2':93313,'merged_H_r2':sum(H.values()),'fixed_target_moment_representation':'independent nine-coordinate Kronecker with exact denominator annihilation'}
print(json.dumps(result,indent=2,default=int))
