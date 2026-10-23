import itertools as I, json
from collections import Counter
import sympy as S
checks=0
def ck(x,why):
 global checks
 checks+=1
 if not x: raise RuntimeError(why)
def parity(x):return x.bit_count()%2
def syn(x):
 r=0
 for k in range(7):
  if x>>k&1:r^=k+1
 return r
C=[x for x in range(128) if not syn(x)]
# Characters derived directly by their 16 values on C, no dual-code coset construction.
characters={tuple((-1)**parity(q&c) for c in C):q for q in range(128)}
ck(len(characters)==16,'16 different code characters')
reps=[0]+[1<<(s-1) for s in range(1,8)]
t,z=S.symbols('t z'); classes=Counter(); polys={}; traces={}
for chi,q in characters.items():
 # f(rep_s+c)=chi(c) f(rep_s). Each flip changes s to s XOR a.
 A=S.zeros(8)
 for s in range(8):
  for a in range(1,8):
   dest=s^a;c=reps[s]^(1<<(a-1))^reps[dest]
   ck(syn(c)==0,'transition cocycle in C')
   A[s,dest]+=(-1)**parity(q&c)
 ck(A==A.T,'signed quotient symmetric')
 H=t*(7*S.eye(8)-A)+S.diag(0,*([1]*7))
 p=S.factor(H.charpoly(z).as_expr())
 key=str(p);classes[key]+=1;polys[key]=p
 # Independent Schur expansion at distinguished zero-potential origin.
 # H00=7t; offdiag=-t v; lower=I+t(7I-A1).
 v=A[1:,0];Q=A[1:,1:]
 # E=7t - t² vT [I+t(7I-Q)-E I]^-1 v.
 # E=7t-7t²-t³ vT Qv+..., giving cubic coefficient -vT Qv.
 ck((v.T*v)[0]==7,'universal second order coefficient')
 cubic=-(v.T*Q*v)[0]
 traces[key]=int(cubic)
 ck(S.expand(p.subs(z,7*t-7*t*t+cubic*t**3)).series(t,0,4).removeO()==0,'Schur cubic agrees with exact polynomial')
ck(sorted((traces[k],n) for k,n in classes.items())==[(-42,1),(-6,7),(6,7),(42,1)],'four low branches and degeneracies')
# Independent simplex audit via GF2 boundary masks, avoiding hyperbolic evaluation.
n=5;edges=list(I.combinations(range(n),2)); ix={e:i for i,e in enumerate(edges)};tri=list(I.combinations(range(n),3))
face=[sum(1<<ix[e] for e in I.combinations(f,2)) for f in tri]
closed=Counter();corr=Counter();c4=sum(1<<ix[tuple(sorted(e))] for e in [(0,1),(1,2),(2,3),(3,0)])
for m in range(1<<len(face)):
 b=0
 for i,f in enumerate(face):
  if m>>i&1:b^=f
 if not b:closed[m.bit_count()]+=1
 if b==c4:corr[m.bit_count()]+=1
ck(dict(closed)=={0:1,4:5,6:10},'n5 closed surface polynomial')
ck(dict(corr)=={2:2,4:5,6:8,8:1},'n5 C4 numerator')
fluxsets=set()
for x in range(1<<len(edges)):
 flux=tuple(parity(x&f) for f in face);fluxsets.add(flux)
 for tetra in I.combinations(range(n),4):
  ck(sum(flux[tri.index(f)] for f in I.combinations(tetra,3))%2==0,'tetrahedron Bianchi')
 rootfix=x
 for i in range(1,n):
  if x>>ix[(0,i)]&1:
   rootfix^=sum(1<<j for j,e in enumerate(edges) if i in e)
 ck(all(not(rootfix>>ix[(0,i)]&1) for i in range(1,n)),'root gauge')
 ck(tuple(parity(rootfix&f) for f in face)==flux,'curvature unchanged')
ck(len(fluxsets)==64,'curvature orbit count')
# Fano restriction, baseline signs recomputed and direct Hamiltonian DP.
n=7;edges=list(I.combinations(range(n),2));ix={e:i for i,e in enumerate(edges)}
lines=[{i,(i+1)%7,(i+3)%7} for i in range(7)]
ell={e:next(i for i,l in enumerate(lines) if set(e)<=l) for e in edges}
baseline=sum(1<<ix[e] for e in edges if (e[1]-e[0])%7 not in {1,2,4})
cycles=set()
for a,b,c,d in I.permutations(range(7),4):
 cycles.add(sum(1<<ix[tuple(sorted(e))] for e in [(a,b),(b,c),(c,d),(d,a)]))
rainbow=[];supports=Counter()
for cyc in cycles:
 labels=[ell[e] for e in edges if cyc>>ix[e]&1]
 if len(set(labels))==4:
  rainbow.append(cyc);supports[sum(1<<l for l in labels)]+=1
ck(len(cycles)==105 and len(rainbow)==21,'rainbow incidence')
ck(sorted(supports.values())==[3]*7,'seven supports multiplicity three')
# Incidence alone establishes support equals nonzero dual affine hyperplanes in given convention.
a=[4,1,3,7,6,5,2]
ck(set(supports)=={sum(parity(b&v)<<i for i,v in enumerate(a)) for b in range(1,8)},'line normal convention')
def hp(x):
 dp={(1<<v,v):1 for v in range(7)}
 for mask in range(1,128):
  for last in range(7):
   val=dp.get((mask,last),0)
   if not val:continue
   for v in range(7):
    if mask>>v&1:continue
    e=tuple(sorted((last,v)));neg=x>>ix[e]&1
    if (last<v) != bool(neg):
     key=(mask|1<<v,v);dp[key]=dp.get(key,0)+val
 return sum(dp.get((127,v),0) for v in range(7))
hist=Counter();rhist=Counter()
for m in range(128):
 rel=sum(1<<ix[e] for e in edges if m>>ell[e]&1);x=baseline^rel
 R=sum((-1)**parity(rel&c) for c in rainbow);H=hp(x)
 ck(4*H==693+3*R,'normalized Wilson Hamiltonian identity')
 hist[H]+=1;rhist[R]+=1
ck(hist=={189:16,171:112},'Fano histogram')
cut=sum(1<<ix[e] for e in edges if 0 in e)
ck(hp(baseline)==189 and hp(baseline^cut)==135,'switching hostile H')
ck(all(parity(c&cut)==0 for c in cycles),'switching preserves Wilson loops')
# Direct binary rank is an independent transversality certificate.
def rank(rows):
 piv={}
 for r in rows:
  while r:
   k=r.bit_length()-1
   if k in piv:r^=piv[k]
   else:piv[k]=r;break
 return len(piv)
linebits=[sum(1<<ix[e] for e in edges if ell[e]==i) for i in range(7)]
cuts=[sum(1<<ix[e] for e in edges if i in e) for i in range(6)]
ck(rank(linebits+cuts)==13,'line slice transverse to cuts')
# Full rainbow energy by enumerating the 11-dimensional image, not 32768 gauge representatives.
cols=[sum(((cyc>>j)&1)<<k for k,cyc in enumerate(rainbow)) for j in range(21)]
image={0}
for col in cols:image|={v^col for v in image}
ck(len(image)==2048 and rank(rainbow)==11,'rainbow evaluation rank')
fullhist=Counter()
for v in image:fullhist[21-2*v.bit_count()]+=1024
ck(dict(fullhist)=={-7:282624,-3:637952,1:688128,5:322560,9:143360,13:21504,21:1024},'full rainbow histogram')
linecode=[]
for m in range(128):
 syndrome=0
 for i,normal in enumerate(a):
  if m>>i&1:syndrome^=normal
 if not syndrome:
  linecode.append(sum(1<<ix[e] for e in edges if m>>ell[e]&1))
ck(rank(cuts+linecode)==10,'candidate full kernel dimension')
ck(all(not parity(v&cyc) for v in cuts+linecode for cyc in rainbow),'kernel contains cuts plus line code')
print(json.dumps({'status':'PASS','checks':checks,'method':'real-space signed syndrome quotients and Schur complement; direct independent incidence/DP','full_sector_charpolys':[{'polynomial':k,'multiplicity':v,'cubic':traces[k]} for k,v in classes.items()],'n5_closed':dict(closed),'n5_c4_numerator':dict(corr),'Fano_H_histogram':dict(hist),'Fano_R_histogram':dict(rhist),'full_rainbow_histogram':dict(fullhist),'rainbow_rank':11,'rainbow_kernel_dimension':10},indent=2))
