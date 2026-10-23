import itertools as I,json
from collections import Counter,defaultdict
import sympy as S
checks=0
def ck(v,msg):
 global checks
 checks+=1
 if not v:raise ValueError(msg)
# S3 full six-edge integral, without root fixing or assumed face independence.
G=list(I.permutations(range(3)));one=G.index((0,1,2));n=6
M=[[G.index(tuple(a[b[i]] for i in range(3))) for b in G] for a in G]
inv=[next(j for j in range(6) if M[i][j]==one) for i in range(6)]
cl=[0 if g==(0,1,2) else 1 if sum(g[i]==i for i in range(3))==1 else 2 for g in G]
terms=Counter();flat=0
for a,b,c,d,e,f in I.product(range(6),repeat=6):
 faces=[M[M[a][d]][inv[b]],M[M[a][e]][inv[c]],M[M[b][f]][inv[c]],M[M[d][f]][inv[e]]]
 powers=tuple(sum(cl[x]==j for x in faces) for j in range(3));terms[powers]+=1
 if powers==(4,0,0):flat+=1
x,y,z=S.symbols('x y z');variables=(x,y,z)
P=S.expand(sum(N*S.prod(v**p for v,p in zip(variables,pow)) for pow,N in terms.items())/6**6)
coeff=[(x+3*y+2*z)/6,(x-3*y+2*z)/6,(x-z)/6]
ck(S.expand(P-sum(d*d*a**4 for d,a in zip((1,1,2),coeff)))==0,'full S3 symbolic class weights')
ck(flat==216,'full S3 flat configurations')
ck(S.expand(P.subs({y:1,z:1})*216-(x**4+30*x*x+80*x+105))==0,'rooted identity polynomial from full integration')
q=S.symbols('q')
Q=S.cancel(((x+q-1)**4+(q-1)*(x-1)**4)/q)
ck(S.expand(Q-(x**4+6*(q-1)*x*x+4*(q-1)*(q-2)*x+(q-1)*(q*q-3*q+3)))==0,'arbitrary q polynomial')
# Independent nonabelian D8 control, elements r^i s^j.
D=list(I.product(range(4),range(2)))
def mul(a,b):return ((a[0]+(-1)**a[1]*b[0])%4,(a[1]+b[1])%2)
def inverse(a):return next(b for b in D if mul(a,b)==(0,0))
hist=Counter()
for a,b,c in I.product(D,repeat=3):
 faces=(a,b,c,mul(mul(a,b),inverse(c)))
 hist[sum(f==(0,0) for f in faces)]+=1
ck(S.expand(sum(N*x**k for k,N in hist.items())-Q.subs(q,8))==0,'D8 rooted identity polynomial')
# Derive rainbow kernel directly as nullspace of geometric cycles, no Fano code input.
edges=list(I.combinations(range(7),2));ix={e:i for i,e in enumerate(edges)}
lines=[{i,(i+1)%7,(i+3)%7} for i in range(7)]
ell={e:next(i for i,l in enumerate(lines) if set(e)<=l) for e in edges}
rows=set()
for a,b,c,d in I.permutations(range(7),4):
 es=[tuple(sorted(e)) for e in ((a,b),(b,c),(c,d),(d,a))]
 if len({ell[e] for e in es})==4:rows.add(sum(1<<ix[e] for e in es))
piv={}
for r in rows:
 for p in sorted(piv,reverse=True):
  if r>>p&1:r^=piv[p]
 if r:piv[r.bit_length()-1]=r
free=[i for i in range(21) if i not in piv];basis=[]
for i in free:
 v=1<<i
 for p in sorted(piv):
  if (v&piv[p]).bit_count()%2:v^=1<<p
 basis.append(v)
ker={0}
for b in basis:ker|={v^b for v in ker}
ck(len(ker)==1024 and all((v&r).bit_count()%2==0 for v in ker for r in rows),'direct nullspace')
cut=[sum(1<<j for j,e in enumerate(edges) if i in e) for i in range(1,7)]
def rootfix(v):
 w=v
 for i in range(1,7):
  if v>>ix[(0,i)]&1:w^=cut[i-1]
 return w
cosets=defaultdict(list)
for v in ker:cosets[rootfix(v)].append(v)
ck(len(cosets)==16 and all(len(v)==64 for v in cosets.values()),'nullspace gauge cosets')
summary=Counter((min(v.bit_count() for v in vs),sum(v.bit_count()==min(w.bit_count() for w in vs) for v in vs)) for vs in cosets.values())
ck(summary=={(0,1):1,(3,1):7,(6,3):7,(9,35):1},'logical minima and multiplicities')
weights=Counter(v.bit_count() for v in ker)
ck(weights=={0:1,3:7,6:28,7:42,8:105,9:182,10:147,11:147,12:182,13:105,14:42,15:28,18:7,21:1},'kernel weight enumerator')
# Coordinate/scaling identities are algebraic checks, not spectral existence certificates.
u,v,g=S.symbols('u v g',positive=True)
ck(S.expand(g*((u+v)/S.sqrt(2)-(v-u)/S.sqrt(2))**4)==4*g*u**4,'quartic orthogonal coordinates')
ck(S.simplify(g*(g**(-S.Rational(1,6)))**4-g**S.Rational(1,3))==0,'quartic potential scaling')
ck(S.simplify((g**(-S.Rational(1,6)))**-2-g**S.Rational(1,3))==0,'kinetic scaling')
# Gauge-invariant single-edge kinetics: exact short virtual paths.
def energy(v):return S.Rational(sum((v&r).bit_count()%2 for r in rows),4)
ck(all(energy(1<<i)==1 for i in range(21)),'all single-edge denominators J')
second=Counter();third=Counter()
for i,j in I.product(range(21),repeat=2):
 v=(1<<i)^(1<<j)
 if rootfix(v) in cosets:second[rootfix(v)]+=1/energy(1<<i)
ck(dict(second)=={0:21},'second-order scalar coefficient')
for i,j,k in I.product(range(21),repeat=3):
 mid=(1<<i)^(1<<j)
 if rootfix(mid) in cosets:continue
 end=mid^(1<<k)
 if rootfix(end) in cosets:
  third[rootfix(end)]+=1/(energy(1<<i)*energy(mid))
ck(len(third)==7 and set(third.values())=={6},'third-order seven matching transitions coefficient6')
ck(all(min(v.bit_count() for v in cosets[k])==3 for k in third),'third-order graph equals minimal matching sectors')
logical=list(cosets);Adj=S.Matrix([[int((a^b) in third) for b in logical] for a in logical])
ck(Adj.eigenvals()=={S.Integer(7):1,S.Integer(1):7,S.Integer(-1):7,S.Integer(-7):1},'edge-kinetic crown adjacency spectrum')
fourth=Counter();paircost=Counter(energy((1<<i)^(1<<j)) for i,j in I.product(range(21),repeat=2) if i!=j)
for i,j,k,l in I.product(range(21),repeat=4):
 s1=1<<i;s2=s1^(1<<j);s3=s2^(1<<k);end=s3^(1<<l)
 if end not in ker or s2 in ker or s3 in ker:continue
 fourth[rootfix(end)]+=1/(energy(s1)*energy(s2)*energy(s3))
ck(dict(fourth)=={0:518},'all four-step Q-intermediate paths')
ck(paircost=={S.Integer(1):42,S.Rational(3,2):168,S.Integer(2):210},'half-integer pair costs')
ck(441-fourth[0]==-77,'fourth-order scalar after energy dependence')
print(json.dumps({'status':'PASS','checks':checks,'edge_kinetic_fourth_path_sum':str(fourth[0]),'edge_kinetic_fourth_order_scalar':str(441-fourth[0]),'edge_kinetic_second_order':21,'edge_kinetic_third_order_transition_count':len(third),'edge_kinetic_third_order_each':6,'S3_full_symbolic_integral':str(P),'S3_full_flat_count':flat,'D8_rooted_histogram':dict(hist),'kernel_weight_enumerator':dict(sorted(weights.items())),'logical_distance_and_minimizer_multiplicities':[[d,m,n] for (d,m),n in sorted(summary.items())]},indent=2))
