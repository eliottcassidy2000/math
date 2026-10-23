import itertools as I,json
from fractions import Fraction as F
from collections import Counter
import sympy as S
checks=0
def ck(x,m):
 global checks
 checks+=1
 if not x:raise RuntimeError(m)
G=list(I.permutations(range(3)));idx={x:i for i,x in enumerate(G)};e=idx[(0,1,2)]
M=[[idx[tuple(x[y[k]] for k in range(3))] for y in G] for x in G]
inv=[next(y for y in range(6) if M[x][y]==e) for x in range(6)]
sign=[(-1)**sum(x[i]>x[j] for i in range(3) for j in range(i+1,3)) for x in G]
std=[sum(x[i]==i for i in range(3))-1 for x in G]
def conv(f,g):return [sum(f[x]*g[M[inv[x]][y]] for x in range(6))/6 for y in range(6)]
# Distribution of commutator products, retaining Haar density normalization.
comm=[F(0)]*6
for a,b in I.product(range(6),repeat=2):comm[M[M[M[a][b]][inv[a]]][inv[b]]]+=F(1,6)
ck(sum(comm)==6,'commutator Haar density normalization')
Ks=[];surface=[]
for eps in [F(1,4),F(-1,4)]:
 w=[1+eps*sign[i]+F(1,4)*std[i] for i in range(6)]
 ck(min(w)>0 and sum(w)==6,'strictly positive normalized density')
 K=S.Matrix([[S.Rational(w[M[inv[x]][y]])/6 for y in range(6)] for x in range(6)])
 ck(K==K.T and K*S.ones(6,1)==S.ones(6,1),'selfadjoint stochastic kernel')
 # Polynomial annihilation and trace recover whole spectrum without eigenvectors.
 ck((K-S.eye(6))*(K-S.Rational(eps)*S.eye(6))*(K-S.eye(6)/8)==S.zeros(6),'minimal polynomial roots')
 ck(K.det()==S.Rational(eps)/8**4,'determinant')
 sig=S.Matrix(sign);ck((sig.T*K*sig)[0]/6==eps,'normalized reflection form')
 Ks.append(K)
 for faces in [2,4,6]:
  wf=w
  for k in range(1,faces):wf=conv(wf,w)
  cg=[F(0)]*6;cg[e]=6
  for genus in range(4):
   if genus:cg=conv(cg,comm)
   Z=sum(cg[x]*wf[x] for x in range(6))/6
   expect=1+eps**faces+F(2)**(2-2*genus)*F(1,8)**faces
   ck(Z==expect,'genus handle gluing normalization')
   surface.append([str(eps),genus,faces,str(Z)])
# Generator built from central spectral projectors, then verify actual jump rates.
P0=S.ones(6)/6;Ps=S.Matrix(sign)*S.Matrix(sign).T/6;Pt=S.eye(6)-P0-Ps
L=-2*Ps-3*Pt
ck(P0+Ps/4+Pt/8==Ks[0],'generator exponential via orthogonal projectors')
ck(L*S.ones(6,1)==S.zeros(6,1),'zero row sum')
for x,y in I.product(range(6),repeat=2):
 if x==y:continue
 rate=S.Rational(1,3) if sign[M[inv[x]][y]]==-1 else S.Rational(2,3)
 ck(L[x,y]==rate,'all CTMC jump rates divided by log2')
# Independent quartic line-model coefficient from real-space Schur moments.
def parity(x):return x.bit_count()%2
def syn(x):
 r=0
 for k in range(7):
  if x>>k&1:r^=k+1
 return r
C=[x for x in range(128) if syn(x)==0]
chars={tuple((-1)**parity(q&c) for c in C):q for q in range(128)}
rep=[0]+[1<<(s-1) for s in range(1,8)];bands=Counter()
for q in chars.values():
 A=S.zeros(8)
 for s in range(8):
  for a in range(1,8):A[s,s^a]=(-1)**parity(q&(rep[s]^(1<<(a-1))^rep[s^a]))
 v=A[1:,0];Q=A[1:,1:]
 c3=int(-(v.T*Q*v)[0]);c4=int(49-(v.T*Q*Q*v)[0])
 bands[c3,c4]+=1
ck(bands=={(-42,-203):1,(-6,-11):7,(6,-11):7,(42,-203):1},'independent real-space quartic coefficients')
ck((-11)-(-203)==192,'line gap quartic coefficient')
print(json.dumps({'status':'PASS','checks':checks,'line_band_cubic_quartic_multiplicities':[[a,b,n] for (a,b),n in sorted(bands.items())],'minus_kernel_determinant':str(Ks[1].det()),'genus_gluing_tests':surface},indent=2))
