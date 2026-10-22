"""Independent flat-class audit: symbolic arithmetic and two-round triangle flips."""
import json, collections
from pathlib import Path
import numpy as np
import sympy as s
z,u=s.symbols('z u');D=lambda f:s.expand(z*s.diff(f,z))
g=z*z*u+(z**4-z*z)*u*u;h=-z*z*u+(z*z+z**4)*u*u
fig=s.series((2*g*h-s.Rational(3,2)*g*g)/u,u,0,3).removeO()
theta=(z**4+z**6/s.Integer(3))*u*u
A=s.expand(fig+theta-z*z*u+4*z**4*u*u)
assert s.expand(A-(-(z*z+s.Rational(7,2)*z**4)*u+(12*z**4-s.Rational(8,3)*z**6)*u*u))==0
X=z*z*(1-u)+z**6*u*u/3
c=[None]+[X]
for j in range(1,5):
 if j==1:c[1]=D(X)
 else:c.append(D(c[-1]))
b1=c[1];b2=c[2]+c[1]**2;b3=c[3]+3*c[1]*c[2]+c[1]**3;b4=c[4]+4*c[1]*c[3]+3*c[2]**2+6*c[1]**2*c[2]+c[1]**4
p1=(b2-b1)/2;p2=(3*b4-2*b3-3*b2+2*b1)/24
shift=(D(D(A))+(2*c[1]-1)*D(A))/2
Y1=s.series(A+u*p1,u,0,3).removeO()
Y2=s.series((7*z**4+26*z**6)*u*u+A*A/2+u*(A*p1+shift)+u*u*p2,u,0,3).removeO()
k=s.series((1-u+s.Rational(5,6)*u*u)*(1+Y1.subs(z,1)+Y2.subs(z,1)),u,0,3).removeO()
assert k==1-s.Rational(5,2)*u+s.Rational(65,8)*u*u
# Independent tensor evaluation, using generic scalar contractions and parity.
cores=json.loads((Path(__file__).resolve().parents[1]/'session30-decorated-engine/symbolic-kernels.json').read_text())
polys=[]
for core in cores:
 if core['weight']!=2:continue
 edges=core['edges'];degrees=collections.Counter(x for edge in edges for x in edge);leaf={x for x in degrees if degrees[x]==1}
 powers=collections.Counter();colors=collections.Counter()
 for a,b in edges:
  if a in leaf:colors[b]+=1
  elif b in leaf:colors[a]+=1
  elif a!=b:powers[tuple(sorted((a,b)))]+=1
 poly=0
 if all(v%2==0 for v in powers.values()) and all(v%2==0 for v in colors.values()):
  for term in core['terms']:
   mon=s.Rational(term['numerator'],core['denominator'])
   for (a,b),j in zip(edges,term['indices']):
    if a==b:mon*=(-z*z,-z*z,z*z,-z*z)[j]
    elif a in leaf:mon*=(-z,z,0,0)[j]
    else:mon*=(z,-z,0,0)[j]
   poly+=mon
 polys.append(s.expand(poly))
assert s.expand(sum(polys))==7*z**4+26*z**6

def initial(q):
 squares={a*a%q for a in range(1,q)}
 return np.array([[0 if i==j else 1 if (i-j)%q in squares else -1 for j in range(q-1)] for i in range(q-1)],dtype=np.int64)
def perturb(S,seed):
 T=S.copy();rng=np.random.default_rng(seed);packs=[]
 for round in range(2):
  remaining=list(rng.permutation(len(T)));triples=[]
  while True:
   chosen=None
   for i,a in enumerate(remaining):
    for j,b in enumerate(remaining[i+1:],i+1):
     for c in remaining[j+1:]:
      if T[a,b]==T[b,c]==T[c,a]:chosen=(a,b,c);break
     if chosen:break
    if chosen:break
   if chosen is None:break
   for a,b in ((chosen[0],chosen[1]),(chosen[1],chosen[2]),(chosen[2],chosen[0])):T[a,b]*=-1;T[b,a]*=-1
   triples.append(chosen);remaining=[a for a in remaining if a not in chosen]
  packs.append(len(triples))
 return T,packs
rows=[]
for q in (11,23,47,83,131,199):
 S=initial(q);T,packs=perturb(S,q);n=len(S);d=T.sum(1);E=T-S;B=T@T+(n-1)*np.eye(n,dtype=np.int64)
 assert np.array_equal(d,S.sum(1));assert np.max(abs(E).sum(1))<=8
 assert np.max(abs(B))<=34;assert np.max(abs(T@d))<=9
 for w in (1,.7+.4j,1.2j):
  Q=w*T/n;Q2=Q@Q;I=np.eye(n);G=-np.linalg.solve(I+Q2,Q2);H=np.linalg.solve(I-Q2,Q2);O=np.linalg.solve(I-Q2,Q)
  diag=np.diag(H-G);delta=d/n
  exact=np.sum(2*np.diag(G)*np.diag(H)-1.5*np.diag(G)**2)+np.sum(G*(H*H-O*O-2*G*G/3))-diag@G@diag-w*w*delta@np.linalg.solve(I+Q2,delta)-2*diag@G@np.ones(n)
  target=-(w*w+3.5*w**4)/n+(12*w**4-8*w**6/3)/n**2
  residual=abs(n**3*(exact-target))
  assert residual<1000,(n,w,residual)
  rows.append({'n':n,'z':str(w),'pack_sizes':packs,'max_B':int(abs(B).max()),'max_Sd':int(abs(T@d).max()),'scaled_A1_remainder_abs':float(residual)})
  print(rows[-1],flush=True)
Path(__file__).with_name('receipt.json').write_text(json.dumps({'A1':str(A),'A2':str(sum(polys)),'kappa2_over_e':str(k),'nonzero_core_polynomials':[str(p) for p in polys if p!=0],'controls':rows},indent=2)+'\n')
