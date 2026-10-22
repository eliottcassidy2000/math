import json,itertools as it
from pathlib import Path
import sympy as sp
import seamed as s
D=Path(__file__).parent
x=json.loads((D/'transfer-results.json').read_text()); letters=[int(t) for t in x['matrices']]
N={}
for t in letters:
 m=x['matrices'][str(t)];v=x['initial'][str(t)]
 # 5-state rowvector -> 4-state columnvector (H,S,E,Q).
 rows=[[sum(m[i][:3]),m[i][0],m[i][1],sum(m[i][3:])] for i in(2,0,1,3)]
 M=sp.Matrix.hstack(sp.Matrix(rows[0]),sp.Matrix([rows[1][j]-rows[0][j] for j in range(4)]),sp.Matrix([rows[2][j]-rows[0][j] for j in range(4)]),sp.Matrix(rows[3]))
 V=sp.Matrix([sum(v[:3]),v[0],v[1],sum(v[3:])])
 N[t]=M.row_join(V).col_join(sp.zeros(1,5))
tau={t:s.mirror(t,5) for t in letters};s.require(set(tau.values())==set(letters),'mirror closedalphabet')
fixed=[t for t in letters if tau[t]==t];s.require(len(fixed)==4,'fixedalphabet4')
def K(Q):return sum((N[tau[t]]*Q*N[t] for t in letters),sp.zeros(5))
F=sum((N[t] for t in fixed),sp.zeros(5));Pe=sp.eye(5);Po=F;even=[];odd=[]
for i in range(61):
 even.append(int(Pe[0,4]));odd.append(int(Po[0,4]));Pe=K(Pe);Po=K(Po)
s.require(even[1]==11522,'r2independentblueHsum')
# Exact rational Berlekamp-Massey, validated with held-out operator coefficients.
from fractions import Fraction
def bm(seq):
 C=[Fraction(1)];B=[Fraction(1)];L=0;m=1;b=Fraction(1)
 for n in range(len(seq)):
  d=Fraction(seq[n])+sum(C[i]*seq[n-i] for i in range(1,L+1))
  if d==0:m+=1;continue
  old=C[:];coef=-d/b
  if len(C)<len(B)+m:C += [Fraction(0)]*(len(B)+m-len(C))
  for j in range(len(B)):C[j+m]+=coef*B[j]
  if 2*L<=n:L=n+1-L;B=old;b=d;m=1
  else:m+=1
 return C[:L+1]
w=sp.symbols('w');rats=[]
for seq in(even,odd):
 c=bm(seq[:40]);d=len(c)-1
 for i in range(40,len(seq)):s.require(sum(c[j]*seq[i-j] for j in range(d+1))==0,'held-outrecurrence')
 den=sum(sp.Rational(v.numerator,v.denominator)*w**i for i,v in enumerate(c))
 num=sum(sum(sp.Rational(c[j].numerator,c[j].denominator)*seq[i-j] for j in range(min(i,d)+1))*w**i for i in range(d))
 rats.append(sp.cancel(num/den))
# These rational formulas are recurrence candidates until operator identity certified.
print('CANDIDATES',rats,flush=True)
# Prove the denominator annihilates the entire matrix orbit: sufficient all-order certificate.
cert=[]
for seed,rat in zip((sp.eye(5),F),rats):
 den=sp.Poly(sp.denom(rat),w);co=list(reversed(den.all_coeffs()));degree=len(co)-1
 powers=[seed]
 for i in range(degree):powers.append(K(powers[-1]))
 residual=sum((co[j]*powers[degree-j] for j in range(degree+1)),sp.zeros(5))
 cert.append(residual==sp.zeros(5))
# If scalar denominator only, certify observability-zero for25 successivesteps.
 if residual!=sp.zeros(5):
  for advance in range(2):residual=K(residual)
  for j in range(25):s.require(residual[0,4]==0,'CayleyHamilton scalarannihilator');residual=K(residual)
# finite exact independent r3blue target DP.
blue3=0
for t in letters:
 for u in fixed:
  target=s.complement(s.join((t,u,tau[t])),15)
  blue3+=s.hdata(target,15)[0][-1]
s.require(blue3==odd[1],'r3all72blueDP')
z=sp.symbols('z');blueGF=sp.cancel(rats[0].subs(w,z*z)+z*rats[1].subs(w,z*z))
result={'fixed_tiles':fixed,'N':{t:[list(map(int,row)) for row in matrix.tolist()] for t,matrix in N.items()},'tau':tau,'even_first8':even[:8],'odd_first8':odd[:8],'even_GF':str(rats[0]),'odd_GF':str(rats[1]),'blue_GF':str(blueGF),'matrix_annihilator':cert,'r3blue_independent_H':blue3}
(D/'blue-transfer-results.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k not in('N','tau')},indent=2))
