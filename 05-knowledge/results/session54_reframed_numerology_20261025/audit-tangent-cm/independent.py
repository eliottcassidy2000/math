import sympy as S, json, math
from pathlib import Path
x,z,X=S.symbols('x z X')
def req(c,m):
 if not c: raise RuntimeError(m)
data=json.loads((Path(__file__).resolve().parent/'inputs'/'producer-results.json').read_text())
rows=data['candidate_conductors']; out=[]
for row in rows:
 n=row['n']; ph=S.Poly(S.cyclotomic_poly(n,z),z); d=ph.degree()//2
 # Recover psi by triangular coefficient subtraction, without resultants.
 rem=ph.as_expr(); psi=0
 for k in range(d,-1,-1):
  c=S.expand(rem).coeff(z,d+k)
  psi+=c*x**k
  rem=S.expand(rem-c*S.expand(z**d*(z+1/z)**k))
 req(rem==0,'reconstruction')
 P=S.Poly(psi,x); den=S.Poly((2-x)**2*(2+x),x)
 jf=S.rem(S.Poly(64*(10+3*x)**3,x)*S.invert(den,P),P)
 cols=[]
 for k in range(d):
  p=S.rem(jf*S.Poly(x**k,x),P)
  cols.append([p.nth(i) for i in range(d)])
 M=S.Matrix.hstack(*[S.Matrix(c) for c in cols]); cp=M.charpoly(X).as_expr()
 expected=S.sympify(row['polynomial'],locals={'J':X})
 req(S.expand(cp-expected)==0,'multiplication charpoly n='+str(n))
 req(S.gcd(cp,S.diff(cp,X))==1,'collision')
 out.append({'n':n,'degree':d,'multiplication_charpoly_matches':True})
# Independent bounded brute screen up to 5000, not used to prove completeness.
ns=[]
for n in range(3,5001):
 v={min(a*a%n,(-a*a)%n) for a in range(n) if math.gcd(a,n)==1}
 if len(v)<=3:ns.append(n)
req(ns==[r['n'] for r in rows],'brute group screen to5000')
# Direct symbolic Legendre substitution, with c=cos(2pi/m), x=4c²-2.
c=S.symbols('c'); lam=(1-c)/(1+c)
j=256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2)
req(S.cancel(j-64*(10+3*(4*c*c-2))**3/((2-(4*c*c-2))**2*(2+(4*c*c-2))))==0,'j map')
print(json.dumps({'status':'PASS','independent_multiplication_matrices':out,'group_brute_limit':5000,'symbolic_legendre':True},indent=2))
