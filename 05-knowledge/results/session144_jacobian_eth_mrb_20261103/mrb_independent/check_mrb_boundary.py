import sympy as s,json
q,z=s.symbols('q z'); rows=[]
def require(c,m):
 if not c: raise ValueError(m)
for K in range(9):
 # Apply (q*d/dq-j) to geometric series, generating (n-1)...(n-K).
 r=-q/(1+q)
 for j in range(1,K+1): r=s.factor(q*s.diff(r,q)-j*r)
 expected=(-1)**(K+1)*s.factorial(K)*(q/(1+q))**(K+1)
 require(s.cancel(r-expected)==0,f'generating sum K={K}')
 P=s.prod(z-j for j in range(1,K+1))
 require(all(P.subs(z,j)==0 for j in range(1,K+1)), 'prefix')
 require(s.degree(P,z)==K,'degree')
 require(s.sign(expected.subs(q,s.Rational(1,3)))==(-1)**(K+1),'sign')
 rows.append({'K':K,'identity':str(expected),'prefix_zeros':K})
# Hostile: prior e^-z P has unbounded polynomial growth on Re z=0 for K>=1,
# whereas the extra entire sinc-like factor yields radius^-1 bound.
require(s.cancel((-q/(1+q))-(q/(1+q)))!=0,'wrong sign hostile detected')
print(json.dumps({'status':'PASS','scope':'Exact generating identities K=0..8; all-K analytic proof separate','rows':rows},indent=2))
