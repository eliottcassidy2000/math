import json,sympy as s
from pathlib import Path
D=Path(__file__).parent
raw=json.loads((D/'compressed-matrices.json').read_text());M={int(t):s.Matrix(m) for t,m in raw.items()}
prev=json.loads((D/'blue-transfer-results.json').read_text());tau={int(t):u for t,u in prev['tau'].items()};G=s.diag(1,-1,1)
def require(ok,msg):
 if not ok:raise RuntimeError(msg)
for t in M:require(M[tau[t]]==G*M[t].T*G,'perletterLorentzianadjoint')
h,ss,r,w=s.symbols('h ss r w');v=s.Matrix([h,ss,r]);mons=[h*h,ss*ss,r*r,h*ss,h*r,ss*r]
def coeff(p):
 p=s.Poly(s.expand(p),h,ss,r)
 return [p.coeff_monomial(m) for m in mons]
L=s.zeros(6)
for A in M.values():
 x,y,z=A*v;out=[x*x,y*y,z*z,x*y,x*z,y*z]
 L+=s.Matrix([coeff(p) for p in out])
e0=s.Matrix([1,0,0,0,0,0]);fixed=sum((M[t] for t in M if tau[t]==t),s.zeros(3));cost_even=s.Matrix([[1,-1,1,0,0,0]]);cost_odd=s.Matrix([coeff((v.T*G*fixed*v)[0])]);res=(s.eye(6)-w*L).inv()*e0
even=s.cancel((cost_even*res)[0]);odd=s.cancel((cost_odd*res)[0])
require(s.cancel(even-1-s.sympify(prev['even_GF']))==0,'even rational exact')
require(s.cancel(odd-s.sympify(prev['odd_GF']))==0,'odd rational exact')
result={'moment_basis':['H^2','S^2','R^2','HS','HR','SR'],'moment_matrix':L.tolist(),'even_cost':list(cost_even),'odd_cost':list(cost_odd),'denominator':str(s.factor((s.eye(6)-w*L).det())),'even_GF_including_empty':str(even),'odd_GF':str(odd),'adjoint_verified_letters':18}
(D/'lorentzian-results.json').write_text(json.dumps(result,indent=2,default=int)+'\n');print(json.dumps(result,indent=2,default=int))
