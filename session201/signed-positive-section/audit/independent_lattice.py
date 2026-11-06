import sympy as s
import json
x=s.symbols('x');a,e,p,b,h,r=s.symbols('a e p b h r')
pol=x**3-x*x-x-1
R=s.Matrix([[0,0,1],[1,1,1],[0,1,0]])
Q=s.Matrix([[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]])
B=s.Matrix([[1,2,1],[0,1,0],[1,2,2],[0,1,1]])
C=s.Matrix([[1,0,0],[0,1,-1],[0,0,1]])
M=s.Matrix([[0,0,1],[1,0,1],[0,1,1]])
J=s.Matrix([[1,-1,1,0],[1,-1,0,0],[1,0,1,0],[1,0,0,0]])
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def ell(v):
 z=s.Poly(s.rem(v,pol,x),x)
 return sum(c*z.nth(i) for i,c in enumerate([1,2,4]))
need(s.Poly(pol,x).is_irreducible,'irreducible cubic')
need(R.charpoly(x).as_expr()==pol and R.det()==C.det()==1,'unimodular conjugacy')
need(C*R==M*C and J*B==B*R,'unit action')
need(s.rem(x*(x*x-x-1),pol,x)==1,'unit inverse')
need(s.rem(x**6*(13-7*x),pol,x)==3+x,'boundary')
v=B*s.Matrix([a,e,p]);w=B*s.Matrix([b,h,r])
need(s.expand((v.T*Q*w)[0]-ell((a+(e-p)*x+p*x*x)*(b+(h-r)*x+r*x*x)))==0,'symbolic pairing')
G=s.simplify(C.inv().T*B.T*Q*B*C.inv())
need(G==s.Matrix([[1,2,4],[2,4,7],[4,7,13]]) and G.det()==-1,'perfect pairing')
q=s.Matrix([5,1,5,1]); inv=B*R.inv()*s.Matrix([3,1,0])
need(inv==s.Matrix([1,0,4,3]) and J*inv==q,'unique lattice inverse')
need(J*s.Matrix([1,0,4,1])==q,'actual off-face inverse')
need(s.Matrix([[ -1,1,0,1]])*s.Matrix([1,0,4,3])==s.Matrix([2]),'negative a equivalent')
print(json.dumps({'status':'PASS','method':'SymPy exact quotient-ring remainders, symbolic polynomial identity, and matrix conjugacy','sympy_version':s.__version__,'gram':list(map(list,G.tolist())),'lattice_inverse':list(inv)},default=str,indent=2))
