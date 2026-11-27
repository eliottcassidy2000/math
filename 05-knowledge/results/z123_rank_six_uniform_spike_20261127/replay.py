import sympy as S
from itertools import product
from functools import lru_cache
q=S.symbols('q')
def gauss(n,k):
 return S.Poly(S.cancel(S.prod((q**(n-i)-1)/(q**(i+1)-1) for i in range(k))),q).as_expr()
def choose(x,k):return S.prod(x-j for j in range(k))/S.factorial(k)
def G(q):
 m=q**3;b=q**9;N=[gauss(6,k) for k in (1,2,3)];a,l,h=N;t=q*q+q+1
 types=[(a*choose(m,2)+choose(m,3),3,3,6),(b,t+1,2*t,(q+3)*t),
        (m*l-b*t,q+2,q+2,3*(q+1)),(h-b,t,t,(q+1)*t)]
 return {():1,(1,):a+m,(2,):l+a*m+choose(m,2),(3,):sum(r[0] for r in types),
     (1,2):(q+1)*l+2*(a*m+choose(m,2)),(1,3):sum(r[0]*r[1] for r in types),
     (2,3):sum(r[0]*r[2] for r in types),(1,2,3):sum(r[0]*r[3] for r in types)}
def uniform(r,n):
 from itertools import combinations
 return {I:S.prod(choose(n-a,b-a) for a,b in zip((0,)+I,I)) for j in range(r) for I in combinations(range(1,r),j)}
def make(r,s,L,R):
 @lru_cache(None)
 def f(*H):
  a=0
  for I in product(range(r+1),repeat=len(H)):
   J=tuple(h-i for h,i in zip(H,I))
   if tuple(sorted(I))!=I or tuple(sorted(J))!=J or not all(0<=j<=s for j in J): continue
   a+=L[tuple(sorted(set(i for i in I if 0<i<r)))]*R[tuple(sorted(set(j for j in J if 0<j<s)))]
  return S.expand(a)
 return f
def z(f,R):
 a=f(R-1)
 b=f(R-2)+f(R-3,R-1)-f(R-3,R-2)
 c=f(R-3)+f(R-4,R-1)-f(R-4,R-3)+f(R-5,R-2)-f(R-5,R-3)+f(R-5,R-3,R-1)-f(R-5,R-3,R-2)-f(R-5,R-4,R-1)+f(R-5,R-4,R-3)
 return [S.expand(x) for x in (a,b,c)]
n,t=S.symbols('n t')
Gq=G(q)
Z6=z(make(4,3,Gq,uniform(3,n)),6)
A,B,C=Z6
D=S.Poly(S.expand(B*B-A*C),n,q)
# Independent rank-6 extraction formulas with a Boolean three-element side.
T=S.symbols('w1 w2 w3 f12 f13 f23 f123')
w1,w2,w3,f12,f13,f23,f123=T
F={():1};F.update(dict(zip([(1,),(2,),(3,),(1,2),(1,3),(2,3),(1,2,3)],T)))
assert z(make(4,3,F,uniform(3,3)),6)==[
 w2+3*w3+3,
 3*f13+3*f23-3*w1+4*w2+5*w3+1,
 -f12+f123+7*f13+5*f23-7*w1+7*w2+w3+1]
assert all(min(S.Poly(D.as_expr().subs(q,i),n).shift(3).all_coeffs())>0 for i in range(2,180))
dh=S.Poly(D.as_expr().subs(q,181),n)
assert all(dh.eval(i)>0 for i in range(3,76))
ds=D.shift_list([76,181])
assert ds.eval({n:0,q:0})<0 and all(c<0 for c in ds.coeffs())
print('uniform mixed family: index 2 passes q <= 179, fails on q >= 181 and n >= 76')
def analyze(n0,q0):
 a,b,c=[S.expand(f.subs(n,n0)) for f in Z6]
 gap=S.Poly(b*b-a*c,q)
 assert max(gap.shift(q0).all_coeffs())<0
 f=t**3+a*t**2+(b-3)*t+c-2*a
 delta=S.Poly(S.discriminant(f,t),q)
 edge=S.Poly(f.subs(t,-2),q)
 assert max(delta.shift(q0).all_coeffs())<0 and min(edge.shift(q0).all_coeffs())>0
 print('rank6 witness',q0,n0,'elements',Gq[(1,)].subs(q,q0)+n0)
 print('z=',[x.subs(q,q0) for x in [a,b,c]],'gap=',gap.eval(q0))
 print('degrees/leading',[(S.degree(x,q),S.Poly(x,q).LC()) for x in [a,b,c]],'gap',gap.degree(),gap.LC())
 print('cubic discriminant shift',delta.degree(),len(delta.shift(q0).all_coeffs()),'all negative; F(-2) all positive')
 return gap
analyze(76,181)
gap=analyze(3,239)
assert all(gap.eval(i)>0 for i in range(2,239))
a,b,c=z(make(4,4,Gq,Gq),7);D7=S.Poly(b*b-a*c,q)
assert all(D7.eval(i)>0 for i in range(2,26))
assert max(D7.shift(26).all_coeffs())<0
print('linked rank7: formal integer crossover 26; smallest field size on failing side is 27')
