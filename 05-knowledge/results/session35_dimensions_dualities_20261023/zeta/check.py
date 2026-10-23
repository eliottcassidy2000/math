from itertools import product
import sympy as S
import mpmath as mp
from pathlib import Path
mp.mp.dps=45
checked=0
for n in range(1,6):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for bits in product([-1,1],repeat=len(pairs)):
        M=S.zeros(n)
        for (i,j),b in zip(pairs,bits): M[i,j]=b; M[j,i]=-b
        Q=S.eye(n)-M*M
        D=M.row_join(M+S.eye(n)).col_join((M-S.eye(n)).row_join(-M))
        assert S.eye(2*n)-D*D==S.diag(2*Q,2*Q)
        assert Q.trace()==n*n
        for k in range(4): assert S.trace(Q**k)>0
        checked+=1
print('FINITE-EXACT labeled tournament universe n=1..5:',checked)
# Paley correction, exact adjacency characteristic polynomials.
x=S.Symbol('x')
for p in [3,7,11,19]:
    residues={i*i%p for i in range(1,p)}
    A=S.Matrix(p,p,lambda i,j:int((j-i)%p in residues))
    expected=(x-S.Rational(p-1,2))*(x*x+x+S.Rational(p+1,4))**((p-1)//2)
    assert S.expand(A.charpoly(x).as_expr()-expected)==0
    print('Paley adjacency',p,':',S.factor(expected))
# Exact hostile: symmetric, real, positive on real axis, off-axis zeros.
s=S.Symbol('s'); a=S.Rational(1,4)
F=S.expand(((s-S.Rational(1,2)-a)**2+1)*((s-S.Rational(1,2)+a)**2+1))
assert S.expand(F.subs(s,1-s)-F)==0
assert S.expand(F.subs(s,S.Rational(1,4)+S.I))==0
print('Reflection+real-axis positivity hostile:',F,'zeros 1/4+-i, 3/4+-i')
# Numerical independent scalar heat-tower test, Q=[1].
for t in [mp.mpf('0.1'),mp.mpf('0.01'),mp.mpf('0.001')]:
    H=sum(mp.power(2,k)*mp.exp(-mp.power(2,k)*t) for k in range(18))
    P=sum(mp.power(2,k)*t*mp.exp(-mp.power(2,k)*t) for k in range(-160,18))
    correction=sum((-t)**j/(mp.factorial(j)*(mp.power(2,j+1)-1)) for j in range(30))
    err=abs(H-P/t+correction)
    assert err<mp.mpf('1e-40')
    print('NUMERICAL heat identity t=',t,'error=',mp.nstr(err,5))
# Exact n=1->2 lattice theta coefficients, Q=[1] to 2I_2.
N=20
a=[0]*(N+1)
for v in range(-4,5):
    if v*v<=N: a[v*v]+=1
b=[0]*(N+1)
for u,v in product(range(-4,5),repeat=2):
    if 2*(u*u+v*v)<=N: b[2*(u*u+v*v)]+=1
rhs=[sum(a[i]*a[j] for i in range(N+1) for j in range(N+1) if 2*(i+j)==k) for k in range(N+1)]
assert b==rhs
print('FINITE-EXACT theta doubling coefficients through exponent20:',b)
