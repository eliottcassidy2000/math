from itertools import product
from collections import Counter
import math
import numpy as np
import sympy as S
import mpmath as mp
mp.mp.dps=35
# Independent Paley identity via adjacency algebra, extending original test set.
for p in [3,7,11,19,23,31,43]:
    residues={j*j%p for j in range(1,p)}
    A=S.Matrix(p,p,lambda i,j:int((j-i)%p in residues))
    J=S.ones(p); I=S.eye(p); M=A-A.T
    assert A+A.T==J-I
    assert M*M==J-p*I
    assert 4*A*A+4*A+(p+1)*I==(p+1)*J
    # on 1-perp: 4 A^2+4 A+(p+1)I=0, hence modulus sqrt(p+1)/2.
print('EXACT Paley adjacency/skew polynomial identities: p=3,7,11,19,23,31,43 PASS')
# General non-scalar Q, exact lattice doubling coefficients with rigorous box bound Q>=I.
M=S.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
Q=S.eye(3)-M*M
assert Q.eigenvals()=={S.Integer(1):1,S.Integer(4):2}
N=16
counts=Counter()
R=math.isqrt(N)
for v in product(range(-R,R+1),repeat=3):
    w=S.Matrix(v); k=int((w.T*Q*w)[0])
    if k<=N: counts[k]+=1
# direct six-dimensional enumeration; 2Q>=2I permits R=floor(sqrt(N/2)).
D=M.row_join(M+S.eye(3)).col_join((M-S.eye(3)).row_join(-M))
QD=S.eye(6)-D*D
assert QD==S.diag(2*Q,2*Q)
qn=np.array(QD).astype(int); actual=Counter(); R=math.isqrt(N//2)
for v in product(range(-R,R+1),repeat=6):
    w=np.array(v); k=int(w@qn@w)
    if k<=N: actual[k]+=1
expected=Counter()
for i,ci in counts.items():
    for j,cj in counts.items():
        if 2*(i+j)<=N: expected[2*(i+j)]+=ci*cj
assert actual==expected
print('EXACT non-scalar 3-to-6 theta doubling through exponent16:',sorted(actual.items()))
# Non-scalar tower: test formula and an explicit error bound, not only cancellation.
qs=[mp.mpf(1),mp.mpf(4),mp.mpf(4)]
for t in map(mp.mpf,['0.0003','0.03','0.7']):
    h=lambda x:sum(mp.exp(-x*q) for q in qs)
    H=sum(2**k*h(2**k*t) for k in range(22))
    P=lambda u:sum(mp.power(2,k)*u*h(mp.power(2,k)*u) for k in range(-140,22))
    correction=sum((-t)**m*sum(q**m for q in qs)/(mp.factorial(m)*(2**(m+1)-1)) for m in range(70))
    assert abs(H-P(t)/t+correction)<mp.mpf('1e-28')
    assert abs(P(2*t)-P(t))<mp.mpf('1e-28')
    remainder=H-P(t)/t+3
    assert 0<=remainder<=t*sum(qs)/3
    assert abs(remainder-t*sum(qs)/3)<=t*t*sum(q*q for q in qs)/14
print('NUMERICAL non-scalar tower identity, log-periodicity, first/second remainder bounds PASS')
# Matrix theta inversion on a genuinely non-diagonal Q; numerical box truncation only.
q=np.array(Q).astype(float); qi=np.linalg.inv(q)
vec=np.array(list(product(range(-13,14),repeat=3)))
def theta(a,t): return np.exp(-math.pi*t*np.einsum('bi,ij,bj->b',vec,a,vec)).sum()
for t in [.15,.4,1.3]:
    lhs=theta(q,t); rhs=np.linalg.det(q)**(-.5)*t**(-1.5)*theta(qi,1/t)
    assert abs(lhs-rhs)<1e-11
print('NUMERICAL non-diagonal Poisson theta inversion PASS')
# Epstein split-integral completion: diagonal positive forms in dimensions 1..5.
def theta_diag(a,t):
    return mp.fprod(mp.jtheta(3,0,mp.exp(-mp.pi*t*x)) for x in a)
def lam(a,s):
    n=len(a); d=mp.fprod(a)**(-mp.mpf('.5')); inv=[1/x for x in a]
    term1=mp.quad(lambda t:(theta_diag(a,t)-1)*t**(s/2-1),[1,2,5,mp.inf])
    term2=mp.quad(lambda t:(theta_diag(inv,t)-1)*t**((n-s)/2-1),[1,2,5,mp.inf])
    return -2/s+2*d/(s-n)+term1+d*term2
for n in range(1,6):
    a=[mp.mpf(1)+mp.mpf(j)/5 for j in range(n)]
    s=mp.mpc('.37','.71'); det=mp.fprod(a); dual=[1/x for x in a]
    err=abs(lam(a,s)-det**(-mp.mpf('.5'))*lam(dual,n-s))
    assert err<mp.mpf('1e-29')
    # Strictly positive completed value implies simple zero after reciprocal Gamma.
    for m in [1,2]: assert lam(a,mp.mpf(-2*m))>0
    b=[x/det**(mp.mpf(1)/n) for x in a]; bi=[1/x for x in b]
    assert abs(lam(b,s)+lam(bi,s)-lam(b,n-s)-lam(bi,n-s))<mp.mpf('1e-28')
    print('NUMERICAL Epstein reflection/negative-even completed positivity/determinant normalization dimension',n,'PASS')
print('ALL INDEPENDENT CHECKS PASS; numerical cases are not a substitute for analytic proofs.')
