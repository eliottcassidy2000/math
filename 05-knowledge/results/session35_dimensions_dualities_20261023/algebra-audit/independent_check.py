"""Independent Session35 checks: direct rational Gram elimination and linear solves.
No hook-content enumeration, adjugate, inverse, or producer imports are used.
"""
from fractions import Fraction as F
from itertools import permutations
import json
import sympy as s

def require(x, label):
    if not x: raise RuntimeError(label)
def cyc(p):
    unseen=set(range(len(p))); result=0
    while unseen:
        a=next(iter(unseen)); result+=1
        while a in unseen:
            unseen.remove(a); a=p[a]
    return result

def direct_gram(n,t):
    ps=list(permutations(range(n))); invs=[]
    for p in ps:
        inv=[0]*n
        for i,j in enumerate(p): inv[j]=i
        invs.append(inv)
    return [[t**cyc(tuple(inv[q[i]] for i in range(n))) for q in ps] for inv in invs]
def psd(A):
    # Symmetric rational Schur elimination; zero diagonal forces zero row in PSD.
    A=[row[:] for row in A]
    while A:
        if any(A[i][i]<0 for i in range(len(A))): return False
        for i in range(len(A)):
            if A[i][i]==0 and any(A[i][j] for j in range(len(A))): return False
        piv=next((i for i in range(len(A)) if A[i][i]>0),None)
        if piv is None: return True
        ids=[i for i in range(len(A)) if i!=piv]; d=A[piv][piv]
        A=[[A[i][j]-A[i][piv]*A[piv][j]/d for j in ids] for i in ids]
    return True

gram_checks=0
for n in range(1,5):
    for j in range((-4 if n==1 else 0),4*(n+1)+1):
        t=F(j,4)
        expected=t>=0 and (t.denominator==1 or t>=n-1)
        require(psd(direct_gram(n,t))==expected,('direct Gram PSD',n,t))
        gram_checks+=1
# Direct sign quadratic, not interpolated binomial: row-sum eigenvalue at order5.
n=5;t=F(7,2)
sign_eigen=sum((-1)**(n-cyc(p))*t**cyc(p) for p in permutations(range(n)))
require(sign_eigen==F(-105,32), 'order5 alternating eigenvalue')

# Build borders by indexed entries, independently of producer block constructors.
def extend(M):
    n=M.rows; N=s.zeros(n+2)
    for i in range(n+2):
        for j in range(i+1,n+2):
            a=M[i-1,j-1] if 0<i<j<n+1 else 1
            N[i,j]=a;N[j,i]=-a
    return N

def double(M):
    n=M.rows; N=s.zeros(2*n)
    for i in range(2*n):
        for j in range(i+1,2*n):
            if j<n: a=M[i,j]
            elif i>=n: a=-M[i-n,j-n]
            else: a=M[i,j-n]+int(i==j-n)
            N[i,j]=a;N[j,i]=-a
    return N

def f(M,x):
    b=s.ones(M.rows,1)
    if isinstance(x,s.Basic) and x.free_symbols:
        unknowns=s.symbols('v:'+str(M.rows))
        solution=s.solve_linear_system((x*s.eye(M.rows)-M).row_join(b),*unknowns)
        return s.factor(sum(solution[v] for v in unknowns))
    return (b.T*(x*s.eye(M.rows)-M).LUsolve(b))[0]
recursion_checks=0
for n in range(1,7):
    for seed in range(4):
        M=s.zeros(n)
        for i in range(n):
            for j in range(i+1,n):
                M[i,j]=((i+2)*(j+3)+seed*7)%9-4;M[j,i]=-M[i,j]
        N=extend(M); W=double(M)
        for x,r in [(1,1),(7,5),(41,29)]:
            fm=f(M,x); fn=f(N,x)
            require(fn==((x*x+1)*fm+2*x)/(x*x+1+2*x*fm),('E marked',n,seed,x))
            require(f(W,x)==s.Rational(x,r)*f(M,r),('D marked',n,seed,x))
            require((x*s.eye(n+2)-N).det()==(x*s.eye(n)-M).det()*(x*x+1+2*x*fm),('E determinant',n,seed,x))
            require((x*s.eye(2*n)-W).det()==((x*x+1)*s.eye(n)-2*M*M).det(),('D determinant',n,seed,x))
            recursion_checks+=4
print(json.dumps({'direct_gram_psd_checks':gram_checks,'order5_sign_eigenvalue':str(sign_eigen),'weighted_skew_recursion_checks':recursion_checks,'dimensions':'1..6','test_points':'(z,r)=(1,1),(7,5),(41,29); 2r^2=z^2+1'},indent=2))
# Pole/cancellation controls with known kernel and marked-orthogonal spectral blocks.
z=s.symbols('z')
C=s.Matrix([[0,1,-1],[-1,0,1],[1,-1,0]])
J=s.Matrix([[0,1],[-1,0]])
require(s.cancel(f(C,z)-3/z)==0,'cycle mark misses nonzero spectrum')
require(s.cancel(f(extend(C),z)-(5*z*z+3)/(z*(z*z+7)))==0,'E canceled old spectral poles')
require(s.cancel(f(double(C),z)-6*z/(z*z+1))==0,'D genuine w=0 poles')
require(s.cancel(f(double(J),z)-4*z/(z*z+3))==0,'D removable w=0 singularity')
require(s.cancel(f(double(J),z)).subs(z,s.I)==2*s.I,'finite branch point value')
print('Five symbolic pole/cancellation controls pass.')
