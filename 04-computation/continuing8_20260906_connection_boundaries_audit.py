"""Independent actual-monomial, finite-field-rank and differential referee."""
from fractions import Fraction as F
from math import comb
from pathlib import Path
from hashlib import sha256
import sys
sys.stdout.reconfigure(newline="\n")
GATES=0
def need(ok,why):
    global GATES
    GATES+=1
    if not ok: raise ArithmeticError(why)
here=Path(__file__).resolve().parent
primary=here/"continuing8_20260906_connection_boundaries.py"
need(sha256(primary.read_bytes()).hexdigest()=="bba0c468e6a42142ec623b1a5bd2deb1a3e08d8398629ea13e5ec4f3dd9cae43","frozen producer source")

def add(A,B):
    C=dict(A)
    for k,v in B.items():
        C[k]=C.get(k,0)+v
        if not C[k]: del C[k]
    return C
def times(A,B):
    C={}
    for a,u in A.items():
        for b,v in B.items():
            k=tuple(x+y for x,y in zip(a,b)) if isinstance(a,tuple) else a+b
            C[k]=C.get(k,0)+u*v
    return {k:v for k,v in C.items() if v}
def scaled(A,c): return {k:c*v for k,v in A.items() if c*v}
def power(A,n,one):
    B=one
    while n:
        if n&1: B=times(B,A)
        A=times(A,A); n//=2
    return B
def derivative(A): return {k-1:k*v for k,v in A.items() if k}
def op(A,index,parameter=6):
    return add(times({0:parameter,2:1},derivative(A)),scaled(times({1:1},A),-2*index))

# Direct sparse coefficient arithmetic, without the producer's symbolic engine.
for m in (-3,1,2,5):
    for a in (-1,1):
        phi={1:a}
        for n in range(10):
            left=op(power(phi,n,{0:1}),m)
            right=scaled(add(scaled(power(phi,max(0,n-1),{0:1}),6*n),
                              scaled(power(phi,n+1,{0:1}),n-2*m)),a)
            need(left==right,"both permitted charts on complete monomial test bank")

charts=[{2:1},{1:2,2:-1},{0:1,2:1},{3:1},{0:1,1:2}]
for phi in charts:
    B=add(times(times({0:6,2:1},phi),derivative(phi)),
          scaled(times({1:1},add(power(phi,2,{0:1}),{0:6})),-1))
    need(bool(B),"nonlinear/shifted chart coefficient obstruction")
    for n in range(1,8):
        # Clear the forced rational scalar x/phi at nonzero index.
        left=times(phi,op(power(phi,n,{0:1}),2))
        response=add(scaled(power(phi,n-1,{0:1}),6*n),scaled(power(phi,n+1,{0:1}),n-4))
        residual=add(left,scaled(times({1:1},response),-1))
        need(residual==times(B,scaled(power(phi,n-1,{0:1}),n)),"full forced-intertwiner residual factor")
        # At m=0, clear the denominator of the freely chosen rational scalar.
        den=add(power(phi,2,{0:1}),{0:6})
        hnum=times({0:6,2:1},derivative(phi))
        target=add(scaled(power(phi,n-1,{0:1}),6*n),scaled(power(phi,n+1,{0:1}),n))
        need(times(den,op(power(phi,n,{0:1}),0))==times(hnum,target),"m0 rational-scalar exception")
for a in (F(2),F(3,2),F(-3)):
    phi={1:a}
    for m in (1,2,5):
        for n in range(8):
            left=op(power(phi,n,{0:1}),m)
            right=scaled(add(scaled(power(phi,max(0,n-1),{0:1}),6*a*a*n),
                              scaled(power(phi,n+1,{0:1}),n-2*m)),1/a)
            need(left==right,"quadratic-parameter covariance with actual rational scale")
need(op({0:F(-1,2)},2)=={1:F(2)},"actual compatible Student target")
for c0 in range(4):
    for c2 in range(4):
        need((12*c2-4*c0-2)%4!=0,"mod4 rejection holds even in unrestricted polynomial operator")
need(op({},2)=={},"mod2 target2x admits the zero source")
print("STUDENT exact chart identities, rational scales, m0 exception and unrestricted mod4 obstruction PASS")

# Bivariate multiplication starts from the actual x,t source expressions.
pi={(0,1):1,(2,2):1}
y={(1,2):1,(3,3):1}
one={(0,0):1}
core=add(power(pi,3,one),scaled(power(y,2,one),-1))
need(core=={(0,3):1,(2,4):2,(4,5):1},"actual integer source identity by sparse multiplication")

def rank_mod(rows,p):
    A=[[x%p for x in row] for row in rows]
    if not A:return 0
    pivot=0
    for col in range(len(A[0])):
        r=next((i for i in range(pivot,len(A)) if A[i][col]),None)
        if r is None: continue
        A[pivot],A[r]=A[r],A[pivot]
        inv=pow(A[pivot][col],-1,p)
        A[pivot]=[(x*inv)%p for x in A[pivot]]
        for i in range(pivot+1,len(A)):
            c=A[i][col]
            if c:A[i]=[(u-c*v)%p for u,v in zip(A[i],A[pivot])]
        pivot+=1
        if pivot==len(A):break
    return pivot

def determinant(A):
    A=[row[:] for row in A]; n=len(A); last=1; parity=1
    for k in range(n-1):
        r=next((i for i in range(k,n) if A[i][k]),None)
        if r is None:return 0
        if r!=k:A[k],A[r]=A[r],A[k]; parity*=-1
        pivot=A[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=A[i][j]*pivot-A[i][k]*A[k][j]
                need(numerator%last==0,"fraction-free determinant exact division")
                A[i][j]=numerator//last
            A[i][k]=0
        last=pivot
    return parity*A[-1][-1]

cases=[(p,a) for p in (2,3,5,7,11,13) for a in (1,2,3) if p**a<=125]+[(2,4),(2,5),(3,4),(17,1)]
need(len(cases)==18,"declared producer and fresh prime-power cases")
for p,a in cases:
    P=p**a; L=P//2; eps=P%2; T=3*L+eps
    solutions=[((2*T-3*b)//2,b) for b in range(2*T//3+1) if (2*T-3*b)%2==0]
    need(solutions==[(T-3*j,2*j) for j in range(L+1)],"complete native monomial intercept universe")
    # Native unsigned diagonal coefficient matrix, without removing a factor.
    native=[[comb(T-j,k-j) if j<=k<=T else 0 for j in range(L+1)] for k in range(P+1)]
    for order,expected in ((P-1,True),(P,False)):
        matrix=native[:order+1]
        augmented=[row+[int(k==0)] for k,row in enumerate(matrix)]
        need((rank_mod(matrix,p)==rank_mod(augmented,p))==expected,"literal full-source projected membership and first failed row")
    # All columns retain an integral unit initial block, independent of p.
    need(all(native[j][j]==1 and all(native[k][j]==0 for k in range(j)) for j in range(L+1)),
         "integral saturation at every finite source projection")
    need(T+L+1==4*T//3+1 and T+P-(T+L+1)==(P+1)//2-1,"both recognition cutoffs and exact delay")
    if P<=32:
        literal=times(power(pi,eps,one),power(core,L,one))
        need(literal=={(2*k,T+k):comb(P,k) for k in range(P+1)},"complete actual-source lift before characteristic reduction")
        reduction={k:v%p for k,v in literal.items() if v%p}
        need(reduction=={(0,T):1,(2*P,T+P):1},"exact two-term Frobenius lift")
        for j,(c,e) in enumerate(solutions):
            column=times(power(pi,c,one),power(y,e,one))
            expected={(2*k,T+k):comb(T-j,k-j) for k in range(j,T+1)}
            need(column==expected,"each actual monomial column preserves its entire diagonal")
        matrix=[[comb(T-j,k-j) if j<=k<=T else 0 for j in range(L+1)]+[int(k==0)] for k in range(L+2)]
        need(abs(determinant(matrix))==comb(P+L,L+1),"independent characteristic-zero first-obstruction determinant")
    print("FROBENIUS",p,a,"P",P,"T",T,"char0",T+L+1,"modp",T+P)
print("SCOPE characteristic p with P=p^a; independent x,t; full native source; inclusive row projection")
print("PASS",GATES,"always-active independent exact gates; raw LF")
