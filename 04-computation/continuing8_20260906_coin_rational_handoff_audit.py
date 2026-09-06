"""Independent algebraic and literal-spine controls; no producer imports."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from math import comb
from pathlib import Path
import json
import sys
import sympy as s
sys.stdout.reconfigure(newline='\n')
G=0
HERE=Path(__file__).resolve().parent
STEM='continuing8_20260906_coin_rational_handoff'
def gate(ok,label):
    global G
    G+=1
    if not ok: raise RuntimeError(label)
def find(name):
    for p in (HERE/name,HERE.parent/'05-knowledge/results'/name,
              Path('C:/w/continuing8_20260906_coin')/name):
        if p.is_file(): return p
    raise FileNotFoundError(name)
for ext,pin in (('.py','359eba70c2e89226e44486441919d6a95a911188054de6fd19f9e09de7ead804'),
                ('.out','94e04c1f7b92b12f70952e5542370fec4ccb427f3e7890ddd550ea3d2e84c86e')):
    raw=find(STEM+ext).read_bytes()
    gate(sha256(raw).hexdigest()==pin,'frozen producer '+ext)
    if ext=='.out': gate(b'\r' not in raw,'actual producer LF')

x,a,b=s.symbols('x a b')
for d in range(7):
    M=s.Matrix(d+1,d+1,lambda row,col:comb(col,row)*a**row*b**(col-row) if col>=row else 0)
    gate(s.expand(M.det()-a**(d*(d+1)//2))==0,'affine coefficient change determinant, arbitrary nonzero unit slope')
for prime in (2,3,5,7):
    for degree in range(1,5):
        gate(all(pow(aa,degree*(degree+1)//2,prime) for aa in range(1,prime)),'all invertible affine slopes preserve nonzero polynomial')
Fd=lambda d,z:z/(1-d*z)
gate(s.cancel(Fd(3,4*x+1)+Fd(6,-x)+s.Rational(1,2))==0,'nonunit slope genuinely allows denominator prime')

# Separate root-proposed extension: substitution need only remain nonconstant.
# Check all degree-at-most-two pairs over three small fields, using coefficient
# convolution rather than symbolic substitution or the producer's affine map.
composition_controls=0
def trim_mod(v,p):
    v=[c%p for c in v]
    while len(v)>1 and not v[-1]: v.pop()
    return v
def mul_mod(v,w,p):
    z=[0]*(len(v)+len(w)-1)
    for i,aa in enumerate(v):
        for j,bb in enumerate(w): z[i+j]=(z[i+j]+aa*bb)%p
    return trim_mod(z,p)
for prime in (2,3,5):
    polynomials=[trim_mod(v,prime) for v in product(range(prime),repeat=3)]
    for outer in polynomials:
        if outer==[0]: continue
        for inner in polynomials:
            if len(inner)==1: continue
            value=[0]
            for coefficient in reversed(outer):
                value=mul_mod(value,inner,prime)
                value[0]=(value[0]+coefficient)%prime
                value=trim_mod(value,prime)
            gate(value!=[0] and len(value)-1==(len(outer)-1)*(len(inner)-1),
                 'nonconstant polynomial composition is injective over a field')
            composition_controls+=1
    for inner in ([0,0,1],[0,2,-1]):
        gate(len(trim_mod(inner,prime))>1,'x^2 and 2x-x^2 retain positive degree')
gate(composition_controls==15546,'complete degree-two composition control universe')

def coeffs(box):
    d=len(box)-1
    return [sum(box[k]*(-1)**(j-d+k)*comb(k,j-d+k)
                for k in range(d+1) if 0<=j-d+k<=k) for j in range(d+1)]
def ev(v,z):
    out=0
    for c in reversed(v): out=out*z+c
    return out
boxes=0
for d in range(5):
    words=list(product((0,1),repeat=d))
    layers=[[w for w in words if sum(w)==k] for k in range(d+1)]
    for box in product(*(range(comb(d,k)+1) for k in range(d+1))):
        chosen=[w for k,layer in enumerate(layers) for w in layer[:box[k]]]
        c=coeffs(box)
        gate(all(isinstance(v,int) for v in c),'integer native Taylor coefficients')
        for z in (F(1,3),F(2,5)):
            literal=sum(z**(d-sum(w))*(1-z)**sum(w) for w in chosen)
            gate(ev(c,z)==literal,'literal suffix leaves equal reconstructed power polynomial')
        boxes+=1
gate(boxes==782,'entire bounded Bernstein box universe')

# A stationary actual polynomial realization, summed by its exact matrix resolvent.
M=s.Matrix([[2*x-1,1-x],[0,1]])
state=s.Matrix([1,1]); R=(s.eye(2)-x*M).inv()
FF=s.cancel(x*(1-x)*(s.Matrix([[1,0]])*R*state)[0])
GG=s.cancel(x*(1-x)*(s.Matrix([[-1,1]])*R*state)[0])
gate(s.cancel(FF-x*(1+x)/(1+2*x))==0 and s.cancel(GG-x*x/(1+2*x))==0,'exact stationary orientation matrix sums')
gate(s.cancel(FF+GG.subs(x,1-x)-s.Rational(1,2))!=0,'legal stationary source is not exactly fair')
gate(FF.subs(x,s.Rational(1,3))+GG.subs(x,s.Rational(2,3))==s.Rational(16,35),'stationary exact bias hostile')
Fvn=x*(1-x)+x*x/2; Gvn=x*x/2
gate(s.expand(Fvn+Gvn.subs(x,1-x))==s.Rational(1,2),'literal first-pair von Neumann orientation decomposition')
gate(s.Poly(Fvn,x).nth(2)==-s.Rational(1,2),'integral-germ hypothesis fails for von Neumann')
for k in range(20):
    word='0011'+'00'*k+'01'
    mismatch=next(i for i,c in enumerate(word) if c!=word[0])
    stop=next(i+2 for i in range(0,len(word),2) if word[i]!=word[i+1])
    gate(mismatch==2 and stop==6+2*k,'same critical value and arbitrarily late von Neumann halt')

K2=(1-2*x)**4-x*x-(1-x)**2
gate(s.cancel(K2/(x*(1-x))+6-16*x*(1-x))==0,'literal first signed handoff load')
gate((K2/(x*(1-x))).subs(x,s.Rational(1,4))==-3,'aggregate cannot be a legal signed row')
for N in range(2,10):
    poly=s.Poly(s.expand((1-2*x)**(2*N)-x**N-(1-x)**N),x)
    par=s.Poly(1-x**N-(1-x)**N,x)
    gate(all(int(v)%2==0 for v in (poly-par).all_coeffs()),'aggregate parity despite illegal first row')
    gate(poly.eval(0)==poly.eval(1)==0,'aggregate endpoint condition')
    for q in (F(1,7),F(1,4),F(2,5)):
        gate(abs(F(poly.eval(s.Rational(q.numerator,q.denominator))))<=q**N+(1-q)**N,'actual aggregate tail bound controls')

def C(n,k): return comb(n,k) if 0<=k<=n else 0
def donor(m):
    """Whole annulus bisection minus actual prescribed heads, not signed defect."""
    M=2*m; W=[0]*m; V=[0]*m
    for weight in range(1,M):
        total=C(M-m,weight)+C(M-m,weight-m)
        preset=sum(C(M-n-1,weight-1)+C(M-n-1,weight-n)
                   for n in range(m+1,M) if (n-m)%2)
        gate(total%2==0,'entire Hamming layer can be bisected')
        needed=total//2-preset
        zero_capacity=C(m-1,weight-1); one_capacity=C(m-1,weight-m)
        gate(0<=needed<=zero_capacity+one_capacity,'literal remaining donor head count')
        take0=min(needed,zero_capacity); take1=needed-take0
        if 0<=weight-1<m: W[weight-1]=take0
        if 0<=M-1-weight<m: V[M-1-weight]=take1
    gate(all(0<=v<=comb(m-1,k) for row in (W,V) for k,v in enumerate(row)),'both native donor orientation boxes')
    return coeffs(W),coeffs(V)
rows={}
def orientation(stop):
    Fout=[0]*(2*stop+1); Gout=Fout.copy()
    for m in range(1,stop):
        for target,row in zip((Fout,Gout),rows[m]):
            for j,c in enumerate(row):
                target[m+j]+=c; target[m+j+1]-=c
    return Fout,Gout
def reflected(v):
    return [sum(v[k]*comb(k,j)*(-1)**j for k in range(j,len(v))) for j in range(len(v))]
for m in (1,2,4,8,16,32):
    rows[m]=donor(m)
    for n in range(m+1,2*m): rows[n]=([int((n-m)%2)],)*2
    f,g=orientation(2*m); gr=reflected(g)
    target=[F((int(k==0)-int(k==2*m))-(-1)**k*C(2*m,k),2) for k in range(len(f))]
    gate([a+b for a,b in zip(f,gr)]==target,'complete exact dyadic-prefix fairness')

def det_mod(matrix,p):
    """Fraction-free determinant recurrence over the finite field."""
    A=[[v%p for v in row] for row in matrix]; sign=1; last=1; n=len(A)
    for k in range(n-1):
        pivot=next((i for i in range(k,n) if A[i][k]),None)
        if pivot is None: return 0
        if pivot!=k: A[k],A[pivot]=A[pivot],A[k]; sign=-sign
        v=A[k][k]; inverse=pow(last,-1,p)
        for i in range(k+1,n):
            for j in range(k+1,n): A[i][j]=(v*A[i][j]-A[i][k]*A[k][j])*inverse%p
            A[i][k]=0
        last=v
    return sign*A[-1][-1]%p
f,g=orientation(64)
determinants=[]
for n in (2,4,8,12,16,24,30):
    ds=[det_mod([[row[i+j+1] for j in range(n)] for i in range(n)],101) for row in (f,g)]
    gate(ds[0]!=0 and (ds[1]!=0 if n>2 else ds[1]==0),'independent Hankel determinant control')
    if n==2: gate(any(g[k] for k in (1,2,3)),'exceptional 2x2 second-orientation rank is one')
    determinants.append([n,*ds])

# Complete original eight-bit spill coloring, classified by native Hamming weights.
heads=[0]*9; total=[0]*9; signed={2:[0]*9,4:[0]*9}
for w in product((0,1),repeat=8):
    if all(bit==w[0] for bit in w): continue
    m=next(i for i in range(1,8) if w[i]!=w[0]); weight=sum(w)
    head=(w[0]==0) if m in (1,4,7) else (w[0]==1) if m in (2,6) else w[4]==0 if m==3 else w[6]==1
    total[weight]+=1; heads[weight]+=head
    if m>=2: signed[2 if m<4 else 4][weight]+=2*head-1
gate(all(2*a==b for a,b in zip(heads,total)),'literal full-eight-bit composition fairness')
gate(signed[2]==[0,0,-2,-4,0,4,2,0,0] and signed[4]==[-x for x in signed[2]],'nonzero shell residuals cancel only after handoff')
print('AFFINE CONTENT: exact substitution determinant and nonunit slope hostile PASS')
print('POLYNOMIAL EXTENSION: independent nonconstant compositions',composition_controls)
print('LITERAL BERNSTEIN BOXES',boxes)
print('STATIONARY orientation sums:',str(FF),str(GG),'bias 1/3 ->16/35')
print('VON NEUMANN: reflected half-integral germs fair; no finite critical deadline')
print('AGGREGATE: correct parity/tails do not restore first-row signed capacity')
print('DYADIC PREFIXES through64 and independent Hankel determinants mod101:',json.dumps(determinants,separators=(',',':')))
print('FULL SPILL Hamming head counts',heads,'opposite shell residuals',signed[2],signed[4])
print('PASS',G,'always-active independent exact gates; raw LF')
print('Scope: rational orientation-resolved polynomial architecture under finite linear deadline; no bound on optimal slope.')
