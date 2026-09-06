"""Independent symbolic fibers/windows and outward-dyadic elimination audit."""
from pathlib import Path
from fractions import Fraction as F
from math import factorial, comb, isqrt
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve().parent
STEM="continuing2_20260906_laurent_sparse_amplitude"
gates=0
def need(ok,label):
    global gates
    gates+=1
    if not ok:
        raise ArithmeticError(label)

pins={".py":"2afc7df14b92e70a700a9a06a17b34ee2331eabde89902dcf408f597dccb9e28",
      ".out":"7b95d9102c18309ce245d9017f3456b0f75447117f1ec6ac363393c63f183777",
      "_certificate.json":"df8471af42220b5b26f7e4ff78dd9611d7999c0578e36603a6a7a33ba3c0afa2"}
for ending,pin in pins.items():
    path=HERE/(STEM+ending)
    if not path.exists() and ending==".out":
        path=HERE.parent/"05-knowledge/results"/(STEM+ending)
    need(hashlib.sha256(path.read_bytes()).hexdigest()==pin,"frozen producer "+ending)
cert=json.loads((HERE/(STEM+"_certificate.json")).read_text())
s,u,v=S.symbols("s u v")
def P(expr):return S.Poly(S.expand(expr),s,domain=S.QQ)
def cs(expr,n):return [P(expr).nth(i) for i in range(n)]
def listed(xs):return [S.Rational(x) for x in xs]
def ffall(n,r):return factorial(n)//factorial(n-r)

# Literal multinomial rows, without prescribed j or any producer source.
rows={}
for mass in range(1,29):
    fiber=[]
    for a in range(mass+1):
        for b in range(mass-a+1):
            c=mass-a-b
            if -27*a+b+15*c==0:
                fiber.append((a,b,c,factorial(mass)//(factorial(a)*factorial(b)*factorial(c))))
    need(bool(fiber)==(mass in (14,28)),"literal first-return clock")
    rows[mass]=fiber
first=sum(w*(-s)**(a-1) for a,b,c,w in rows[14])
target=sum(w*(-s)**(a-2)*s for a,b,c,w in rows[28])
first,target=map(S.expand,(first,target))
need(len(rows[14])==5 and len(rows[28])==10,"complete channels")
need(all((a,b,c)==(1+j,12-3*j,1+2*j) for j,(a,b,c,w) in enumerate(rows[14])),"first original phase monomial")
need(all((a,b,c)==(2+j,24-3*j,2+2*j) for j,(a,b,c,w) in zip(range(-1,9),rows[28])),"doubled original carry")
p=P(first).monic()
need(S.expand(first-2002*p.as_expr())==0,"original first normalization")
need(cs(p.as_expr(),5)==listed(cert["P"]),"quartic first equation")
need(cs(target,10)==listed(cert["full_target"]),"literal full target scale")
need(P(target).degree()+1==len(cert["full_target"]),"target comparison includes the full degree")
need(P(target).nth(0)<0,"negative carried endpoint retained")

# Literal north/east dynamic-program counts for the complete Laurent beta path.
def pathcount(x,y):
    grid=[[0]*(y+1) for _ in range(x+1)]
    grid[0][0]=1
    for a in range(x+1):
        for b in range(y+1):
            if a or b:grid[a][b]=(grid[a-1][b] if a else 0)+(grid[a][b-1] if b else 0)
    return grid[x][y]
kernel=sum((-s)**(j+1)*pathcount(12-3*j,1+j)*u**(2*(4-j)) for j in range(-1,5))
H=S.expand((1+u)**14*kernel)
need(S.Poly(H,u).nth(9)==S.expand(-s*first),"actual selected zero")
beta=S.Poly(sum(S.Poly(kernel.subs(s,1),u).nth(2*i)*v**i for i in range(6)),v)
need(beta.as_expr()==v**5-13*v**4+55*v**3-84*v**2+35*v-1,"complete beta kernel")
need(beta.count_roots(0,S.oo)==5,"all actual beta roots positive")
windows=[]
for r in range(9):
    derivative=S.Poly(S.diff(H,u,r),u)
    deg=18-2*r
    raw=sum(derivative.nth(i)*derivative.nth(deg-i) for i in range(deg+1))
    jr=S.cancel(raw/s)
    need(not S.denom(jr).has(s),"exact positive phase removal")
    need(P(jr).degree()+1==len(cert["full_windows"][r]),"window comparison includes the full degree")
    need(cs(jr,len(cert["full_windows"][r]))==listed(cert["full_windows"][r]),"literal full derivative window")
    remainder=P(jr).rem(p)
    need(cs(remainder.as_expr(),4)==listed(cert["residues"][r]),"exact window remainder")
    windows.append(remainder)
tr=P(target).rem(p)
need(cs(tr.as_expr(),4)==listed(cert["target"]),"exact target remainder")
need(cert["labels"]==[[0,0],[1,7],[1,8],[3,2]],"specified four atoms")
norms=[ffall(9,r)**2 for d,r in cert["labels"]]
need(norms==listed(cert["normalizers"]),"positive stated normalizations")

# Outward dyadic interval arithmetic. Unlike producer's permutation expansion,
# determinants are obtained by interval Gaussian elimination with legal pivots.
SCALE=2**240
def down(x):return F((x*SCALE).numerator//(x*SCALE).denominator,SCALE)
def up(x):return -down(-x)
def interval(a,b=None):
    a=F(a);b=a if b is None else F(b)
    need(a<=b,"ordered interval")
    return (down(a),up(b))
def add(a,b):return (down(a[0]+b[0]),up(a[1]+b[1]))
def neg(a):return (-a[1],-a[0])
def sub(a,b):return add(a,neg(b))
def mul(a,b):
    vals=[x*y for x in a for y in b]
    return (down(min(vals)),up(max(vals)))
def div(a,b):
    need(b[0]*b[1]>0,"interval divisor avoids zero")
    return mul(a,(down(1/b[1]),up(1/b[0])))
def power(a,n):
    out=(F(1),F(1))
    for _ in range(n):out=mul(out,a)
    return out
def eval_interval(poly,root):
    out=(F(0),F(0))
    for coeff in poly.all_coeffs():out=add(mul(out,root),interval(F(coeff)))
    return out
def sqrt_enclosure(a,b):
    scale=2**210
    lo=isqrt((a*scale*scale).numerator//(a*scale*scale).denominator)
    hi=isqrt((b*scale*scale).numerator//(b*scale*scale).denominator)+1
    result=(F(lo,scale),F(hi,scale))
    need(result[0]**2<=a and result[1]**2>=b,"independent dyadic square-root enclosure")
    return result
def determinant(matrix):
    a=[list(row) for row in matrix]
    det=(F(1),F(1))
    for i in range(len(a)):
        legal=[r for r in range(i,len(a)) if a[r][i][0]*a[r][i][1]>0]
        need(bool(legal),"nonzero interval elimination pivot")
        row=max(legal,key=lambda r:min(abs(x) for x in a[r][i]))
        if row!=i:
            a[row],a[i]=a[i],a[row]
            det=neg(det)
        pivot=a[i][i]
        det=mul(det,pivot)
        for r in range(i+1,len(a)):
            ratio=div(a[r][i],pivot)
            for c in range(i+1,len(a)):
                a[r][c]=sub(a[r][c],mul(ratio,a[i][c]))
            a[r][i]=(F(0),F(0))
    return det

brackets=[tuple(F(x) for x in row) for row in cert["intervals"]]
need(all(0<a<b for a,b in brackets) and all(brackets[i][1]<brackets[i+1][0] for i in range(3)),"four disjoint positive root brackets")
matrix=[];rhs=[]
for i,(a,b) in enumerate(brackets):
    need(p.eval(S.Rational(a))*p.eval(S.Rational(b))<0,"one simple root per bracket, exhaustive degree")
    # Fresh refinement compensates for the independent elimination path's
    # deliberately looser dependency enclosure; not a change to producer data.
    for _ in range(64):
        mid=(a+b)/2
        value=p.eval(S.Rational(mid))
        need(value!=0,"nonrational refinement midpoint")
        if p.eval(S.Rational(a))*value<0:b=mid
        else:a=mid
    root=interval(a,b)
    sqrt=sqrt_enclosure(a,b)
    js=[eval_interval(w,root) for w in windows]
    tv=eval_interval(tr,root)
    need(all(j[1]<0 for j in js) and tv[1]<0,"all nine windows and actual target strictly negative")
    row=[]
    for (d,r),norm in zip(cert["labels"],norms):
        if d==0 and r==0:
            atom=(F(1),F(1))
        else:atom=div(mul(power(sqrt,d),js[r]),mul(interval(norm),js[0]))
        need(atom[0]>0,"each normalized matrix entry positive")
        row.append(atom)
    matrix.append(row)
    rhs.append(div(tv,js[0]))

det=determinant(matrix)
need(det[1]<0,"main determinant strictly negative")
det_bounds=[(-371065,-371064),(-365451,-365450),(-156491,-156490),(-1845480,-1845479),(-1916,-1915)]
coeff_bounds=[(984869,984870),(421733,421734),(4973466,4973467),(5162,5163)]
need(F(det_bounds[0][0],10**6)<det[0]<=det[1]<F(det_bounds[0][1],10**6),"published determinant outer bounds")
for col in range(4):
    replaced=[[rhs[r] if c==col else matrix[r][c] for c in range(4)] for r in range(4)]
    numerator=determinant(replaced)
    need(numerator[1]<0,"negative Cramer numerator")
    value=div(numerator,det)
    need(value[0]>0,"strict positive common algebraic coefficient")
    dl,dh=det_bounds[col+1]
    cl,ch=coeff_bounds[col]
    need(F(dl,10**6)<numerator[0]<=numerator[1]<F(dh,10**6),"published replacement determinant bounds")
    need(F(cl,10**6)<value[0]<=value[1]<F(ch,10**6),"published coefficient bounds independently enclosed")
    print("CRAMER",col,"coefficient in",F(cl,10**6),F(ch,10**6),"numerator in",F(dl,10**6),F(dh,10**6))

# Irreducibility and the rational/full-eight-root firewall.
mod=S.Poly(11*p.as_expr(),s, modulus=3).monic()
need(mod.monic().as_expr()==s**4-s-1,"correct monic reduction modulo3")
need(all(mod.eval(j)%3 for j in range(3)),"no linear factor overF3")
for a in range(3):
    for b in range(3):
        q=S.Poly(s*s+a*s+b,s,modulus=3)
        need(not mod.rem(q).is_zero,"no monic quadratic divisor overF3")
need(p.nth(0)==S.Rational(1,11),"quartic root norm is nonsquare rational")
need(isqrt(11)**2!=11,"norm obstruction to square root in quartic field")
need(S.Poly(p.as_expr().subs(s,u*u),u).is_irreducible,"independent full-eight degree check")
print("FIELD p irreducible; norm r=1/11 nonsquare; p(z^2) irreducible degree8")
print("SCOPE fixed actual row; common positive algebraic coefficients on four positive branches; rational odd-atom/all-eight identities excluded")
print("PASS",gates,"always-active gates; literal fibers, symbolic windows, outward-dyadic elimination; no producer import")
