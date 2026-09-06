"""Independent shared-root wall audit: quotient-field arithmetic and signs."""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve().parent
STEM="continuing6_20260906_shared_wall"
GATES=0
r,x,y,z,u,s,v,t=S.symbols("r x y z u s v t")


def gate(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise ArithmeticError(label)


def same(a,b,label):
    gate(S.cancel(a-b)==0,label)


def locate(name):
    for path in (HERE/name,HERE.parent/"05-knowledge"/"results"/name):
        if path.exists():
            return path
    raise FileNotFoundError(name)


for suffix,expected in ((".py","ef604a3c9276171b844f6e260bc788a17160ad6e85a70f3fbd131ee4f3c362c4"),
                        (".out","4803b7577514138cbe0fca4ef339b0e2cb690ba421e2276eb823846ef0d33572"),
                        ("_certificate.json","6adab7b43f30296cdac428596ce7c046d163751b6d48fea4db02fe17a100dce4")):
    gate(hashlib.sha256(locate(STEM+suffix).read_bytes()).hexdigest()==expected,"producer raw pin "+suffix)
certificate=json.loads(locate(STEM+"_certificate.json").read_bytes())

B=v**5-13*v**4+55*v**3-x*v*v+y*v-z
C=v**4-12*v**3+45*v*v-S.Rational(2,3)*x*v+S.Rational(3,7)*y
D=v**4-11*v**3+36*v*v-S.Rational(5,12)*x*v+S.Rational(1,7)*y
phase=S.Rational(1,11)-10*s+x*s*s-S.Rational(12,7)*y*s**3+z*s**4
RC=705672*r**4-15079284*r**3+115450055*r*r-375696750*r+441317250
RD=120434688*r**4-2310536448*r**3+15142183583*r*r-39453584808*r+35022385200
charts={}
HC=44*r*r*u**4-132*r*r*u**3-495*r*u**4+1584*r*u**3-3*r+1650*u**4-5940*u**3+330*u
HD=66*r*r*u**4-132*r*r*u**3-704*r*u**4+1452*r*u**3-r+2167*u**4-4752*u**3+110*u
for name,interlacer,quad,expected in (("C",C,5*u*u-24*u+9,RC),("D",D,23*u*u-60*u+12,RD)):
    sol=S.solve([B.subs(v,r),interlacer.subs(v,r)],[y,z])
    same(B.subs(v,r).subs(sol),0,"exact shared beta chart "+name)
    same(interlacer.subs(v,r).subs(sol),0,"exact shared interlacer chart "+name)
    expression=S.factor(phase.subs(sol).subs(s,u/r))
    denominator=9 if name=="C" else 12
    same(S.diff(expression,x),u*u*quad/(denominator*r*r),"vanishing native x coefficient "+name)
    H=HC if name=="C" else HD
    same(expression,x*u*u*quad/(denominator*r*r)-H/((33 if name=="C" else 11)*r),"complete constant and coefficient chart "+name)
    constant=S.together(expression.subs(x,0)).as_numer_denom()[0]
    linear=S.rem(constant,quad,u)
    a=S.expand(linear).coeff(u,1)
    b=S.expand(linear).coeff(u,0)
    result=S.Poly(S.resultant(linear,quad,u),r).clear_denoms()[1].primitive()[1]
    if result.LC()<0:
        result=-result
    same(result.as_expr(),expected,"primitive singular-wall eliminant "+name)
    gate(S.gcd(expected,a)==1,"no singular rational reconstruction denominator "+name)
    gate(S.Poly(expected,r).count_roots(0,S.oo)==4,"all four wall roots are positive "+name)
    gate(S.gcd(expected,S.diff(expected,r))==1,"all four wall roots are simple "+name)
    recovered=S.cancel(-b/a)
    numerator=S.together(quad.subs(u,recovered)).as_numer_denom()[0]
    same(S.rem(numerator,expected,r),0,"no resultant extraneous branch "+name)
    gate(S.Poly(quad,u).count_roots(0,S.oo)==2,"both possible u branches positive "+name)
    charts[name]=(sol,recovered)

gate(len(certificate["singular_pairs"])==8,"all eight singular pair isolations retained")
for name in ("C","D"):
    entries=[row for row in certificate["singular_pairs"] if row[0]==name]
    gate(len(entries)==4,"four complete interval owners "+name)
    previous=S.Integer(0)
    expected=RC if name=="C" else RD
    for _,left,right,ulo,uhi in entries:
        a,b=S.Rational(left),S.Rational(right)
        gate(0<=previous<a<b,"all singular r isolations positive and disjoint "+name)
        gate(S.Poly(expected,r).count_roots(a,b)==1,"each published r isolation owns exactly one root "+name)
        gate(0<S.Rational(ulo)<S.Rational(uhi),"recorded positive u enclosure "+name)
        previous=b
    gate(S.Rational(entries[-1][1])>S.Rational(71,10),"largest singular root outside beta cap "+name)
gate((5*r*r-26*r-67).subs(r,S.Rational(71,10))>0,"anchor Cauchy bound excludes r>=71/10")

# The selected algebraic point is evaluated in QQ[r]/RC, without trusting
# decimal coordinates or substituting rounded coefficients into root tests.
R=S.Poly(RC,r,domain=S.QQ)


def red(a):
    return S.rem(S.expand(a),RC,r)


def mul(a,b):
    return red(a*b)


def inv(a):
    return S.invert(a,RC,r)


def algebraic(a):
    num,den=S.cancel(a).as_numer_denom()
    return mul(red(num),inv(red(den)))


lo=F(3528453826758,10**12)
hi=F(3528453826759,10**12)
gate(R.count_roots(S.Rational(lo),S.Rational(hi))==1,"unique selected algebraic r")
gate(R.eval(S.Rational(lo))*R.eval(S.Rational(hi))<0,"opposite exact isolating signs")
for _ in range(40):
    midpoint=(lo+hi)/2
    mvalue=R.eval(S.Rational(midpoint))
    gate(mvalue!=0,"exact bisection midpoint is not the algebraic root")
    if R.eval(S.Rational(lo))*mvalue<0:
        hi=midpoint
    else:
        lo=midpoint
gate(R.eval(S.Rational(lo))*R.eval(S.Rational(hi))<0,"refined exact isolating signs")
xx=red(S.Rational(3,5)*r*(4*r*r-45*r+150)+S.Rational(1,10))
yy=red(charts["C"][0][y].subs(x,xx))
zz=red(charts["C"][0][z].subs(x,xx))
uu=algebraic(charts["C"][1])
ss=mul(uu,inv(r))
same(zz,r*r/18,"strict positive z construction")
same(red(5*uu*uu-24*uu+9),0,"selected u is exactly a coefficient-zero root")
same(algebraic(phase.subs({x:xx,y:yy,z:zz,s:ss})),0,"original phase exactly vanishes in the field")
for key,expr in (("x",xx),("y",yy),("z",zz)):
    same(red(S.sympify(certificate[key],locals={"r":r})-expr),0,"frozen algebraic coordinate "+key)


def iadd(A,B):
    return A[0]+B[0],A[1]+B[1]


def imul(A,B):
    values=[a*b for a in A for b in B]
    return min(values),max(values)


def enclosure(poly):
    out=(F(0),F(0))
    for c in S.Poly(poly,r,domain=S.QQ).all_coeffs():
        q=F(int(c.p),int(c.q))
        out=iadd(imul(out,(lo,hi)),(q,q))
    return out


for name,expr in (("x",xx),("y",yy),("z",zz),("s",ss)):
    gate(enclosure(expr)[0]>0,"strict positive field coordinate "+name)
uI=enclosure(uu)
gate(0<uI[0]<uI[1]<1,"selected u is the lower root (12-3sqrt11)/5")


def evaluate(coefficients,point):
    out=S.Integer(0)
    for c in reversed(coefficients):
        out=red(point*out+c)
    return out


def remove_root(coefficients):
    # Synthetic division in the quotient ring, ascending input/output.
    descending=[coefficients[-1]]
    for c in reversed(coefficients[1:-1]):
        descending.append(red(c+r*descending[-1]))
    same(red(coefficients[0]+r*descending[-1]),0,"known shared root divides natively")
    return list(reversed(descending))


BB=[-zz,yy,-xx,55,-13,1]
CC=[S.Rational(3,7)*yy,-S.Rational(2,3)*xx,45,-12,1]
DD=[S.Rational(1,7)*yy,-S.Rational(5,12)*xx,36,-11,1]
Bred=remove_root(BB)
Cred=remove_root(CC)
# Broad independent rational brackets, distinct from producer narrow decimals.
brackets={"Bred":[(18,20),(669,671),(2450,2453),(6330,6334)],
          "Cred":[(387,390),(1950,1954),(6129,6133)],
          "D":[(181,184),(1492,1495),(3382,3386),(5938,5942)]}
polys={"Bred":Bred,"Cred":Cred,"D":DD}
sign_rows=[]
for name,intervals in brackets.items():
    gate(len(intervals)==len(polys[name])-1,"sign brackets exhaust polynomial degree "+name)
    for a,b in intervals:
        left=enclosure(evaluate(polys[name],S.Rational(a,1000)))
        right=enclosure(evaluate(polys[name],S.Rational(b,1000)))
        sa=1 if left[0]>0 else -1 if left[1]<0 else 0
        sb=1 if right[0]>0 else -1 if right[1]<0 else 0
        gate(sa*sb==-1,"uniform rational sign bracket "+name)
        sign_rows.append([name,a,b,sa,sb])
# The full ordering includes the exact known r root of both B and C.
bints=[(F(a,1000),F(b,1000)) for a,b in brackets["Bred"]]
bints.insert(3,(lo,hi))
cints=[(F(a,1000),F(b,1000)) for a,b in brackets["Cred"]]
cints.insert(2,(lo,hi))
dints=[(F(a,1000),F(b,1000)) for a,b in brackets["D"]]
gate(all(bints[i][1]<bints[i+1][0] for i in range(4)),"five distinct positive beta roots")
for i in range(4):
    gate(bints[i][1]<dints[i][0]<dints[i][1]<bints[i+1][0],"strict D interlacing in native order")
    if i==2:
        gate(cints[i]==bints[3],"one exact shared C root")
        gate(bints[i][1]<cints[i][0],"shared C root lies at the correct weak interval endpoint")
    else:
        gate(bints[i][1]<cints[i][0]<cints[i][1]<bints[i+1][0],"other C roots strictly interlace")

# Independent original Hadamard convolution, evaluated directly in the field.
O={j:S.Integer(comb(14,2*j+1)) for j in range(7)}
E={j:S.Integer(comb(14,2*j)) for j in range(8)}
beta=dict(zip(range(-1,5),[1,13,55,xx,yy,zz]))
cr=dict(zip(range(-1,4),[1,12,45,S.Rational(2,3)*xx,S.Rational(3,7)*yy]))
dr=dict(zip(range(-1,4),[1,11,36,S.Rational(5,12)*xx,S.Rational(1,7)*yy]))


def convolution(A,B):
    out={}
    for i,a in A.items():
        for j,b in B.items():
            out[i+j]=red(out.get(i+j,0)+mul(a,b))
    return out


qa=convolution(O,O)
for j,a in convolution(E,E).items():
    qa[j-1]=red(qa.get(j-1,0)+a)
qb=convolution(beta,beta)
for j,b in convolution(cr,dr).items():
    qb[j+1]=red(qb.get(j+1,0)+2*b)
qcoeff={j:mul(qa[j],qb[j]) for j in qa.keys()&qb.keys()}
same(qcoeff[-1],28,"complete inverse response carry")
original_P=S.Integer(0)
original_Q=S.Integer(0)
power=inv(-ss)
for j in range(-1,9):
    if j>=0 and j in O and j in beta:
        original_P=red(original_P+mul(mul(O[j],beta[j]),power))
    if j in qcoeff:
        original_Q=red(original_Q+mul(qcoeff[j],power))
    power=mul(power,-ss)
same(original_P,0,"literal original Hadamard P vanishes")
qI=enclosure(original_Q)
gate(qI[0]>-5634725 and qI[1]<-5634723,"independent quotient-field original response enclosure")
print("Both exact shared charts reconstructed; each has four simple positive singular r values, with no extraneous reconstruction branch.")
print("Selected r uniquely isolated in(3528453826758,3528453826759)/10^12; lower u branch, x,y,z,s positive.")
print("The referee refines r by40 exact sign bisections before evaluating the quotient-field response.")
for row in sign_rows:
    print("SIGN",json.dumps(row,separators=(",",":")))
print("B has five simple positive roots; C shares exactly one and otherwise strictly interlaces; D strictly interlaces.")
print("Original native phase coefficient and phase vanish; Q(-s) lies strictly between-5634725 and-5634723.")
print("Division by the shared-chart phase coefficient would discard an admissible full-model wall. General shared sign remains open.")
print("PASS",GATES,"always-active exact gates.")
