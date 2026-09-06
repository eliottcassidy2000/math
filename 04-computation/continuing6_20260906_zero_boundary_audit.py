"""Independent zero-beta referee: dense native carriers and Sturm charts."""
from pathlib import Path
from math import comb
import hashlib
import json
import sys
import sympy as S

sys.stdout.reconfigure(newline="\n")
HERE=Path(__file__).resolve().parent
STEM="continuing6_20260906_zero_boundary"
GATES=0
x,y,s,t,w,v,u,X=S.symbols("x y s t w v u X")


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


for suffix,expected in ((".py","dd28188852e5fb955d19d886687f2cb791e25cbf1a6f1c3b2bf66158cc467e6b"),
                        (".out","6fe1ddf5e2182828405968b34c8236bc8ccea8f5b046e843781dcf26362fef0c"),
                        ("_certificate.json","f6c169d4190f96465f93d195fe4517537a9c40e85370ae6a736e0e2886e9b12e")):
    gate(hashlib.sha256(locate(STEM+suffix).read_bytes()).hexdigest()==expected,"frozen producer pin "+suffix)


# Dense Laurent expressions reconstruct each original coefficient, including
# the odd/even binomial split and the exponent-minus-one inverse carry.
binomial=S.Poly((1+u)**14,u)
O=sum(binomial.nth(2*j+1)*t**j for j in range(7))
E=sum(binomial.nth(2*j)*t**j for j in range(8))
beta=1/t+13+55*t+x*t**2+y*t**3
cr=1/t+12+45*t+S.Rational(2,3)*x*t**2+S.Rational(3,7)*y*t**3
dr=1/t+11+36*t+S.Rational(5,12)*x*t**2+S.Rational(1,7)*y*t**3
left=S.expand(O**2+E**2/t)
right=S.expand(beta**2+2*t*cr*dr)
Q=sum(S.expand(left).coeff(t,j)*S.expand(right).coeff(t,j)*t**j for j in range(-1,14))
P=sum(S.expand(O).coeff(t,j)*S.expand(beta).coeff(t,j)*t**j for j in range(7))
p=S.expand(P.subs(t,-s)/2002)
T=S.expand(s*Q.subs(t,-s))
same(p,-S.Rational(12,7)*y*s**3+x*s**2-10*s+S.Rational(1,11),"literal original phase normalization")
same(S.expand(Q).coeff(t,-1),28,"inverse carry stays28")
gate(S.Poly(T,s).degree()==8,"original mixed top response survives zero beta")

B=v**5-13*v**4+55*v**3-x*v**2+y*v
C=v**4-12*v**3+45*v**2-S.Rational(2,3)*x*v+S.Rational(3,7)*y
D=v**4-11*v**3+36*v**2-S.Rational(5,12)*x*v+S.Rational(1,7)*y
series=S.series((C/B).subs(v,1/w)/w,w,0,5).removeO().expand()
mu=[series.coeff(w,j) for j in range(5)]
H2=S.Matrix(2,2,lambda i,j:mu[i+j+1]).det()
H3=S.Matrix(3,3,lambda i,j:mu[i+j]).det()
cap=S.Rational(7,72)*(x-75)*(135-x)
same(H2,(x-75)/3,"independent shifted C Gram constraint")
same(H3,(x-75)*(135-x)/9-S.Rational(8,7)*y,"independent ordinary C Gram constraint")
same(cap,S.Rational(175,2)-S.Rational(7,72)*(x-105)**2,"packet rectangle and maximum cap")

# Solve the original phase as an affine equation in y, independently from the
# producer's displayed eliminant, and keep s>0 in the later analytic proof.
Y=S.solve(p,y)[0]
same(p.subs(y,Y),0,"original phase is exactly retained")
F=S.cancel(-968*T.subs(y,Y)).expand()
J=S.cancel(S.Rational(72,7)*s**3*(Y-cap)).expand()
gate(S.Poly(F,x,s).total_degree()==8 and len(S.Poly(F,x,s).terms())==11,"complete restricted response support")
same(J,s**3*(x*x-210*x+10125)+6*x*s*s-60*s+S.Rational(6,11),"native packet inequality after phase elimination")
A=S.Poly(F,x).nth(2)
same(A,264385*s**5*(19601*s+31920),"positive quadratic leading coefficient")
same(S.diff(T,x,2),397670*s**5*(66-85*s),"first-branch x convexity")
same(S.diff(T,y,2),S.Rational(4011660,7)*s**7*(140-13*s),"first-branch y convexity")

cert=json.loads(locate(STEM+"_certificate.json").read_bytes())
symbols={"x":x,"y":y,"s":s}
for key,expr in (("response_F",F),("response_T",T),("phase_p",p),("necessary_J",J)):
    same(S.sympify(cert[key],locals=symbols),expr,"frozen full expression "+key)

L1,R1=S.Rational(1,110),S.Rational(1,90)
L2,M2,R2=S.Rational(63,1000),S.Rational(11,100),S.Rational(1,8)
L3,R3=S.Rational(9,16),S.Rational(3,5)
gate(p.subs({x:75,y:S.Rational(175,2),s:L1})>0,"uniform first left endpoint")
gate(p.subs({x:135,y:0,s:R1})<0,"uniform first right endpoint")
gate(p.subs({x:135,y:0,s:L2})<0,"uniform second left endpoint")
for a in (R2,L3):
    q=S.expand(p.subs({y:cap,s:a}))
    centre=105-3/a
    low=q.subs(x,centre)
    same(q,low+a**3*(x-centre)**2/6,"independent completed-square endpoint")
    gate(low>0,"uniform positive middle/third endpoint")
gate(0<L1<R1<L2<R2<L3<R3,"three disjoint phase branches")

charts={}
for xx in (75,135):
    for yy in (0,S.Rational(175,2)):
        charts[f"raw first corner x={xx},y={yy}"]=(-T.subs({x:xx,y:yy}),L1,R1)
disc=S.discriminant(F,x)
charts["negative discriminant"]=(-S.cancel(disc/s**4),L2,M2)
charts["J at x=105"]=(J.subs(x,105),M2,R2)
charts["F at x=105"]=(F.subs(x,105),M2,R2)
charts["negative Fx at x=105"]=(-S.diff(F,x).subs(x,105),M2,R2)
charts["J at x=84"]=(J.subs(x,84),L3,R3)


def variation(seq,point):
    signs=[]
    for poly in seq:
        value=poly.eval(point)
        if value:
            signs.append(1 if value>0 else -1)
    return sum(a!=b for a,b in zip(signs,signs[1:]))


sturm_rows=[]
coefficient_count=0
for item in cert["certificates"]:
    label=item["label"]
    if "interval" in item:
        raw,lo,hi=charts[label]
        poly=S.Poly(raw,s,domain=S.QQ)
        gate([str(lo),str(hi)]==item["interval"],"exact certificate chart domain")
        # A separate proof of each sign: exact rational Sturm variation count.
        seq=[S.Poly(z,s,domain=S.QQ) for z in S.sturm(poly.as_expr(),s)]
        va,vb=variation(seq,lo),variation(seq,hi)
        gate(poly.eval(lo)>0 and poly.eval(hi)>0,"positive chart endpoint signs")
        gate(va-vb==0,"independent Sturm exclusion of every interior root")
        sturm_rows.append([label,poly.degree(),va,vb])
        degree=item["degree"]
        gate(degree>=poly.degree(),"homogeneous degree bound")
        same(poly.as_expr(),sum(S.Rational(c)*s**i for i,c in enumerate(item["power_coefficients"])),"complete native power coefficients")
        transformed=S.Poly(sum(poly.nth(i)*(lo+hi*u)**i*(1+u)**(degree-i) for i in range(degree+1)),u)
        for i,c in enumerate(item["homogeneous_coefficients"]):
            gate(S.Rational(c)==transformed.nth(i)>0,"independent positive homogeneous coefficient")
            coefficient_count+=1
    else:
        x0,s0=map(S.Rational,item["quadrant_base"])
        gate((x0,s0) in ((84,L3),(75,R3)),"only the two necessary quadrant bases")
        shifted=S.Poly(F.subs({x:X+x0,s:u+s0}).expand(),X,u,domain=S.QQ)
        gate(len(item["coefficients"])==21,"complete quadratic-by-sextic rectangle")
        seen=set()
        for record in item["coefficients"]:
            i,j=record["x_degree"],record["s_degree"]
            gate((i,j) not in seen and 0<=i<=2 and 0<=j<=6,"each quadrant monomial checked once")
            seen.add((i,j))
            gate(shifted.coeff_monomial(X**i*u**j)==S.Rational(record["coefficient"])>0,"positive native quadrant coefficient")
            coefficient_count+=1
        gate(len(seen)==21,"quadrant retains every monomial including constant")
gate(len(sturm_rows)==9 and coefficient_count==107,"entire finite sign certificate covered")

same(S.diff(J,x),s*s*(2*s*(x-105)+6),"phase packet derivative")
same(S.diff(F,x),2*A*(x-105)+S.diff(F,x).subs(x,105),"response derivative on the admissible left side")
gate(66-85*R1>0 and 140-13*R1>0,"separate convexity holds throughout first branch")
gate(6-42*L3<0,"third branch J strictly decreases at every x<=84")

# Actual full model and the intentionally larger C-only boundary.
numeric={x:84,y:35}
gate(S.Poly((B/v).subs(numeric),v).count_roots(0,S.oo)==4,"four positive nonzero beta nodes")
for label,numerator in (("C",C),("D",D)):
    numeric_series=S.series((numerator/B).subs(numeric).subs(v,1/w)/w,w,0,9).removeO().expand()
    mm=[numeric_series.coeff(w,j) for j in range(9)]
    H=S.Matrix(5,5,lambda i,j:mm[i+j])
    for size in range(1,6):
        gate(H[:size,:size].det()>0,"literal full-model residue Gram positivity "+label)
same(B.subs({x:75,y:0}),v*v*(v-3)*(v-5)**2,"C-only repeated beta hostile")
same((C/B).subs({x:75,y:0}),S.Rational(2,3)/v+S.Rational(1,3)/(v-3),"positive C-only residue measure")
gate(D.subs({x:75,y:0,v:5})!=0,"C-only extension really exceeds two-interlacer scope")
hostile={x:84,y:S.Rational(1050,11),s:S.Rational(1,6)}
same(p.subs(hostile),0,"without-packet hostile retains original phase")
same(T.subs(hostile),S.Rational(228261457669,313632),"positive raw response outside packet")
same(H3.subs(hostile),-S.Rational(639,11),"precise failed C Gram inequality")
off={x:84,y:35,s:S.Rational(1,4)}
same(T.subs(off),S.Rational(412444713751,32768),"off-phase response may be positive inside full geometry")
same(p.subs(off),S.Rational(335,176),"original phase equation cannot be discarded")

print("Independent dense binomial14 carriers: original phase and inverse carry28 retained.")
print("Two C/B Gram inequalities give the complete packet rectangle; original phase elimination matches F,J.")
for row in sturm_rows:
    print("STURM",json.dumps(row,separators=(",",":")))
print("All107 coefficient signs and both quadrant arrays verified independently; nine charts have no roots by Sturm.")
print("All-phase coverage and derivative orientations checked; full-model, C-only repeated, off-phase and no-packet controls retained.")
print("Audit scope: original positive phases on zero beta boundary; shared-root and general sign remain open.")
print("PASS",GATES,"always-active exact gates.")
