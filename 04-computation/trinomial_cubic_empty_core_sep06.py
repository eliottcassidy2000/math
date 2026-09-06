#!/usr/bin/env python3
"""Exact three-minor certificate for (-21,2g-21,3g-21).

Scope: integers g>=11, gcd(g,21)=1. SymPy computes the full anchored
remainder and symbolic Hermite minors. Named literal controls are separate.
No repository producer or theorem import; explicit gates survive -O.
"""
import hashlib
import json
from math import factorial, gcd

import sympy as S

T,G,Z = S.symbols("tau g s")
CHECKS = 0
# Coefficients in descending powers of s for F(g)=F(s+11).
SHIFTED = {
    "F5": [21956511799903,447248962924222,3643877376899891,
           14842834323259268,30227960317081116,24622377792469920],
    "F11": [140892841323838612672637,5965004072995999656704026,
            114714746974061256056387574,1322756664440858664210176156,
            10161118159501638625772346617,54598896731247155548982886362,
            209397347882224306724476576476,573182745817688552468969130664,
            1097399763351738257252121426176,1399541151141843962715007684992,
            1070005363538327691749959345920,371517273613263339950411923200],
    "F3": [74863,860377,3285600,4167636],
    "F15": [141387288101294140031271809,7376544414686784800082667655,
            179216442364749171385042322678,2689550519731785864607802436098,
            27881104094183316787745471091341,211468678134040753645990344857219,
            1212233732617430261960470627128172,5347770140720395948826740811577460,
            18303547297713040708451344127658688,48600649295038820551605884973566416,
            99288474821888816842823372182321600,153249236436278889420201068735488192,
            172971501480743281449391277285037312,134765740119228340806711738891048960,
            64801772706710330583381523716710400,14495117169229267635561205094400000]
}
DENOMINATORS = [28510570408320000,119489335876142491574207692800000000,
                41207553558341744647619275991171138077065216000000000000]


def require(test,label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def falling(value,degree):
    return S.prod(value-k for k in range(degree))


def from_descending(coefficients,variable):
    return sum(S.Integer(value)*variable**(len(coefficients)-1-i)
               for i,value in enumerate(coefficients))


def channels(a,b,c,mass):
    answer = []
    for x in range(mass+1):
        numerator = a*x-b*(mass-x)
        if numerator % (c-b):
            continue
        z = numerator//(c-b)
        y = mass-x-z
        if min(y,z)>=0:
            answer.append(((x,y,z),factorial(mass)//
                           (factorial(x)*factorial(y)*factorial(z))))
    return sorted(answer,key=lambda item:item[0][2])


def literal_rows(a,b,c,maximum):
    state = {(0,0):1}
    answer = []
    for _ in range(maximum):
        following = {}
        for (charge,z),value in state.items():
            for delta,dz in ((-a,0),(b,0),(c,1)):
                key = charge+delta,z+dz
                following[key] = following.get(key,0)+value
        state = following
        answer.append({z:value for (charge,z),value in state.items() if charge==0})
    return answer


def main():
    u = G-7
    p = S.Poly(72*T**3+504*u*T**2+84*u*(u-1)*T+u*(u-1)*(u-2),T)
    first = S.Poly(sum(falling(G,10-j)*T**j/
                       (factorial(9-3*j)*factorial(1+2*j)) for j in range(4)),T)
    content = falling(G,7)/factorial(9)
    require(S.cancel(first.as_expr()-content*p.as_expr())==0,"complete first normalization")
    positive_factor = falling(2*G,14)/128
    q = S.Poly(sum(falling(2*G,21-j)*T**j/
                   (factorial(21-3*j)*factorial(2*j)) for j in range(8)),T)
    q_normalized = S.Poly(sum(S.cancel(q.nth(j)/positive_factor)*T**j for j in range(8)),T)
    inverse = -(72*T**2+504*u*T+84*u*(u-1))/(u*(u-1)*(u-2))
    require(S.cancel(S.rem(T*inverse-1,p.as_expr(),T))==0,"inverse lower carry")
    raw_remainder = S.Poly(S.rem(q_normalized.as_expr()*inverse,p.as_expr(),T),T)
    r = [S.factor(raw_remainder.nth(j)) for j in range(3)]
    require(raw_remainder.degree()==2,"complete anchored residue degree")
    expected = [
        -(G-8)*(G-7)*(981739413*G**4-29190405089*G**3+323857306851*G**2
                     -1589986674295*G+2916051235200)/28510570408320000,
        -(G-7)*(1715049346*G**4-49273046168*G**3+530071519257*G**2
                -2530936512265*G+4525872382400)/593970216840000,
        -(70405697381*G**4-1950797880338*G**3+20268905917687*G**2
          -93593992728910*G+162061152680400)/4157791517880000]
    require(all(S.cancel(a-b)==0 for a,b in zip(r,expected)),"symbolic normalized remainder")
    # Newton sums from the monic first cubic give the six needed power traces.
    monic = [7*u,S.Rational(7,6)*u*(u-1),u*(u-1)*(u-2)/72]
    traces = [S.Integer(3)]
    for n in range(1,7):
        value = sum(monic[k-1]*traces[n-k] for k in range(1,min(n,3)+1) if n-k>0)
        if n<=3:
            value += n*monic[n-1]
        traces.append(S.expand(-value))
    matrix = S.Matrix(3,3,lambda i,j:S.expand(sum(r[k]*traces[i+j+k] for k in range(3))))
    factors = {name:S.expand(from_descending(values,G-11)) for name,values in SHIFTED.items()}
    expected_minors = [
        (G-7)*factors["F5"]/DENOMINATORS[0],
        (G-8)*(G-7)**2*factors["F11"]/DENOMINATORS[1],
        (G-8)**2*(G-7)**4*factors["F3"]*factors["F15"]/DENOMINATORS[2]]
    for k,expected_minor in enumerate(expected_minors,1):
        actual = (-1)**k*matrix[:k,:k].det(method="domain-ge")
        require(S.cancel(actual-expected_minor)==0,"symbolic signed minor "+str(k))
    for name,values in SHIFTED.items():
        require(all(value>0 for value in values),"all shifted coefficients positive "+name)
    require(S.cancel(S.discriminant(p.as_expr(),T)-15552*(G-8)*(G-7)**2*factors["F3"])==0,
            "positive first discriminant")
    # Actual real-core theorem sector: at most one of three first roots lies here.
    boundary = -S.Rational(4,27)
    require(S.cancel(27*p.diff().eval(boundary)-(2268*(G-11)**2+11844*(G-11)+11216))==0,
            "strict sector derivative at left boundary")
    require(S.cancel(p.diff().diff().as_expr()-(432*T+1008*u))==0,
            "sector derivative increases")
    print("P_normalized",p.as_expr())
    print("positive_factor K=(2g)_14/128; first content=(g)_7/9!")
    print("normalized_anchored_remainder",[str(value) for value in r])
    print("signed_minor_factorization",["(g-7)F5/d1","(g-8)(g-7)^2 F11/d2",
                                         "(g-8)^2(g-7)^4 F3 F15/d3"])
    print("positive_denominators",DENOMINATORS)
    print("positive_shifted_factors",json.dumps(SHIFTED,sort_keys=True))
    manifest = []
    for g in (11,13,16,17,19,20,22,29):
        a,b,c = 21,2*g-21,3*g-21
        require(gcd(g,21)==1 and gcd(a,gcd(b,c))==1,"primitive named control")
        native = literal_rows(a,b,c,2*g)
        require(all(not row for row in native[:g-1]),"every earlier mass empty")
        raw = [channels(a,b,c,mass) for mass in (g,2*g)]
        for mass,row in zip((g,2*g),raw):
            require(native[mass-1]=={v[2]:weight for v,weight in row},"native/arithmetic full row")
        require([v for v,_ in raw[0]]==[(g-10+j,9-3*j,1+2*j) for j in range(4)],"four first channels")
        require([v for v,_ in raw[1]]==[(2*g-21+j,21-3*j,2*j) for j in range(8)],"eight doubled channels")
        pf,qf = [S.Poly(sum(weight*T**j for j,(_,weight) in enumerate(row)),T) for row in raw]
        require(S.cancel(pf.as_expr()-first.as_expr().subs(G,g))==0 and
                S.cancel(qf.as_expr()-q.as_expr().subs(G,g))==0,"literal symbolic row specialization")
        rf = S.Poly(S.rem(qf.as_expr()*S.invert(T,pf.as_expr(),T),pf.as_expr(),T),T)
        scale = positive_factor.subs(G,g)
        require(all(S.cancel(rf.nth(k)/scale-r[k].subs(G,g))==0 for k in range(3)),"literal anchored residue")
        multiplication = S.zeros(3)
        for column in range(3):
            reduced = S.Poly(S.rem(T**(column+1),pf.as_expr(),T),T)
            for row in range(3):
                multiplication[row,column] = reduced.nth(row)
        response = sum((rf.nth(k)*multiplication**k for k in range(3)),S.zeros(3))
        hermite = S.Matrix(3,3,lambda i,j:S.trace(response*multiplication**(i+j)))
        minors = [(-1)**k*hermite[:k,:k].det() for k in range(1,4)]
        require(all(value>0 for value in minors),"three positive signed literal minors")
        require(all(S.cancel(value/scale**k-expected_minors[k-1].subs(G,g))==0
                    for k,value in enumerate(minors,1)),"multiplication/Newton matrix minors")
        require(pf.count_roots(-S.oo,0)==3 and S.gcd(pf,pf.diff()).degree()==0,"three simple negative first roots")
        require(S.gcd(pf,qf).degree()==0,"coprime complete rows")
        sector_count = pf.count_roots(boundary,0)
        require(sector_count==(1 if g<=20 else 0),"named real-core sector count")
        if g==11:
            require(minors==[S.Integer(752350432547692),S.Rational(21236578848830718128306804,3),
                             S.Rational(942083106929885721679784497035400,27)],"inherited cubic Hermite control")
        print("named_control",(-a,b,c),"first_mass",g,"channels",[len(row) for row in raw],
              "signed_minors_positive",True,"real_core_sector_roots",sector_count)
        manifest.append([g,[str(v) for v in minors],int(sector_count)])
    # A nonprimitive member retains the formal mass rows but not first mass g.
    hostile = literal_rows(21,3,15,12)
    hostile_mass = next(i for i,row in enumerate(hostile,1) if row)
    require(hostile_mass==4,"nonprimitive g12 first-return hostile")
    print("nonprimitive_hostile (-21,3,15), g=12, first support return=4")
    print("PROVED: integral g>=11, gcd(g,21)=1; first nonzero mass g or 2g, both attainable")
    print("SCOPE: this fixed cubic-first-row family; general higher-channel signs remain OPEN")
    encoded = json.dumps({"shifted":SHIFTED,"denominators":DENOMINATORS,"controls":manifest},
                         sort_keys=True,separators=(",",":")).encode()
    print("semantic_sha256",hashlib.sha256(encoded).hexdigest())
    print("explicit_checks",CHECKS)


if __name__=="__main__":
    main()
