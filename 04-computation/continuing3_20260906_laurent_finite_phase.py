"""Exact uniform first-phase and improved tail certificates for the two-anchor model.

Standard library only. All source carries are reconstructed. The finite Bernstein
bank is a polynomial identity certificate on a whole box, not a sample census.
"""
from fractions import Fraction as F
from itertools import combinations, product
from math import comb
from pathlib import Path
import argparse
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(label)


def tidy(p):return {e:F(v) for e,v in p.items() if v}
def const(c,n):return tidy({(0,)*n:F(c)})
def var(i,n):return {tuple(int(j==i) for j in range(n)):F(1)}
def scale(p,c):return tidy({e:v*c for e,v in p.items()})


def add(*rows):
    out={}
    for row in rows:
        for e,v in row.items():out[e]=out.get(e,F(0))+v
    return tidy(out)


def mul(a,b):
    out={}
    for e,v in a.items():
        for f,w in b.items():
            key=tuple(x+y for x,y in zip(e,f))
            out[key]=out.get(key,F(0))+v*w
    return tidy(out)


def power(a,k,n):
    out=const(1,n)
    for _ in range(k):out=mul(out,a)
    return out


def evalp(p,point):
    return sum(v*prod(point[i]**e[i] for i in range(len(e))) for e,v in p.items())


def prod(items):
    out=F(1)
    for v in items:out*=v
    return out


def conv(a,b,n):
    out={}
    for i,v in a.items():
        for j,w in b.items():out[i+j]=add(out.get(i+j,{}),mul(v,w))
    return out


def numeric_conv(a,b):
    out=[F(0)]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):out[i+j]+=v*w
    return trim(out)


def trim(a):
    a=list(map(F,a))
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def divrem(a,b):
    a=trim(a);b=trim(b);q=[F(0)]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b) and any(a):
        offset=len(a)-len(b);c=a[-1]/b[-1];q[offset]=c
        for j,v in enumerate(b):a[offset+j]-=c*v
        a=trim(a)
    return trim(q),a


def horner(a,x):
    out=F(0)
    for c in reversed(a):out=out*x+c
    return out


def variations(values):
    signs=[1 if v>0 else -1 for v in values if v]
    return sum(a!=b for a,b in zip(signs,signs[1:]))


def sturm_positive(a):
    chain=[trim(a),trim([i*a[i] for i in range(1,len(a))])]
    while len(chain[-1])>1:
        r=divrem(chain[-2],chain[-1])[1]
        if not any(r):break
        chain.append([-v for v in r])
    return variations([p[0] for p in chain])-variations([p[-1] for p in chain])


def source():
    x,y,z=(var(i,3) for i in range(3))
    B={j-1:v for j,v in enumerate([const(1,3),const(13,3),const(55,3),x,y,z])}
    C={j-1:v for j,v in enumerate([const(1,3),const(12,3),const(45,3),scale(x,F(2,3)),scale(y,F(3,7))])}
    D={j-1:v for j,v in enumerate([const(1,3),const(11,3),const(36,3),scale(x,F(5,12)),scale(y,F(1,7))])}
    B2=conv(B,B,3);CD=conv(C,D,3)
    grouped={j:add(B2.get(j,{}),scale(CD.get(j-1,{}),2)) for j in range(-2,9)}
    O={j:comb(14,2*j+1) for j in range(7)}
    E={j:comb(14,2*j) for j in range(8)}
    alpha={}
    for a,b,shift in [(O,O,0),(E,E,-1)]:
        for i,v in a.items():
            for j,w in b.items():alpha[i+j+shift]=alpha.get(i+j+shift,0)+v*w
    for j in range(-1,14):need(alpha[j]==comb(28,2*j+2),'independent even/odd square multiplier')
    need(alpha.get(-2,0)==0,'only unavailable lower index removed')
    q={j:scale(grouped[j],alpha[j]) for j in range(-1,9)}
    need(q[-1]==const(28,3),'genuine lower carry retained')
    T={}
    for j,row in q.items():
        for e,v in row.items():T[e+(j+1,)]=v*(-1)**(j%2)
    # Substitute e5=(12/7)e4*u-e3*u^2+10u^3-u^4/11 and Q/s^7=T/s^8.
    f={(0,1,1):F(12,7),(1,0,2):F(-1),(0,0,3):F(10),(0,0,4):F(-1,11)}
    R={}
    for (a,b,c,j),v in T.items():
        R=add(R,mul({(a,b,8-j):v},power(f,c,3)))
    expected={
      (0,2,0):F(-26075790,7),
      (1,1,1):F(153780300,7),(0,2,1):F(-179344800,7),
      (2,0,2):F(-16900975),(1,1,2):F(647843760,7),(0,1,2):F(-1329865290),
      (2,0,3):F(-53986980),(1,0,3):F(1122025905),(0,1,3):F(-4282905900),
      (1,0,4):F(3467704710),(0,1,4):F(5521932000,11),(0,0,4):F(-10070260200),
      (1,0,5):F(-3690469830,11),(0,1,5):F(-9902880),(0,0,5):F(-30313505040),
      (1,0,6):F(6175260),(0,0,6):F(3654364350),
      (0,0,7):F(-1022439600,11),(0,0,8):F(565082)}
    need(R==expected,'complete formal original-root elimination identity')
    return q,T,R


def elementary(rows,k,n):
    return add(*[prodpoly([rows[i] for i in c],n) for c in combinations(range(len(rows)),k)])


def prodpoly(rows,n):
    out=const(1,n)
    for row in rows:out=mul(out,row)
    return out


def anchor_bounds():
    a=[var(i,5) for i in range(5)]
    e1,e2,e3,e4=[elementary(a,k,5) for k in range(1,5)]
    p2=add(*[power(v,2,5) for v in a]);p3=add(*[power(v,3,5) for v in a])
    need(add(power(e1,3,5),scale(mul(e1,p2),-3),scale(p3,2))==scale(e3,6),'third-power anchor identity')
    S={};W={}
    for i,j in combinations(range(5),2):
        rest=[a[k] for k in range(5) if k not in (i,j)]
        S=add(S,mul(power(add(a[i],scale(a[j],-1)),2,5),elementary(rest,2,5)))
    for i,j,k,l in combinations(range(5),4):
        for p,q,r,t in [(i,j,k,l),(i,k,j,l),(i,l,j,k)]:
            W=add(W,power(add(mul(a[p],a[q]),scale(mul(a[r],a[t]),-1)),2,5))
    need(add(power(e2,2,5),scale(e4,-20))==add(S,scale(W,F(1,3))),'explicit fourth-coefficient sum of squares')
    positive={}
    for i in range(5):
        positive=add(positive,mul(power(a[i],2,5),elementary([a[j] for j in range(5) if j!=i],2,5)))
    need(add(mul(e1,e3),scale(e4,-4))==positive,'4e4<=e1e3 coefficient identity')
    need(5*F(71,10)**2-26*F(71,10)-67==F(9,20),'largest-root rational bound')
    need((F(71,10)*59-52)/3==F(1223,10)<123<130,'third-coefficient upper bound')
    need(F(55**2,20)==F(605,4)<152,'fourth-coefficient upper bound')
    need(F(55,10)**5<F(72)**2,'pair-product AM-GM fifth-coefficient bound')
    need(13**2>2*59,'at least three roots are positive')
    need(2002-F(11154,110)>0,'strict positivity of first equation up to one-over-110')
    upper=182-F(20020,90)+F(2002*130,90**2)+F(2002*72,90**4)
    slope=-20020+F(4004*130,90)+F(8008*72,90**3)
    need(upper<0,'negative first equation at one-over-90')
    need(slope<0,'unique simple first phase below one-over-90')


def bernstein(T):
    bounds=[F(130),F(152),F(72)];lo=F(1,120);width=F(1,80)-lo
    cube={}
    for (a,b,c,j),v in T.items():
        factor=v*prod(bounds[i]**e for i,e in enumerate((a,b,c)))
        for k in range(j+1):
            key=(a,b,c,k)
            cube[key]=cube.get(key,F(0))+factor*comb(j,k)*lo**(j-k)*width**k
    cube=tidy(cube);degrees=(2,2,2,9)
    need(tuple(max(e[j] for e in cube) for j in range(4))==degrees,'complete multidegree of box polynomial')
    bank=[]
    for idx in product(*[range(d+1) for d in degrees]):
        val=sum(v*prod(F(comb(i,k),comb(d,k)) for i,k,d in zip(idx,e,degrees))
                for e,v in cube.items() if all(k<=i for k,i in zip(e,idx)))
        need(val<-400,'uniform strict Bernstein coefficient margin')
        bank.append((idx,val))
    need(len(bank)==270,'complete tensor Bernstein universe')
    need(max(v for i,v in bank)==F(-22645374245632441,52254720000000),'published exact largest coefficient')
    # A different direction checks the change of basis by coefficient expansion.
    rebuilt={}
    for idx,value in bank:
        factors=[]
        for i,d in zip(idx,degrees):
            factors.append([(k,comb(d,i)*comb(d-i,k-i)*(-1)**(k-i)) for k in range(i,d+1)])
        for terms in product(*factors):
            exponent=tuple(term[0] for term in terms)
            rebuilt[exponent]=rebuilt.get(exponent,F(0))+value*prod(term[1] for term in terms)
    need(tidy(rebuilt)==cube,'independent power-basis reconstruction of entire Bernstein identity')
    return bank


def tail(R):
    c0=F(26075790,7)
    envelope={
      1:F(12300*153780300,7),2:F(12300*647843760,7),
      3:F(1230000*1122025905),
      4:F(1230000*3467704710)+F(100*5521932000,11),
      6:F(1230000*6175260+10000*3654364350),8:F(10000*565082)}
    # Verify that precisely every positive monomial of the formal remainder
    # contributes to the envelope; the sole retained negative term is -c0 e4^2.
    reconstructed={}
    for (a,b,k),v in R.items():
        if v<=0:continue
        need(a<=1 and b<=1,'positive monomial lies in stated reciprocal-floor tariff')
        coefficient=v*F(123)**a*F(100)**(2-b)
        reconstructed[k]=reconstructed.get(k,F(0))+coefficient
    need(reconstructed==envelope,'no positive eliminated monomial omitted from tail')
    need(all(v>0 and k>0 for k,v in envelope.items()),'positive envelope decreases with phase')
    at=-c0+sum(v/F(75000)**k for k,v in envelope.items())
    need(at<-120000,'entire original-root tail starts by 75000')
    return envelope,at


def controls(q,T,R):
    roots=[F(v,5) for v in (1,3,9,22,30)]
    es=[sum(prod(c) for c in combinations(roots,k)) for k in range(1,6)]
    need(es[:2]==[13,55],'strict model control has both exact anchors')
    e3,e4,e5=es[2:]
    C=[3*e4/7,-2*e3/3,45,-12,1];D=[e4/7,-5*e3/12,36,-11,1]
    for i,a in enumerate(roots):
        need((-1)**i*horner(C,a)>0 and (-1)**i*horner(D,a)>0,'both strict interlacings in rational positive control')
    PP=[182,-20020,2002*e3,-3432*e4,2002*e5]
    need(horner(PP,F(1,110))>0>horner(PP,F(1,90)),'literal first-phase control')
    need(evalp(T,(e3,e4,e5,F(75000)))>0 and horner(PP,F(75000))!=0,
         'tail root-equation omission hostile within the valid shape class')
    need(evalp(T,(F(130),F(0),F(0),F(1,60)))==F(34734566083,93312000)>0,
         'coefficient-box interval cannot be extended blindly')
    hx,hy,hz,hs=F(104),F(50),F(37435088,3898125),F(15,2)
    first=[182,-20020,2002*hx,-3432*hy,2002*hz]
    need(horner(first,hs)==0,'exact original first root in Newton-only hostile')
    value=evalp(T,(hx,hy,hz,hs))/hs
    need(value==F(78541969368658673,18480)>0,'positive full carried response at hostile first root')
    margins=[169-F(5,2)*55,55**2-2*13*hx,hx**2-2*55*hy,hy**2-F(5,2)*hx*hz]
    need(margins==[F(63,2),F(321),F(5316),F(2437924,779625)] and min(margins)>0,
         'all four explicit Newton inequalities hold strictly')
    need(0<hx<130 and 0<hy<152 and 0<hz<72,'hostile is inside the same coefficient box')
    intervals=[(F(1,99),F(1,98)),(F(3,32),F(5,53)),(F(70,53),F(37,28))]
    for a,b in intervals:need(horner(first,a)*horner(first,b)<0,'three other positive simple roots exhaust the quartic')
    beta=[-hz,hy,-hx,55,-13,1]
    need(sturm_positive(beta)==1,'hostile beta polynomial has only one positive real root')
    need(all(v*(-1)**i<0 for i,v in enumerate(beta)),'hostile beta has no negative real root')
    for point in [(F(84),F(35),F(1)),(hx,hy,hz),(e3,e4,e5)]:
        for phase in [F(1,10),F(1,2),F(2)]:
            raw=sum(evalp(row,point)*(-phase)**j for j,row in q.items())
            need(phase*raw==evalp(T,point+(phase,)),'ordinary Laurent evaluation matches full polynomial carrier')
    return value


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--certificate');args=parser.parse_args()
    q,T,R=source();anchor_bounds();bank=bernstein(T);envelope,tailvalue=tail(R);hostile=controls(q,T,R)
    print('UNIFORM first phase: exactly one root in (1/110,1/90), no earlier positive root; no interlacer needed')
    print('BOX 0<=e3<=130, 0<=e4<=152, 0<=e5<=72, 1/120<=s<=1/80: s Q(-s)<-400')
    print('BERNSTEIN count=270 multidegree=(2,2,2,9) maximum=',max(v for i,v in bank))
    print('TAIL original P(-s)=0 and e4>1/100: Q(-s)/(s^7 e4^2)<-120000 for s>=75000')
    print('TAIL_ENVELOPE_AT_75000',tailvalue)
    print('NEWTON_HOSTILE e3=104 e4=50 e5=37435088/3898125 s=15/2 Q(-s)=',hostile)
    print('HOSTILE original P has four simple positive roots; beta B has only one real root; not a model counterexample')
    print('REMAINING any unresolved original phase lies in (1/80,75000); full model sign remains OPEN')
    print('PASS',GATES,'always-active exact gates; no numerical optimization or external algebra package')
    if args.certificate:
        def pack(p):return [[list(e),str(v)] for e,v in sorted(p.items())]
        obj={'full_T':pack(T),'eliminated_Q_over_s7':pack(R),'bernstein':[[list(i),str(v)] for i,v in bank],
             'tail_envelope':[[k,str(v)] for k,v in sorted(envelope.items())],'tail_at_75000':str(tailvalue)}
        Path(args.certificate).write_bytes((json.dumps(obj,sort_keys=True,indent=2)+'\n').encode())


if __name__=='__main__':main()
