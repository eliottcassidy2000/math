"""Exact two-residue moment floor and original-root phase tail.

Standard-library symbolic polynomials, full series division and carried-row
reconstruction; no producer imports and no numerical search.
"""
from fractions import Fraction as F
from itertools import combinations,permutations
from math import comb
from pathlib import Path
import json,sys,argparse
sys.stdout.reconfigure(newline='\n')
GATES=0


def need(ok,why):
    global GATES
    GATES+=1
    if not ok:raise ArithmeticError(why)


def clean(p):return {e:F(v) for e,v in p.items() if v}
def constant(v,n=3):return clean({(0,)*n:F(v)})
def variable(i,n=3):return {tuple(int(j==i) for j in range(n)):F(1)}
def scale(p,c):return clean({e:v*c for e,v in p.items()})


def add(*ps):
    out={}
    for p in ps:
        for e,v in p.items():out[e]=out.get(e,F(0))+v
    return clean(out)


def multiply(a,b):
    out={}
    for e,v in a.items():
        for f,w in b.items():
            k=tuple(x+y for x,y in zip(e,f));out[k]=out.get(k,F(0))+v*w
    return clean(out)


def power(p,k,n=3):
    out=constant(1,n)
    for _ in range(k):out=multiply(out,p)
    return out


def product(ps,n):
    out=constant(1,n)
    for p in ps:out=multiply(out,p)
    return out


def evaluate(p,a):
    out=F(0)
    for e,v in p.items():
        for x,k in zip(a,e):v*=x**k
        out+=v
    return out


def determinant(matrix):
    n=len(matrix);out={}
    for order in permutations(range(n)):
        sign=(-1)**sum(order[i]>order[j] for i in range(n) for j in range(i+1,n))
        out=add(out,scale(product([matrix[i][order[i]] for i in range(n)],3),sign))
    return out


def convolution(a,b):
    out={}
    for i,v in a.items():
        for j,w in b.items():out[i+j]=add(out.get(i+j,{}),multiply(v,w))
    return out


def moments():
    x,y,z=(variable(i) for i in range(3));one=constant(1)
    denominator=[one,constant(-13),constant(55),scale(x,-1),y,scale(z,-1)]
    inverse=[one]
    for k in range(1,6):
        inverse.append(scale(add(*[multiply(denominator[j],inverse[k-j]) for j in range(1,k+1)]),-1))
    cn=[one,constant(-12),constant(45),scale(x,F(-2,3)),scale(y,F(3,7))]
    dn=[one,constant(-11),constant(36),scale(x,F(-5,12)),scale(y,F(1,7))]
    rows=[]
    for numerator in [cn,dn]:
        row=[add(*[multiply(numerator[j],inverse[k-j]) for j in range(min(k,4)+1)]) for k in range(6)]
        for k in range(6):
            rebuilt=add(*[multiply(denominator[j],row[k-j]) for j in range(k+1)])
            need(rebuilt==(numerator[k] if k<5 else {}),'full rational-series division through fifth moment')
        rows.append(row)
    C,D=rows
    expectedC=[one,one,constant(3),add(scale(x,F(1,3)),constant(-16)),
       add(scale(x,F(16,3)),constant(-373),scale(y,F(-4,7))),
       add(scale(x,54),scale(y,F(-59,7)),z,constant(-3969))]
    expectedD=[one,constant(2),constant(7),add(scale(x,F(7,12)),constant(-19)),
       add(scale(x,F(115,12)),constant(-632),scale(y,F(-6,7))),
       add(scale(x,F(199,2)),scale(y,F(-92,7)),z,constant(-7171))]
    need(C==expectedC and D==expectedD,'all published C and D moments')
    delta=add(x,constant(-75))
    Cshift2=determinant([[C[i+j+1] for j in range(2)] for i in range(2)])
    need(Cshift2==scale(delta,F(1,3)),'C shifted second determinant')
    Cshift3=determinant([[C[i+j+1] for j in range(3)] for i in range(3)])
    wantedC=add(scale(multiply(x,power(delta,2)),F(-1,27)),scale(multiply(delta,y),F(15,7)),
                scale(power(y,2),F(-16,49)),scale(multiply(delta,z),F(1,3)))
    need(Cshift3==wantedC,'exact C shifted third Gram determinant')
    Dord3=determinant([[D[i+j] for j in range(3)] for i in range(3)])
    wantedD=add(scale(add(constant(-333),scale(delta,2334),scale(power(delta,2),-49)),F(1,144)),scale(y,F(-18,7)))
    need(Dord3==wantedD,'exact D ordinary third Gram determinant')
    return C,D,Cshift2,Cshift3,Dord3


def root_and_floor_arithmetic():
    roots=[variable(i,6) for i in range(5)];M=variable(5,6)
    es=[add(*[product([roots[i] for i in c],6) for c in combinations(range(5),k)]) for k in range(1,6)]
    lhs=add(multiply(M,es[3]),scale(es[4],-5))
    rhs=add(*[multiply(add(M,scale(roots[i],-1)),product([roots[j] for j in range(5) if j!=i],6)) for i in range(5)])
    need(lhs==rhs,'five-node bounded-support identity M e4 minus5 e5')
    p2=add(*[power(a,2,6) for a in roots]);p3=add(*[power(a,3,6) for a in roots])
    need(add(power(es[0],3,6),scale(multiply(es[0],p2),-3),scale(p3,2))==scale(es[2],6),'third-power coefficient identity')
    need(5*F(71,10)**2-26*F(71,10)-67==F(9,20)>0,'strict largest-beta-root bound')
    need((F(71,10)*59-52)/3==F(1223,10)<123,'third coefficient bound for tail')
    K=F(15,7)+F(71,150)
    c=F(75)/(27*K)
    need(K==F(2747,1050) and c==F(8750,8241),'exact shifted-Gram domination slope')
    denominator=F(2334)/c-F(2592,7)
    need(denominator>0,'final delta elimination keeps inequality orientation')
    floor=F(333)/denominator
    need(floor==F(161875,888583),'published strict fourth-coefficient floor')
    need(floor-F(9,50)==F(96503,44429150)>0,'simple strict floor nine-fiftieths')
    return floor


def original_tail():
    x,y,z=(variable(i) for i in range(3))
    GB=[constant(1),constant(13),constant(55),x,y,z]
    GC=[constant(1),constant(12),constant(45),scale(x,F(2,3)),scale(y,F(3,7))]
    GD=[constant(1),constant(11),constant(36),scale(x,F(5,12)),scale(y,F(1,7))]
    beta={i-1:v for i,v in enumerate(GB)};crow={i-1:v for i,v in enumerate(GC)};drow={i-1:v for i,v in enumerate(GD)}
    beta2=convolution(beta,beta);cd=convolution(crow,drow)
    O={j:comb(14,2*j+1) for j in range(7)};E={j:comb(14,2*j) for j in range(8)}
    alpha={}
    for p,q,shift in [(O,O,0),(E,E,-1)]:
        for i,v in p.items():
            for j,w in q.items():alpha[i+j+shift]=alpha.get(i+j+shift,0)+v*w
    for j in range(-1,14):need(alpha[j]==comb(28,2*j+2),'full independent alpha parity convolution')
    raw={j:scale(add(beta2.get(j,{}),scale(cd.get(j-1,{}),2)),alpha[j]) for j in range(-1,9)}
    need(raw[-1]==constant(28),'genuine lower carry unchanged')
    P={j:scale(beta[j],O[j]) for j in range(5)}
    need(P=={0:constant(182),1:constant(20020),2:scale(x,2002),3:scale(y,3432),4:scale(z,2002)},'same original first-row normalization')
    f={(0,1,1):F(12,7),(1,0,2):F(-1),(0,0,3):F(10),(0,0,4):F(-1,11)}
    R={}
    for j,row in raw.items():
        for (a,b,c),v in row.items():
            R=add(R,multiply({(a,b,7-j):v*(-1)**(j%2)},power(f,c)))
    expected={
      (0,2,0):F(-26075790,7),
      (1,1,1):F(153780300,7),(0,2,1):F(-179344800,7),
      (2,0,2):F(-16900975),(1,1,2):F(647843760,7),(0,1,2):F(-1329865290),
      (2,0,3):F(-53986980),(1,0,3):F(1122025905),(0,1,3):F(-4282905900),
      (1,0,4):F(3467704710),(0,1,4):F(5521932000,11),(0,0,4):F(-10070260200),
      (1,0,5):F(-3690469830,11),(0,1,5):F(-9902880),(0,0,5):F(-30313505040),
      (1,0,6):F(6175260),(0,0,6):F(3654364350),(0,0,7):F(-1022439600,11),(0,0,8):F(565082)}
    need(R==expected,'unchanged exact original-root response after eliminating e5')
    inverse_floor=F(50,9);envelope={}
    for (a,b,k),v in R.items():
        if v<=0:continue
        need(a<=1 and b<=1 and k>0,'every positive monomial has the stated coefficient tariff')
        envelope[k]=envelope.get(k,F(0))+v*F(123)**a*inverse_floor**(2-b)
    at=F(-26075790,7)+sum(v/F(4100)**k for k,v in envelope.items())
    need(at==F(-5961609014706655321683472522194343,99603957308055354000000000000)<-59000,'new entire tail cutoff4100')
    need(all(v>0 for v in envelope.values()),'tail envelope monotonic in inverse phase')
    ratio1=F(75)*inverse_floor+F(8241,8750)
    ratio2=ratio1*inverse_floor
    need(ratio1==F(10962223,26250) and ratio2==F(10962223,4725),'coupled x-over-y and x-over-y-squared bounds from strict slope')
    tariff={(1,1):ratio1,(1,0):ratio2,(0,1):inverse_floor,(0,0):inverse_floor**2}
    coupled={}
    for (a,b,k),v in R.items():
        if v>0:coupled[k]=coupled.get(k,F(0))+v*tariff[(a,b)]
    at2500=F(-26075790,7)+sum(v/F(2500)**k for k,v in coupled.items())
    need(all(v>0 for v in coupled.values()),'coupled envelope monotonic in inverse phase')
    need(at2500==F(-653528391359305169367452997041401,13323669433593750000000000000)<-49000,'coupled entire tail cutoff2500')
    return raw,R,envelope,at,coupled,at2500


def controls(C,D,C2,C3,D3,raw):
    hostile=(F(104),F(50),F(37435088,3898125))
    Co=determinant([[C[i+j] for j in range(3)] for i in range(3)])
    Ds=determinant([[D[i+j+1] for j in range(2)] for i in range(2)])
    need([evaluate(p,hostile) for p in (Co,C2,D3,Ds)]==[F(2693,63),F(29,3),F(3338,63),F(103,3)],'old hostile passes both lower Hankel packets')
    Cvi=evaluate(C[4],hostile)-F(71,10)*evaluate(C[3],hostile)
    Dvi=evaluate(D[4],hostile)-F(71,10)*evaluate(D[3],hostile)
    need(Cvi==F(2159,105)>0 and Dvi==F(1091,42)>0,'bounded-support endpoint rejects both order-four packets')
    need(evaluate(C3,hostile)==F(-70052935886,81860625)<0,'order-five C shifted Gram also rejects hostile')
    phase=F(15,2);target=sum(evaluate(v,hostile)*(-phase)**j for j,v in raw.items())
    need(target==F(78541969368658673,18480)>0,'hostile full carried response preserved')
    x,y,z=hostile
    need(182-20020*phase+2002*x*phase**2-3432*y*phase**3+2002*z*phase**4==0,'hostile original first-root condition preserved')
    a=[F(v,5) for v in (1,3,9,22,30)]
    es=[sum(product_numbers(c) for c in combinations(a,k)) for k in range(1,6)]
    need(es[:2]==[13,55],'strict model control preserves anchors')
    coeff=tuple(es[2:]);x,y,z=coeff
    cv=[3*y/7,-2*x/3,45,-12,1];dv=[y/7,-5*x/12,36,-11,1]
    def at(poly,t):return sum(v*t**j for j,v in enumerate(poly))
    for numerator,m in [(cv,C),(dv,D)]:
        weights=[]
        for i,t in enumerate(a):
            weight=at(numerator,t)/product_numbers(t-b for j,b in enumerate(a) if j!=i)
            need(weight>0,'literal positive residue at every strict beta root')
            weights.append(weight)
        for j in range(6):need(sum(w*t**j for w,t in zip(weights,a))==evaluate(m[j],coeff),'direct residue moments match polynomial division')
    need(y>F(161875,888583),'strict positive model satisfies new floor')
    need(5*z<=max(a)*y,'literal bounded-support coefficient inequality')
    boundary=(F(75),F(0),F(0))
    for j in range(6):need(evaluate(C[j],boundary)==(F(1) if j==0 else F(3)**j/3),'weak C-only repeated-root measure at0 and3')
    need(evaluate(C3,boundary)==0 and evaluate(D3,boundary)==F(-37,16),'D Gram excludes exactly the old weak C boundary')
    need(at([0,F(-5,12)*75,36,-11,1],5)==F(-25,4),'literal D fails repeated beta root5')
    need(sum(evaluate(v,coeff)*F(-4100)**j for j,v in raw.items())>0
         and 182-20020*4100+2002*x*4100**2-3432*y*4100**3+2002*z*4100**4!=0,
         'large-phase root-condition omission remains a hostile')


def product_numbers(a):
    out=F(1)
    for v in a:out*=v
    return out


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--certificate');args=parser.parse_args()
    C,D,C2,C3,D3=moments();floor=root_and_floor_arithmetic();raw,R,envelope,tail,coupled,tail2500=original_tail();controls(C,D,C2,C3,D3,raw)
    print('MOMENTS full C/B and D/B division through order5, ordinary and shifted determinant identities exact')
    print('FLOOR e4>',floor,'>9/50; uses both interlacers, actual nonnegative beta roots, and both anchors')
    print('DELTA y>(8750/8241)(e3-75); repeated C-only boundary rejected by D Gram=-37/16')
    print('HOSTILE lower Hankel determinants PASS; order4 upper-support violations C=2159/105 D=1091/42')
    print('CONTROL independent-cap original-root envelope: s>=4100 implies Q(-s)/(s^7 e4^2)<-59000')
    print('TAIL_ENVELOPE_AT_4100',tail)
    print('COUPLED_RATIOS x/y<10962223/26250; x/y^2<10962223/4725')
    print('TAIL same original-root full response: s>=2500 implies Q(-s)/(s^7 e4^2)<-49000')
    print('COUPLED_TAIL_ENVELOPE_AT_2500',tail2500)
    print('SCOPE candidate until independent audit; remaining finite phases and general Laurent noncancellation OPEN')
    print('PASS',GATES,'always-active exact gates; no numerical search or producer imports')
    if args.certificate:
        def pack(p):return [[list(e),str(v)] for e,v in sorted(p.items())]
        bank={'C_moments':[pack(p) for p in C],'D_moments':[pack(p) for p in D],
              'C_shifted3':pack(C3),'D_ordinary3':pack(D3),'original_root_response':pack(R),
              'tail_envelope':[[k,str(v)] for k,v in sorted(envelope.items())],
              'coupled_tail_envelope':[[k,str(v)] for k,v in sorted(coupled.items())],
              'coupled_ratios':['10962223/26250','10962223/4725']}
        Path(args.certificate).write_bytes((json.dumps(bank,sort_keys=True,indent=2)+'\n').encode())


if __name__=='__main__':main()
