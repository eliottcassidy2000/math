"""Exact fixed-core inverse; no unbounded odd-perfect-number claim."""
from pathlib import Path
from math import gcd,isqrt
from fractions import Fraction
import json
from sympy import factorint,isprime,divisor_sigma,primerange
HERE=Path(__file__).resolve().parent
def need(p,msg):
    if not p:raise ValueError(msg)
def sigma(n):return int(divisor_sigma(n))
def completion(C,odd=True):
    A=sigma(C);D=2*C-A
    out=dict(C=C,A=A,D=D)
    if D<=0:return dict(out,rejection='nonpositive deficiency')
    q=(A+D-1)//D;B=A-D*(q-1)
    out.update(q=q,B=B)
    if q<2:return dict(out,rejection='q below two')
    if not isprime(q):return dict(out,rejection='composite forced q',factors={str(p):int(e) for p,e in factorint(q).items()})
    if gcd(C,q)!=1:return dict(out,rejection='support collision')
    if odd and (C%2==0 or isqrt(C)**2!=C or q%4!=1):return dict(out,rejection='Euler shape')
    if A%B:return dict(out,rejection='nonintegral power ratio')
    z=A//B;a=0
    while z%q==0:z//=q;a+=1
    out.update(exponent=a,residual=z)
    if z!=1 or a<1:return dict(out,rejection='not a positive q power')
    if odd and a%4!=1:return dict(out,rejection='Euler exponent')
    N=C*q**a;need(sigma(N)==2*N,'independent sigma replay')
    return dict(out,N=N,status='realized perfect completion')
def main():
    # The inverse itself also applies to arbitrary positive cores when Euler
    # shape is disabled. Exhaustion here is an independent bounded control.
    checked=0;hits=[]
    for C in range(1,2001):
        z=completion(C,False)
        expected=[]
        for q in primerange(2,50):
            if gcd(C,q)>1:continue
            for a in range(1,8):
                if sigma(C)*sum(q**i for i in range(a+1))==2*C*q**a:expected.append((q,a))
        actual=[] if 'N' not in z else [(z['q'],z['exponent'])]
        need(expected==[p for p in actual if p[0]<50 and p[1]<=7],'independent finite universe')
        checked+=1
        if actual:hits.append(z)
    base=3**2*7**2*11**2;A=sigma(base);D=2*base-A
    # For q>=17, geometric domination gives the rigorous p<=53 cutoff.
    L=32*base-17*A;need(L>0 and L*54**2>17*A*55,'uniform bound beyond53')
    need(13*D>A and A>12*D,'forced q>=13')
    # At q=13, alpha>=5 already exceeds the limiting required abundancy.
    need(Fraction(sum(13**i for i in range(6)),13**5)>Fraction(2*base,A),'exclude all alpha>=5 at q13')
    # alpha=1 gives K p^2=M(p+1), with nonsquare discriminant.
    K=26*base-14*A;M=14*A;disc=M*M+4*K*M
    need(K>0 and isqrt(disc)**2!=disc,'q13 alpha1 has no integer p')
    rows=[dict(p=p,**completion(base*p*p)) for p in primerange(3,54) if base%p]
    need(all('N' not in z for z in rows),'complete fixed-core family excluded')
    C=base*13**2;spoof=completion(C);need(spoof['q']==22021 and spoof['rejection']=='composite forced q','Descartes hostile')
    need(sigma(C)*(22021+1)==2*C*22021,'formal balance')
    excess=sigma(C*22021)-2*C*22021
    need(excess==1600*sigma(C),'actual nonprime correction')
    report=dict(status='PROVED_INVERSE_WITH_FINITE_EXACT_CONTROLS',bounded_control_cores=checked,bounded_realizations=hits,
        fixed_base=base,fixed_base_ratio=str(Fraction(A,D)),q13_discriminant=disc,q_ge17_p_bound=53,
        exhaustive_core_extensions=rows,descartes=dict(completion=spoof,formal_number=C*22021,actual_excess=excess),
        scope='All primes p,q and all Euler exponents for base3^2*7^2*11^2; not a new global odd-perfect bound.')
    (HERE/'perfect_completion.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS fixed-core inverse,2000 bounded controls,all-exponent base family and Descartes hostile')
if __name__=='__main__':main()
