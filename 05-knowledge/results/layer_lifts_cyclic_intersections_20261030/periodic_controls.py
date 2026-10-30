#!/usr/bin/env python3
"""Independent integer coefficient checks for primitive-controller phase bounds."""
from itertools import product
from math import gcd
from pathlib import Path
import argparse,json

checks=0
def need(ok,msg):
    global checks
    checks+=1
    if not ok:raise RuntimeError(msg)
def coeff(w):
    P=Q=1;B=0
    for a in w:P,Q,B=3*P,2**a*Q,3*B+Q
    return P,Q,B
def content(w):
    P,Q,B=coeff(w);return gcd(abs(P-Q),B)
def primitive(w):return all(w!=w[:j]*(len(w)//j) for j in range(1,len(w)) if len(w)%j==0)
def vp(x,p):
    need(x!=0,'finite valuation');r=0
    while x%p==0:x//=p;r+=1
    return r
def main():
    words=phases=0;examples=[]
    for d in range(1,6):
        for w in product(range(1,5),repeat=d):
            if not primitive(w):continue
            words+=1;P,Q,B=coeff(w);g=gcd(abs(P-Q),B)
            constants=[]
            for j in range(1,d):
                v=w[:j];p,q,b=coeff(v);C=b*(P-Q)-B*(p-q)
                need(C!=0,'noncommuting proper primitive prefix');constants.append(C);phases+=1
                for k in range(21):need(C%content(w*k+v)==0,'full-integer phase divisor')
            for k in range(1,21):
                S=(P**k-Q**k)//(P-Q)
                need(content(w*k)==g*S,'exact period boundary content')
                for ell in (5,7,11):
                    if (P-Q)%ell==0:
                        expected=vp(g,ell)+vp(k,ell)
                    else:
                        ratio=P*pow(Q,-1,ell)%ell;t=1;x=ratio
                        while x!=1:t+=1;x=x*ratio%ell
                        expected=0 if k%t else vp(P**t-Q**t,ell)+vp(k//t,ell)
                    need(vp(g*S,ell)==expected,'independent boundary valuation law')
            if w in ((1,2),(1,3),(1,2,3),(2,3,1)):
                examples.append(dict(word=w,commutator_constants=constants,base_content=g))
    marked_cases=0;resonant=0
    entries=[()]+[w for d in range(1,4) for w in product(range(1,4),repeat=d)]
    controllers=[w for d in range(1,4) for w in product(range(1,4),repeat=d)]
    for u in entries:
        for w in controllers:
            P,Q,B=coeff(w)
            for j in range(1,len(w)):
                v=w[:j];p,q,b=coeff(v+u);C=b*(P-Q)-B*(p-q)
                marked_cases+=1
                if C==0:
                    resonant+=1
                    need((v+u)+w==w+(v+u),'zero commutator is word resonance')
                else:
                    for k in range(6):need(C%content(u+w*k+v)==0,'entry-preserving phase divisor')
    # The excluded nonprimitive boundary really has a commuting proper prefix.
    P,Q,B=coeff((1,1));p,q,b=coeff((1,));need(b*(P-Q)-B*(p-q)==0,'nonprimitive hostile')
    return dict(checks=checks,primitive_words=words,proper_phases=phases,marked_entry_cases=marked_cases,
                resonant_marked_entries=resonant,alphabet=[1,2,3,4],max_length=5,
                repetitions_checked=list(range(21)),valuation_primes=[5,7,11],examples=examples,
                scope='full-integer theorem is proved in COLLATZ.md; these are bounded exact controls')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).with_name('periodic.json'));a=p.parse_args()
    report=json.dumps(main(),indent=2,sort_keys=True)+'\n';a.output.write_text(report);print(report,end='')
