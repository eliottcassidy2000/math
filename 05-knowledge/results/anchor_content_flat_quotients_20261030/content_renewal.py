#!/usr/bin/env python3
"""Exact word/content identities and finite affine return clocks; no orbit conjecture."""
import argparse
from collections import deque
from heapq import heappop, heappush
from itertools import product
from math import gcd
from fractions import Fraction
from pathlib import Path
import json

checks=0
def need(ok,msg):
    global checks
    checks+=1
    if not ok: raise RuntimeError(msg)

def clock(w):
    P,Q,B=1,1,0
    for a in w:
        P,Q,B=3*P,(1<<a)*Q,3*B+Q
    return P,Q,B

def closed(w):
    d=len(w)
    return 3**d,2**sum(w),sum(3**(d-i-1)*2**sum(w[:i]) for i in range(d))

def content(w):
    P,Q,B=clock(w)
    return gcd(abs(P-Q),B)

def order_affine(q,a):
    inv=pow(2**a,-1,q)
    m,b=1,0
    for t in range(1,q*q+1):
        m,b=3*inv*m%q,(3*b+1)*inv%q
        if (m,b)==(1,0):return t
    raise RuntimeError('order bound')

def shortest(q,alphabet,weighted=False):
    """Dijkstra on actual affine coefficient states, tracking a nonempty return."""
    identity=(1,0)
    dist={identity:0}; word={identity:()}; queue=[(0,identity)]
    best=None; bestword=None; settled=0
    while queue:
        cost,(m,b)=heappop(queue)
        if cost != dist[(m,b)]:continue
        if best is not None and cost>=best:break
        settled+=1
        for a in alphabet:
            inv=pow(2**a,-1,q)
            nxt=(3*inv*m%q,(3*b+1)*inv%q)
            c=cost+(a if weighted else 1)
            if nxt==identity:
                if best is None or c<best:best,bestword=c,word[(m,b)]+(a,)
            elif c<dist.get(nxt,10**100):
                dist[nxt]=c;word[nxt]=word[(m,b)]+(a,);heappush(queue,(c,nxt))
    need(best is not None,'finite group has return')
    need(content(bestword)%q==0,'shortest return word has full content')
    return dict(cost=best,word=list(bestword),states_settled_below_cost=settled)

def native(n,steps):
    w=[]
    for _ in range(steps):
        z=3*n+1;a=(z&-z).bit_length()-1;n=z>>a;w.append(a)
    return tuple(w),n

def compositions(A,d):
    if d==1:
        yield (A,)
    else:
        for a in range(1,A-d+2):
            for tail in compositions(A-a,d-1):yield (a,)+tail

def minimal_cost_certificates():
    # The large alphabet is only a discovery route. Minimality is independently
    # checked from scalar clock congruences and at most 357 literal compositions.
    out=[]
    for q,cap in ((5,4),(25,20),(125,100),(625,40)):
        answer=shortest(q,tuple(range(1,cap+1)),True)
        C=answer['cost'];w=tuple(answer['word'])
        need(cap>=C or cap==4,'alphabet contains every cheaper letter')
        candidates=[];tested=0
        for A in range(1,C):
            for d in range(1,A+1):
                if pow(3,d,q)!=pow(2,A,q):continue
                words=list(compositions(A,d));tested+=len(words)
                need(all(clock(v)[2]%q for v in words),'independent lower-cost exclusion')
                candidates.append(dict(binary_clock=A,odd_steps=d,words=len(words)))
        P,Q,B=clock(w);modulus=2*Q
        n=((Q-B)*pow(P,-1,modulus))%modulus
        realized,y=native(n,len(w))
        need(n>0 and realized==w,'positive source realizes shortest return witness')
        need((y-n)%q==0,'marked endpoint returns modulo q')
        weight=Fraction(4,(n+1)*(n+3))
        out.append(dict(q=q,minimum_binary_clock=C,word=w,source=n,endpoint=y,
                        source_weight=[weight.numerator,weight.denominator],
                        lower_clock_candidates=candidates,literal_lower_words=tested))
    return out

def debt_control():
    # A deterministic hostile corpus: every word through total binary clock 12,
    # and actual positive-source prefixes through 48 steps.
    cases=[w for A in range(1,13) for d in range(1,A+1) for w in compositions(A,d)]
    cases += [native(n,48)[0] for n in range(1,2002,2)]
    for w in cases:
        A=sum(w);gs=[content(w[:j]) for j in range(1,len(w)+1)]
        for E in range(1,5):
            score=sum(sum(g%5**e==0 for e in range(1,E+1)) for g in gs)
            need(8*5**(E-1)*score <= (5**E-1)*A,'four-layer sharp linear debt')
        for p in (5,7,11):
            q=p;valuation_sum=0;packing=0;e=1
            while q<=3**A-2**A:
                c=1
                while 3**c<=q:c+=1
                if p==5 and e<=4:c=max(c,(2,8,26,33)[e-1])
                visits=sum(g%q==0 for g in gs)
                need(visits<=A//c,'prime-power first-return packing')
                valuation_sum+=visits;packing+=A//c
                q*=p;e+=1
            need(valuation_sum<=packing,'uniform fixed-prime debt bound')
    return len(cases)

def run():
    for d in range(1,8):
        for w in product(range(1,5),repeat=d):
            P,Q,B=clock(w);need((P,Q,B)==closed(w),'two carry constructions')
            g=content(w);need(gcd(g,6)==1,'content prime support')
            need(B<=3**(d-1)+2**(sum(w)-d+1)*(3**(d-1)-2**(d-1)),'carry budget')
            for j in range(1,d):
                u,v=w[:j],w[j:]
                need(gcd(content(u),g)==gcd(content(u),content(v)),'concatenation renewal')
            if d>1:need(gcd(content(w[:-1]),g)==1,'adjacent prefix contents coprime')
    sources=0
    for n in range(1,10002,2):
        w,y=native(n,24);P,Q,B=clock(w)
        need(P*n+B==Q*y,'native source endpoint')
        need(sum(w)<=48+n.bit_length()-1,'positive-source native clock bound')
        need(content(w)==gcd(abs(P-Q),abs(y-n)),'source return content')
        sources+=1
    phases=[];orders=[];returns=[]
    for e in range(1,5):
        q=5**e;period=4*5**(e-1)
        powers=[pow(2,a,q) for a in range(period)]
        need(len(set(powers))==period,'2 primitive for tested five powers')
        phase=powers.index((-3)%q)
        phases.append(dict(q=q,phase=phase,modulus=period))
        for a in range(1,period+1):
            # Test all b modulo the period through e=3; boundary sweep at e=4.
            bs=range(1,period+1) if e<=3 else (1,phase,phase+period,period)
            for b in bs:
                need((content((a,b))%q==0)==(a%period==phase and b%period==phase),'two-letter five-adic gate')
        actual=[order_affine(q,a) for a in (1,2,3)]
        need(actual==[2*5**(e-1),4*5**max(e-2,0),5**e],'constant letter full clocks')
        orders.append(dict(q=q,letters=[1,2,3],orders=actual))
        if e<=3:
            for alphabet in ((1,2),(1,2,3),(1,2,3,4)):
                returns.append(dict(q=q,alphabet=alphabet,odd_steps=shortest(q,alphabet),
                                    binary_clock=shortest(q,alphabet,True)))
    # Positive actual cycle: periodic content is not a positive descent certificate.
    for t in range(1,21):
        w=(2,)*t;P,Q,B=clock(w)
        need(B==Q-P and content(w)==Q-P,'content inflates on positive 1-cycle')
    # Reduced expanding anchor can have arbitrarily many return contents too.
    for t in range(1,21):
        w=(1,)*t;P,Q,B=clock(w)
        need(B==P-Q and content(w)==P-Q,'content inflates at negative fixed anchor')
    costs=minimal_cost_certificates();debt_cases=debt_control()
    return dict(status='FINITE-EXACT controls for the proved algebra; no universal Collatz bound',
                checks=checks,word_universe='all positive valuation words in {1,2,3,4} of lengths 1..7',
                actual_sources=sources,native_steps_per_source=24,phases=phases,orders=orders,returns=returns,
                unrestricted_cost_certificates=costs,debt_control_words=debt_cases,
                hostile='constant words (2)^t and (1)^t have unbounded raw content at fixed reduced anchor')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=Path(__file__).with_name('content_receipt.json'))
    args=ap.parse_args();answer=run();args.output.write_text(json.dumps(answer,indent=2,sort_keys=True)+'\n');print(json.dumps(answer,indent=2))
