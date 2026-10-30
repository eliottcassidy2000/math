#!/usr/bin/env python3
"""Bounded controls for the phase budget, refuted universally by the successor packet."""
from collections import deque
from fractions import Fraction
from pathlib import Path
import argparse,json
from defect_controls import compositions,valuation

def budget(word):
    P=Q=1;B=A=d=0;epochs=[]
    for a in word:
        P,Q,B=3*P,2**a*Q,3*B+Q;A+=a;d+=1
        if a>1:epochs.append([a,valuation(B+Q-P),valuation(2**a-2),(A+d)%2==0,0])
        elif epochs:epochs[-1][-1]+=1
    debt=sum(max(h-k,0) for a,h,k,active,L in epochs if active)
    paid=sum((Fraction(5*a,8) for a,h,k,active,L in epochs),Fraction())
    for a,h,k,active,L in epochs:
        paid+=Fraction(5*L,8*5**h)-k if active else Fraction(5*L,8)
    return debt,paid

def finite_search(E,cap=12):
    q=5**E;scale=8*q;coef=5*q
    def h(m,b):
        z=(b+1-m)%q
        return next((e for e in range(E) if z%5**(e+1)),E)
    def step(s,a):
        m,b=s;inv=pow(2**a,-1,q);return 3*m*inv%q,(3*b+1)*inv%q
    states={};queue=deque();queued=set();m,b=1,0
    for L in range(2*5**(E-1)):
        for a in range(2,cap+1):
            t=step((m,b),a);kappa=min(valuation(2**a-2),E)
            cost=scale*max(h(*t),kappa)*(t[0]%5 in (1,4))-coef*a
            if t not in states or cost>states[t]:
                states[t]=cost
                if t not in queued:queue.append(t);queued.add(t)
        m,b=step((m,b),1)
    popped=0
    while queue:
        s=queue.popleft();queued.remove(s);score=states[s];popped+=1
        if popped>2000000:raise RuntimeError('finite exploration guard reached')
        for a in range(1,cap+1):
            t=step(s,a)
            if a==1:delta=-coef//5**h(*s) if s[0]%5 in (1,4) else -coef
            else:
                kappa=min(valuation(2**a-2),E)
                delta=scale*max(h(*t),kappa)*(t[0]%5 in (1,4))-coef*a
            score2=score+delta
            if score2>0:raise RuntimeError('positive truncated path; extract exact word before claiming a counterexample')
            if t not in states or score2>states[t]:
                states[t]=score2
                if t not in queued:queue.append(t);queued.add(t)
    return dict(E=E,letters_through=cap,states=len(states),popped=popped,
                maximum_scaled_debt=max(states.values()),scale=scale,
                scope='finite truncated-height graph only; no unbounded-letter or layer claim')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path(__file__).with_name('epoch_budget.json'));a=p.parse_args()
    total=0
    for A in range(1,19):
        for w in compositions(A):
            debt,paid=budget(w)
            if debt>paid:raise RuntimeError(f'exact sufficient-budget failure: {w}: {debt}>{paid}')
            total+=1
    dark=[]
    for r in (0,1,21,46,296):
        w=(2,)+(1,)*(2*r)+(3,);debt,paid=budget(w)
        if debt>paid:raise RuntimeError('phase-aware dark family control')
        dark.append(dict(r=r,debt=debt,budget=str(paid)))
    report=dict(exact_words=total,max_clock=18,dark_controls=dark,finite_searches=[finite_search(E) for E in (1,2,3)],
                status='FINITE-EXACT bounded controls; universal budget REFUTED in phase_budget_shared_orbits_20261030; no universal coverage')
    text=json.dumps(report,indent=2,sort_keys=True)+'\n';a.output.write_text(text);print(text,end='')
