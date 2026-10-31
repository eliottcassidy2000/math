#!/usr/bin/env python3
"""Exact E13/E14 arithmetic and literal native-cut semantics; bounded controls."""
from fractions import Fraction as F
from pathlib import Path
from math import gcd
import json


def need(x, msg):
    if not x: raise ValueError(msg)


def height(x, cap=4):
    if x==0:return cap
    j=0
    while j<cap and x%5==0:j+=1;x//=5
    return j


def lam(e):
    return F(5,8)*(1-F(1,5**e))


audit=json.loads(Path('adaptive10_audit.json').read_text())
need(audit['pass'] and audit['all_positive_letters'], 'missing all-letter audit')
eta=F(-audit['terminal'],audit['D'])
beta=F(-audit['max'],audit['D'])
gamma=min(eta-2,beta)
need(gamma>1 and gamma==F(61035156,48828125),'E12 non-one margin')
margin13=gamma-1
native13=eta-3
need(native13>F(7,10),'native13 lower floor')
remaining=[(a,b) for b in range(1,5) for a in range(5) if b*native13+F(a,5)<1]
need(remaining==[(0,1),(1,1)],'exact residual families')
T13=2*5**12
ones_charge=sum(T13//(2*5**(e-1)) for e in range(1,14))
need(ones_charge==lam(13)*T13,'sharpness checksum')

words=0
cuts={e:0 for e in range(2,5)}
transitions=0
def dfs(w, P, B, Q, A, prefix_heights):
    global words, transitions
    for a in range(1,17-A):
        pp,bb,qq=3*P,3*B+Q,Q*2**a
        c=height(gcd(pp-qq,bb))
        row=w+(a,);cs=prefix_heights+[c];clock=A+a
        words+=1
        for E in range(1,5):
            q=5**E
            z=(B-P+Q)*pow(Q,-1,q)%q
            z2=(bb-pp+qq)*pow(qq,-1,q)%q
            m=P*pow(Q,-1,q)%q;m2=pp*pow(qq,-1,q)%q
            for R in range(1,E+1):
                oldmod=5**max(R,height(z,E));newmod=5**max(R,height(z2,E))
                common=min(oldmod,newmod)
                need(m2%common==3*(m%oldmod)*pow(2,-a,common)%common,'actual multiplier retained')
                need(min(height(z2,E),height(m2%newmod-1,E))==min(c,E),'actual content retained')
                transitions+=1
        for e in range(2,5):
            if c<e or any(x>=e for x in prefix_heights):continue
            # Rebuild each lower-return block literally, starting from identity.
            p,b,q,lower,block_clock=1,0,1,0,0
            for letter,global_c in zip(row,cs):
                p,b,q=3*p,3*b+q,q*2**letter;block_clock+=letter
                local=height(gcd(p-q,b))
                need(min(local,e-1)==min(global_c,e-1),'entered source lost')
                lower+=min(local,e-1)
                if global_c>=e-1:p,b,q=1,0,1
            need((p,b,q)==(1,0,1),'native block did not finish')
            need(sum(min(x,e) for x in cs)==lower+1,'native extra charge not one')
            need(block_clock==clock,'integer clock lost')
            cuts[e]+=1
        dfs(row,pp,bb,qq,clock,cs)
dfs((),1,0,1,0,[])
print(json.dumps(dict(status='FINITE-EXACT proof arithmetic and semantic controls; not enumeration of all native13 words',
    eta10=str(eta),gamma12=str(gamma),E13_coefficient=str(lam(13)),E13_nonone_margin=str(margin13),
    native13_floor=str(native13),remaining_E14=[{'ones':a,'nonone':b} for a,b in remaining],
    E13_ones_clock=T13,E13_ones_charge=ones_charge,words_clock_at_most_16=words,
    all_minimum_precision_transition_checks=transitions,literal_native_cut_checks=cuts,pass_checks=True),indent=2))
