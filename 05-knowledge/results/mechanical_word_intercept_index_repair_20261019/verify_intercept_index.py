"""Focused exact regression for THM-4716 proof 1; no whole-theorem audit."""
from math import gcd

def need(claim, message):
    if not claim:
        raise RuntimeError(message)

def word(q,p,t):
    return tuple(((k+1)*p+t)//q-(k*p+t)//q for k in range(q))

checks=0
wrong_original=0
first=None
for q in range(2,65):
    for p in range(1,q):
        if gcd(p,q)!=1:
            continue
        inv=pow(p,-1,q)
        for t in range(q-1):
            # Direct floor increments are independent of the modular-inverse formula.
            changed=[k for k in range(q+1) if (k*p+t+1)//q != (k*p+t)//q]
            k1=(inv*(q-1-t))%q
            need(changed==[k1],('floor-change index',q,p,t,changed,k1))
            need(1<=k1<=q-1,('nonwrap',q,p,t,k1))
            before=word(q,p,t);after=word(q,p,t+1)
            changed_letters=[k for k in range(q) if before[k]!=after[k]]
            need(changed_letters==[k1-1,k1],('adjacent support',q,p,t,changed_letters,k1))
            need(before[k1-1:k1+1]==(0,1) and after[k1-1:k1+1]==(1,0),('swap direction',q,p,t))
            if k1!=(q-1-t)%q:
                wrong_original+=1
                if first is None:first=(q,p,t,k1,(q-1-t)%q)
            checks+=1
# Mandatory hostile, plus p=1 and both intercept ends as controls.
need((pow(2,-1,5)*4)%5==2,'(5,2,0) hostile')
need(word(5,2,0)==(0,0,1,0,1),'hostile old word')
need(word(5,2,1)==(0,1,0,0,1),'hostile new word')
for q in range(2,65):
    for t in range(q-1):need((pow(1,-1,q)*(q-1-t))%q==q-1-t,'p=1 control')
print('Universe: 2 <= q <= 64, 1 <= p < q, gcd(p,q)=1, 0 <= t <= q-2.')
print(f'Exact intercept transitions checked: {checks}.')
print(f'Incorrect original indices: {wrong_original}. First witness (q,p,t,actual,original): {first}.')
print('Specified hostile (q,p,t)=(5,2,0): actual k1=2; original formula gives 4.')
print('PASS: modular-inverse index, direct floor increments, nonwrapping support, and 01 -> 10 swap agree.')
print('Scope: proof 1 index and adjacent-swap lemma only; no whole-theorem audit.')
