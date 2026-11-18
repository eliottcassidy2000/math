"""Replay a switch balanced at two slots but not at one slot."""
from functools import cache
from fractions import Fraction as F
from math import comb
from padding_replay import slots, rooted
from frame_replay import join, add, cut
from replay import pad

U=(24,109,121,80,2,25,49)
L=pad((6,292,88,19,163,333,155,47,221),1)
R=(0,161,35,7,15,217,31,349,127)
I=3

def templates(U,I,q,flip):
    m=len(U)
    J=join((0,),join(U,(0,)))
    Q=add(J,I<<1)
    if flip:Q=cut(Q,1<<(m+2))
    j=m+1
    @cache
    def count(a,b,c,t,p):
        if not (a or b or c):return 1
        w=count(a,b,c,t-1,-1) if t>1 else 0
        T=(1<<(m+3))-1 if p<0 else Q[p]
        for i in range(m+3):
            if not T>>i&1:continue
            if i==0:
                if a:w+=count(a-1,b,c,t,i)
            elif i==j:
                if b:w+=count(a,b-1,c,t,i)
            else:
                if c>>i&1:w+=count(a,b,c^(1<<i),t,i)
        return w
    mask=((1<<(m+3))-1)^1^(1<<j)
    return {(a,b):k for a in range(1,q+2) for b in range(1,q+2)
            if (k:=count(a,b,mask,q,-1))}

def eval_slots(weights,x,y):
    C=lambda v,k:sum((-1)**(k-h)*comb(k,h)*v[h] for h in range(k+1))
    return sum(w*C(x,a)*C(y,b) for (a,b),w in weights.items())

def shell(q,flip):
    x=[z[-1] for z in slots(L,3)]
    y=[z[-1] for z in slots(R,3)]
    return eval_slots(templates(U,I,q,flip),x,y)

def profiles():
    a,b=([z[-1] for z in slots(P,3)] for P in (L,R))
    assert a==[0,1659,58560,777555] and b==[0,207,7920,105165]
    E,S=rooted(U,I)[1:3]
    assert E==(259,211) and S==(1619,1619)
    h,H=shell(1,False),shell(1,True)
    B,BB=shell(2,False),shell(2,True)
    assert H-h==367189200 and B==BB
    assert F(S[0]*(a[3]*b[2]-b[3]*a[2]),a[2]*b[2])==F(h-H,a[2]*b[2])
    print((h,H,B))
    return h,H,B

if __name__=="__main__":profiles()
