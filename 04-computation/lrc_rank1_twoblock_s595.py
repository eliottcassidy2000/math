#!/usr/bin/env python3
"""The large-owner cover as a 2x2 determinant (TWO-BLOCK) + bounded CRT automaton.
Component C_i=(a,b) of G(S'), owners (u_a,k_a,+1),(u_b,k_b,-1); v=nw arc index j.
Endpoint congruences r_a=w(k_a n+1)-j u_a (|r_a|<u_a/n), r_b=w(k_b n-1)-j u_b.
ELIMINATE j: det[[u_a, r_a],[u_b, r_b]] = u_a r_b - u_b r_a = w * n * u_a u_b * (b-a).
So cover <=> the slack 2-vector (r_a,r_b) makes the OWNER|SLACK 2x2 have det = w n u_a u_b l_i
(the obstruction: rank-1 [det=0, slacks parallel to owners] can't reach the nonzero target).
opus-2026-06-03-S595. Verify the determinant identity + the bounded automaton emptiness."""
from fractions import Fraction as F
from math import gcd
def dist(x): x%=1; return min(x,1-x)
def G_components(Sp,n):
    THR=F(1,n); pts={}
    for u in Sp:
        for k in range(u+1):
            for eps in (1,-1): pts.setdefault(F(k*n+eps,n*u)%1,[]).append((u,k,eps))
    order=sorted(pts); comps=[]; L=len(order)
    for i in range(L):
        a=order[i]; b=order[(i+1)%L]; ln=(b-a) if b>a else (b-a+1); mid=(a+ln/2)%1
        if all(dist(u*mid)>THR for u in Sp): comps.append((a,b,ln,pts[a],pts[b]))
    return comps
def main():
    print("(1) Verify the 2x2-determinant (two-block) identity for a COVERED component:")
    # use a multiple-of-n config; find a component covered by some v-arc and check the identity
    n=8; w=3; v=n*w   # v=24
    Sp=(1,2,3,5,11,13,15)  # n-1=7 other runners, no multiple of 8
    for (a,b,ln,oa,ob) in G_components(Sp,n)[:6]:
        la=[o for o in oa if o[2]==1]; rb=[o for o in ob if o[2]==-1]
        if not la or not rb: continue
        ua,ka,_=la[0]; ub,kb,_=rb[0]
        # is this component covered by v? find j = round(v*a)
        import math
        j=round(float(v*a))
        ra=w*(ka*n+1)-j*ua; rbv=w*(kb*n-1)-j*ub
        lhs=ua*rbv-ub*ra; rhs=w*n*ua*ub*(b-a)
        print(f"   comp len={float(ln):.4f} owners u_a={ua} u_b={ub}: det[[u_a,r_a],[u_b,r_b]]={lhs}; w*n*u_a*u_b*(b-a)={rhs}; identity={lhs==rhs}")
    print()
    print("(2) Bounded CRT automaton (per-component allowed w) for mult-of-n residual rows; intersection empty => loose")
    import random; rng=random.Random(1)
    for n in [10,12,14]:
        m=n-1; B=2*n+8; emptyc=0; tot=0
        for _ in range(400):
            others=set(rng.sample([x for x in range(1,B+1) if x%n!=0],m-1))
            ww=rng.randint(1,3); v=n*ww
            if v in others: continue
            V=tuple(sorted(others|{v}))
            g=0
            for x in V: g=gcd(g,x)
            V=tuple(sorted(x//g for x in V))
            if len(V)!=m or not any(x%n==0 for x in V): continue
            mults=[x for x in V if x%n==0]; vv=mults[0]; Sp=tuple(x for x in V if x!=vv)
            comps=G_components(Sp,n)
            if not comps: continue
            tot+=1
            # bounded automaton: allowed w in [1, wmax] s.t. EVERY component coverable by v=nw
            wmax=( (n-1)*max(Sp) )//n + 1
            allowed=set(range(1,min(wmax,200)+1))
            for (a,b,ln,oa,ob) in comps:
                mid=(a+ln/2)%1
                allowed={w for w in allowed if dist(n*w*mid)<=F(1,n)-F(n*w,2)*ln}
                if not allowed: break
            if not allowed: emptyc+=1
        print(f"   n={n}: residual rows={tot}; bounded-CRT-automaton EMPTY (=> loose)={emptyc}/{tot}")
if __name__=='__main__': main()
