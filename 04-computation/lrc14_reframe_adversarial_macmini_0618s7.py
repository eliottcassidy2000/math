#!/usr/bin/env python3
"""
lrc14_reframe_adversarial_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

ADVERSARIAL verification of the tournament-reframe dictionary claims before canonizing.

CLAIMS TO STRESS:
 (V1) R3 universal: maxgap(x)>1/7 <=> exists empty 1/7-arc just-after some phase.
      Test on HARD shapes: large dissociated, near-AP, random, including k up to 10.
 (V2) AP-orbit invariance is EXACT: meas(S7), measN, E[c3], E[H] are identical for E and
      dilated dE (gcd-coprime). Test E vs 3E vs 5E with EXACT Fraction equality.
 (V3) R2 theorem: maxgap>1/2 <=> Condorcet winner. Stress with k up to 10.
 (V4) The 'mu_{1/7} = P[scale-1/7 local sink]' reformulation matches the existing
      good_set engine's measN EXACTLY (cross-check two independent implementations).
 (V5) HONEST: is the co-monotonicity meas(S7) ~ E[H] a THEOREM or just a trend?
      Count rank-inversions on a big exhaustive batch (k=7). Report exactly.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def round_tournament(E,x):
    n=len(E); A=[[0]*n for _ in range(n)]; tf=True
    for i in range(n):
        for j in range(n):
            if i==j: continue
            rel=((E[i]-E[j])*x)%1
            if 0<rel<F(1,2): A[i][j]=1
            elif rel>F(1,2): A[i][j]=0
            else: A[i][j]=1 if E[i]<E[j] else 0; tf=False
    return A,tf
def scores(A): return [sum(r) for r in A]
def c3f(A):
    n=len(A); return comb(n,3)-sum(comb(s,2) for s in scores(A))
def Hpaths(A):
    n=len(A); h=0
    for p in itertools.permutations(range(n)):
        if all(A[p[t]][p[t+1]] for t in range(n-1)): h+=1
    return h
def maxgap(E,x):
    ps=sorted(set((e*x)%1 for e in E))
    if len(ps)==1: return F(1)
    g=F(0)
    for i in range(len(ps)):
        nxt=ps[(i+1)%len(ps)]+(F(1) if i+1==len(ps) else F(0)); g=max(g,nxt-ps[i])
    return g
def empty17(E,x):
    P=sorted(set((e*x)%1 for e in E))
    for a in range(len(P)):
        lo=P[a]; cnt=0
        for b in range(len(P)):
            if b==a: continue
            d=(P[b]-lo)%1
            if F(0)<d<=F(1,7): cnt+=1
        if cnt==0: return True
    return False
def condorcet(E,x):
    A,_=round_tournament(E,x); s=scores(A); return (len(E)-1) in s
def sectors_hit(E,x): return set(int(((e*x)%1)*7) for e in E)
def breakpoints(E):
    bps=set([F(0),F(1)]); Es=sorted(set(E))
    for e in Es:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
        for m in range(0,2*e+1): bps.add(F(m,2*e))
    diffs=set()
    for a in range(len(Es)):
        for b in range(a+1,len(Es)):
            diffs.add(Es[b]-Es[a]); diffs.add(Es[a]+Es[b])
    for d in diffs:
        if d==0: continue
        for m in range(0,2*d+1): bps.add(F(m,2*d))
    return sorted(x for x in bps if 0<=x<=1)

def profile(E, do_H=True):
    k=len(E); bps=breakpoints(E)
    S7=F(0); N=F(0); Ec3=F(0); EH=F(0)
    bad_R3=0; bad_R2=0
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        w=x1-x0; xm=(x0+x1)/2
        A,tf=round_tournament(E,xm)
        mg=maxgap(E,xm)
        if (mg>F(1,7))!=empty17(E,xm): bad_R3+=1
        if tf and (mg>F(1,2))!=condorcet(E,xm): bad_R2+=1
        if len(sectors_hit(E,xm))==7: S7+=w
        if mg<=F(1,7): N+=w
        Ec3+=w*c3f(A)
        if do_H: EH+=w*Hpaths(A)
    return S7,N,Ec3,EH,bad_R3,bad_R2

print("="*92)
print("ADVERSARIAL verification of the tournament-reframe dictionary")
print("="*92)

# V1 + V3: R3 and R2 on hard shapes
print("\n[V1/V3] R3 (1/7-gap=local sink) and R2 (1/2-gap=Condorcet) on hard shapes:")
hard=[("k7 dissoc",[0,1,3,7,15,31,63]),
      ("k8 nearAP",[0,2,3,4,5,6,8,10]),
      ("k9 mixed",[0,1,2,5,11,17,23,30,41]),
      ("k10 dissoc",[0,1,3,7,15,31,63,127,255,511]),
      ("k8 random",[0,5,13,27,41,58,79,97])]
tot_bad3=tot_bad2=0
for name,E in hard:
    doH = len(E)<=8
    S7,N,Ec3,EH,b3,b2=profile(E,do_H=doH)
    tot_bad3+=b3; tot_bad2+=b2
    print(f"  {name:<14} R3_bad={b3} R2_bad={b2}  measN={float(N):.5f} meas(S7)={float(S7):.5f}")
print(f"  TOTAL R3 failures={tot_bad3}  R2 failures={tot_bad2}  (want 0,0)")

# V2: AP-orbit / scale invariance EXACT
print("\n[V2] EXACT scale invariance: profile(E) == profile(dE) for d coprime?")
for base in [[0,1,2,3,4,5,6,7],[0,1,3,7,15,31,63],[0,2,3,4,5,6,8]]:
    p0=profile(base, do_H=(len(base)<=8))[:4]
    ok=True
    for d in (3,5,11):
        pd=profile([d*e for e in base], do_H=(len(base)<=8))[:4]
        if pd!=p0: ok=False; break
    print(f"  base {base}:  scale-invariant EXACT? {ok}  (meas(S7)={p0[0]}, measN={p0[1]})")

# V4: cross-check measN against the canonical good_set engine
print("\n[V4] cross-check measN (our breakpoint engine) vs canonical good_set engine:")
def good_set_canon(E, thr=F(1,7)):
    E=sorted(set(E)); k=len(E); diffs=set()
    for a in range(k):
        for b in range(a+1,k): diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1); good=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        pts=sorted(((e*xm)%1,e) for e in E); order=[e for _,e in pts]
        fl=[int((e*xm)//1) for e in order]
        for idx in range(k):
            ec=order[idx]; fc=fl[idx]
            if idx<k-1: en=order[idx+1]; fn=fl[idx+1]; wrap=F(0)
            else: en=order[0]; fn=fl[0]; wrap=F(1)
            A=F(en-ec); Cc=F(fc-fn)+wrap
            if A==0:
                if Cc>thr: good.append((x0,x1))
                continue
            xb=(thr-Cc)/A
            lo,hi=(max(x0,xb),x1) if A>0 else (x0,min(x1,xb))
            if lo<hi: good.append((lo,hi))
    good=sorted(good); out=[]
    for a,b in good:
        if out and a<=out[-1][1]: out[-1]=(out[-1][0],max(out[-1][1],b))
        else: out.append((a,b))
    return F(1)-sum((b-a for a,b in out),F(0))
for E in [[0,1,2,3,4,5,6,7],[0,1,2,3,4,5,6,8],[0,2,3,4,5,6,8,10]]:
    mine=profile(E,do_H=False)[1]
    canon=good_set_canon(E)
    print(f"  E={E}:  ours measN={mine}  canon measN={canon}  match? {mine==canon}")

# V5: co-monotonicity meas(S7) vs E[H] — count rank inversions, exhaustive k=7
print("\n[V5] HONEST: meas(S7) vs E[H] co-monotone? rank inversions over primitive k=7 box:")
rows=[]
for rest in itertools.combinations(range(1,9),6):
    E=[0]+list(rest)
    if reduce(gcd,E)!=1: continue
    S7,N,Ec3,EH,_,_=profile(E,do_H=True)
    rows.append((E,S7,EH))
inv=0; tot=0
for i in range(len(rows)):
    for j in range(i+1,len(rows)):
        tot+=1
        a=rows[i]; b=rows[j]
        # inversion: S7 strictly orders one way but H the other
        if (a[1]-b[1])*(a[2]-b[2])<0: inv+=1
print(f"  {len(rows)} shapes, {tot} pairs, {inv} S7/H rank inversions "
      f"({100*inv/tot:.2f}% discordant) => co-monotone is a STRONG TREND, not exact.")
print("="*92)
print("DONE.")
