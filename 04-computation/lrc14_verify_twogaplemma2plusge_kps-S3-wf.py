"""
Adversarial verification of the "two-gap-lemma2plus" S3 advance for LRC(14).

EXACT Fractions throughout. We verify each load-bearing LEMMA from scratch and
hunt for counterexamples to (a) the lemmas, (b) the claim "Lemma A fires => M>=1/14",
and (c) LRC(14) itself on S3 (any covering S3 13-set with exact M<1/14 would REFUTE).
"""
import sys
from fractions import Fraction as F
from math import gcd, ceil, floor
import random

def flush(*a):
    print(*a); sys.stdout.flush()

# ---------------- EXACT M / W TOOLS ----------------
def nrm(x):
    r = x - int(x); r = r + 1 if r < 0 else r
    return r if r <= F(1,2) else 1 - r

def safe_components(A, h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

def Wwidth(A):
    sc=safe_components(A)
    if not sc: return F(0)
    ws=[b-a for a,b in sc]
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1: ws.append((sc[0][1])+(1-sc[-1][0]))
    return max(ws)

def cand(S):
    S=sorted(set(S)); C=set()
    for v in S:
        k=0
        while F(2*k+1,2*v)<=F(1,2): C.add(F(2*k+1,2*v)); k+=1
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            for d in (S[i]+S[j],S[j]-S[i]):
                if d>0:
                    k=1
                    while F(k,d)<=F(1,2): C.add(F(k,d)); k+=1
    C.add(F(1,2)); return C

def Mval(S):
    b=F(0)
    for t in cand(S):
        v=min(nrm(x*t) for x in S)
        if v>b: b=v
    return b

def is_covering(S):
    return all(any(v%q==0 for v in S) for q in range(2,15))

H = F(1,14)

def gcd_list(xs):
    g=0
    for x in xs: g=gcd(g,x)
    return g

# ================= T1: structural monotonicity =================
def test_structural(trials=4000, seed=1):
    rnd=random.Random(seed); bad=0
    for _ in range(trials):
        tau=F(rnd.randint(1,10000), rnd.randint(1,10000))
        us=sorted(rnd.sample(range(1,2000), rnd.randint(2,12)))
        prev=None
        for u in us:
            ut=u*tau
            fl=ut.numerator//ut.denominator
            if prev is not None and fl < prev:
                bad+=1; break
            prev=fl
    return bad

# ================= LA: Lemma A =================
def lemmaA_window(K, Vmin, Vmax):
    return F(14*K+1, 14*Vmin), F(14*K+13, 14*Vmax)

def test_lemmaA_safety(trials=4000, seed=2):
    rnd=random.Random(seed); bad=[]
    for _ in range(trials):
        Vmin=rnd.randint(14,3000); Vmax=Vmin+rnd.randint(0,200); K=rnd.randint(0,30)
        lo,hi=lemmaA_window(K,Vmin,Vmax)
        if not (lo<hi): continue
        for _ in range(3):
            fr=F(rnd.randint(1,999),1000); tau=lo+fr*(hi-lo)
            for u in (Vmin, Vmax, rnd.randint(Vmin,Vmax)):
                if nrm(u*tau) < H:
                    bad.append((K,Vmin,Vmax,u,float(tau),float(nrm(u*tau))))
    return bad

def test_lemmaA_nonempty_identity(trials=20000, seed=3):
    rnd=random.Random(seed); bad=0
    for _ in range(trials):
        Vmin=rnd.randint(1,5000); Vmax=rnd.randint(1,5000); K=rnd.randint(0,50)
        lo,hi=lemmaA_window(K,Vmin,Vmax)
        nonempty = (lo<hi)
        ineq = (13*Vmin - Vmax) > 14*K*(Vmax-Vmin)
        if nonempty != ineq: bad+=1
    return bad

# ================= LB: Lemma B tooth count =================
def test_lemmaB_toothcount(trials=4000, seed=4):
    rnd=random.Random(seed); bad=0
    for _ in range(trials):
        u=rnd.randint(14,5000)
        a=F(rnd.randint(0,1000),1001); Wp=F(rnd.randint(1,500),5000); b=a+Wp
        lo=ceil(a*u); hi=floor(b*u)
        cnt=hi-lo+1 if hi>=lo else 0
        bound=(b-a)*u
        bf=bound.numerator//bound.denominator
        if cnt> bf+1: bad+=1
    return bad

# ================= Lemma A fires => M>=1/14 =================
def lemmaA_fires(P, L):
    Vmin=min(L); Vmax=max(L); Psafe=safe_components(P,H); s=Vmax-Vmin
    if s==0:
        Kmax=2 if (13*Vmin-Vmax)>0 else -1
    else:
        if 13*Vmin-Vmax<=0: return None
        Kmax=(13*Vmin-Vmax-1)//(14*s)
    for K in range(0,Kmax+1):
        lo,hi=lemmaA_window(K,Vmin,Vmax)
        if not (lo<hi): continue
        for (a,b) in Psafe:
            ov_lo=max(lo,a); ov_hi=min(hi,b)
            if ov_lo<ov_hi: return (ov_lo+ov_hi)/2
    return None

# ================= adversarial S3 generator =================
def gen_S3_sets(n=1500, seed=7, maxV=12000):
    rnd=random.Random(seed); out=[]; attempts=0
    while len(out)<n and attempts<n*400:
        attempts+=1
        psz=rnd.randint(2,9)
        P=sorted(rnd.sample(range(1,14), psz))
        csz=rnd.randint(2, max(2,13-len(P)))
        if csz<2: continue
        V0=rnd.randint(20, maxV); spread=rnd.randint(1, 60)
        L=sorted(set(V0+rnd.randint(0,spread) for _ in range(csz*4)))[:csz]
        if len(L)<2: continue
        S=sorted(set(P)|set(L))
        if len(S)<13: continue
        S=S[:13]
        if len(set(S))!=13: continue
        if gcd_list(S)!=1: continue
        if not is_covering(S): continue
        Vmin=min(S); Vmax=max(S)
        k=sum(1 for v in S if v>13)
        if k<2: continue
        if Vmax < 13*Vmin: continue
        Pp=[v for v in S if v<=13]; Ll=[v for v in S if v>13]
        if len(Ll)<2: continue
        out.append((Pp,Ll))
    return out

if __name__=="__main__":
    flush("=== T1 structural monotonicity (4000 trials) ===")
    flush("  violations:", test_structural())

    flush("=== LA safety (4000 trials) ===")
    bs=test_lemmaA_safety(); flush("  safety violations:", len(bs), bs[:3])

    flush("=== LA nonemptiness identity (20000 trials) ===")
    flush("  identity mismatches:", test_lemmaA_nonempty_identity())

    flush("=== LB tooth count (4000 trials) ===")
    flush("  tooth-count violations:", test_lemmaB_toothcount())

    flush("=== generating S3 sets ===")
    sets=gen_S3_sets(n=1500)
    flush("  S3 sets generated:", len(sets))

    flush("=== LA fire => M>=1/14, and global min M ===")
    fired=0; fire_fail=[]; below=[]; mn=None
    for (P,L) in sets:
        S=sorted(set(P)|set(L))
        w=lemmaA_fires(P,L)
        m=Mval(S)
        if mn is None or m<mn: mn=m
        if m<H: below.append((P,L,m))
        if w is not None:
            fired+=1
            wsafe=all(nrm(v*w)>=H for v in S)
            if (not wsafe) or m<H:
                fire_fail.append((P,L,float(w),wsafe,m))
    flush("  Lemma A fired on:", fired, "of", len(sets))
    flush("  fire-implies-M FAILURES:", len(fire_fail), fire_fail[:3])
    flush("  min exact M over all S3:", mn, float(mn) if mn else None, " (1/14=",float(H),")")
    flush("  S3 sets with M<1/14 (REFUTE LRC14):", len(below), below[:5])
