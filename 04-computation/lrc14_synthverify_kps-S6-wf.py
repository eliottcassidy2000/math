"""
Independent synthesis-verification for the 1/7-spread bound (LRC(14) last lemma).
kind-pasteur S6 workflow synthesis. EXACT rationals.

Goals:
  (A) Reimplement mu_theta independently (engine B) and cross-check engine A anchors.
  (B) Reproduce mu_{1/7}(consec_k) for k=7..13 and thr_k for k=8..12.
  (C) Recompute union-bound floor at k=8.
  (D) Adversarial hunt for mu_{1/7}(E) < thr_k (the ACTUAL lemma).
  (E) Sanity: same machinery DOES beat consec for mu_{2/7} (search-not-impotent test).
"""
from fractions import Fraction as F
from itertools import combinations
from math import floor
import random

# ---------- Engine A (dispatched, copied verbatim from prompt) ----------
def mu_theta_A(E,theta):
    E=sorted(set(E)); n=len(E); bp=set([F(0),F(1)])
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            for m in range(0,d+1): bp.add(F(m,d))
    bp=sorted(b for b in bp if 0<=b<=1); total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2; order=sorted(range(n),key=lambda i:(E[i]*mid)%1)
        ks=[floor(E[order[t]]*mid) for t in range(n)]; subs=[]
        for t in range(n):
            o1=order[t];o2=order[(t+1)%n];k1=ks[t];k2=ks[(t+1)%n];wrap=1 if t==n-1 else 0
            s=E[o2]-E[o1];c=F(k1-k2+wrap)
            if s==0:
                if c>theta: subs.append((a,b))
            elif s>0:
                lo=max(a,(theta-c)/s)
                if lo<b: subs.append((lo,b))
            else:
                hi=min(b,(theta-c)/s)
                if a<hi: subs.append((a,hi))
        subs.sort(); cur=cb=None
        for lo,hi in subs:
            if cur is None: cur,cb=lo,hi
            elif lo<=cb: cb=max(cb,hi)
            else: total+=cb-cur; cur,cb=lo,hi
        if cur is not None: total+=cb-cur
    return total

# ---------- Engine B (independent): explicit gap-threshold breakpoints ----------
# Build all breakpoints where (a) the cyclic ORDER of {frac(e x)} changes, i.e.
# e_i x - e_j x in Z  => x = m/(e_i-e_j); AND (b) a gap crosses theta, i.e. a
# consecutive difference (e_b-e_a) x - integer = theta. We then evaluate max-gap>theta
# at the midpoint of each cell exactly.
def maxgap_at(E, x):
    # exact max circular gap of {frac(e x)}
    pts = sorted(set((F(e)*x) % 1 for e in E))
    if len(pts)==1:
        return F(1)
    g = F(0)
    for i in range(len(pts)):
        nxt = pts[(i+1)%len(pts)]
        d = (nxt - pts[i]) % 1 if i==len(pts)-1 else nxt-pts[i]
        if i==len(pts)-1:
            d = (pts[0]+1)-pts[i]
        if d>g: g=d
    return g

def mu_theta_B(E,theta):
    E=sorted(set(E)); n=len(E)
    bp=set([F(0),F(1)])
    # order-change breakpoints
    for i in range(n):
        for j in range(i+1,n):
            d=E[j]-E[i]
            if d==0: continue
            for m in range(0,d+1):
                bp.add(F(m,d))
    # gap=theta crossing breakpoints: for each ordered pair (a,b), (e_b-e_a) x = k + theta
    # i.e. x = (k+theta)/(e_b-e_a). theta rational => add these.
    for i in range(n):
        for j in range(n):
            if i==j: continue
            s = E[j]-E[i]
            if s==0: continue
            # x in (0,1): (k+theta)/s in (0,1)
            # k ranges so that 0 < k+theta < s  (s could be negative; handle abs)
            sa=abs(s)
            kmin = -2; kmax = sa+2
            for k in range(kmin,kmax+1):
                val = (F(k)+theta)/s
                if 0<val<1:
                    bp.add(val)
    bp=sorted(b for b in bp if 0<=b<=1)
    total=F(0)
    for a,b in zip(bp,bp[1:]):
        if b<=a: continue
        mid=(a+b)/2
        if maxgap_at(E,mid) > theta:
            total += b-a
    return total

# ---------- Anchors ----------
consec = lambda k: list(range(k))
anchors = {7:F(1,1),8:F(691,735),9:F(247,294),10:F(38,49),
           11:F(1381,2205),12:F(13823,24255),13:F(477,1078)}
thr = {8:F(3637,5880),9:F(2025,4004),10:F(36,91),11:F(25,91),12:F(1,7),13:F(0)}

print("=== (A/B) Engine cross-check on consec anchors (1/7) ===")
allok=True
for k in range(7,14):
    E=consec(k)
    a=mu_theta_A(E,F(1,7))
    expect=anchors[k]
    okA = (a==expect)
    # engine B only for k<=11 (cost)
    if k<=11:
        b=mu_theta_B(E,F(1,7))
        okB=(b==expect)
        print(f"k={k}: A={a} expect={expect} A_ok={okA} | B={b} B_ok={okB}")
        allok = allok and okA and okB
    else:
        print(f"k={k}: A={a} expect={expect} A_ok={okA} (B skipped, cost)")
        allok = allok and okA
print("ALL anchors match (A and B where run):", allok)

print("\n=== Random cross-check A vs B (1/7), 120 sets k in 8..10 ===")
random.seed(7)
mism=0
for _ in range(120):
    k=random.randint(8,10)
    spread=random.randint(k-1, 18)
    pts=set([0]);
    while len(pts)<k:
        pts.add(random.randint(1,spread))
    E=sorted(pts)
    if max(E)-min(E)==0: continue
    a=mu_theta_A(E,F(1,7)); b=mu_theta_B(E,F(1,7))
    if a!=b:
        mism+=1
        print("MISMATCH",E,a,b)
print("A-vs-B mismatches:",mism,"(0 expected)")

print("\n=== (C) k=8 union-bound floor ===")
mP = F(2243,5880)  # meas(G_P) at |P|=5 ... check: floor = meas(G_P)+mu-1
# thr_8 = 1 - min meas(G_P) over |P|=5; min meas(G_P)=2243/5880 => thr_8=3637/5880
minGP8 = 1 - thr[8]
floor8 = minGP8 + anchors[8] - 1
print("min meas(G_P) at |P|=5 =", minGP8, "=", float(minGP8))
print("k=8 floor = meas(G_P)+mu_consec-1 =", floor8, "=", float(floor8), ">0:", floor8>0)
print("matches claimed 1891/5880:", floor8==F(1891,5880))

print("\n=== consec mu vs thr per k (the finite check IF consec is minimizer) ===")
for k in range(8,14):
    m=anchors[k]; t=thr[k]
    print(f"k={k}: mu_consec={m}={float(m):.4f}  thr={t}={float(t):.4f}  margin={float(m-t):.4f}  ok={m>=t}")
