#!/usr/bin/env python3
"""
lrc14_gauss_renorm_catmap_mac-mini-2026-06-20-S6.py   (mac-mini-2026-06-20-S6)

THREAD: Arnold cat map / hyperbolic torus dynamics  <->  three-gap / continued-fraction
(Gauss-map) RENORMALIZATION of the LRC(14) seven-sector cover (THM-536 Sturmian reframe).

HYPOTHESIS (falsifiable):
 (H1) RENORM SELF-SIMILARITY.  The cover set C_k(a) = { a in [0,1) : the AP_k partial-sum
      Sturmian walk S_e = floor(e*theta) mod 7 (theta=integer-part + a) surjects Z/7 }
      has a self-similar structure under the continued-fraction induction (Gauss map
      G(a)={1/a}) -- i.e. the cover near a->0 is a renormalized COPY of the global cover.
 (H2) CONSEC = RENORM FIXED/PERIODIC POINT.  consec_k is invariant (or periodic) under the
      slope-renormalization in a way that spread sets E are NOT, which would explain its
      extremality as renorm-invariance.
 (H3) LYAPUNOV / CONTRACTION.  The wide->iid decoupling (5/7 decoupling; measS7 -> the iid
      surjection prob 7! S(k,7)/7^k) is governed by the Gauss-map mixing/Lyapunov rate
      lambda = pi^2/(6 ln 2) (Lyapunov exp of CF) -- spreading a cluster pushes its slope
      orbit toward equidistribution at this rate.

All EXACT (Fraction). stdlib only.  Connects: THM-536 (Sturmian reframe), THM-555 (cut/cycle
wall, c3), HYP-2703 (7-band slope identity), HYP-2704 (irreducibly aggregate), the 5/7 decouple.
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, log, pi, factorial
sys.stdout.reconfigure(line_buffering=True)

# ---------- exact cover engine (matches lrc14_sector_majorize_macmini_0618s7b.measS7) ----------
def measS7(E):
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if len(res)==7: total+=x1-x0
    return total

def cover_intervals_x(E):
    """Cover intervals in x in [0,1) (theta=7x)."""
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps); ivs=[]
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        res=set(int(7*e*xm)%7 for e in E)
        if len(res)==7: ivs.append((x0,x1))
    merged=[]
    for iv in ivs:
        if merged and merged[-1][1]==iv[0]: merged[-1]=(merged[-1][0],iv[1])
        else: merged.append([iv[0],iv[1]])
    return [tuple(m) for m in merged]

def span(E): return max(E)-min(E)
def primitive(E):
    g=0
    for e in E: g=gcd(g,e)
    return g==1

print("="*78)
print("PART 0: integer-part decomposition of the partial-sum walk (the 7-band split)")
print("="*78)
print("S_e = floor(e*theta) mod 7 with theta=j+a  ==>  S_e = (e*j + floor(e*a)) mod 7.")
print("AP part (step j mod 7) twisted by Sturmian defect floor(e*a). Cover is 7-periodic in j.")
def covered(E, num, den):
    return len(set(((e*num)//den)%7 for e in E))==7
E8=tuple(range(8))
ok7=True
for a_num,a_den in [(1,3),(1,5),(2,5),(3,7),(2,9)]:
    row=[]
    for j in range(0,14):
        th_num=j*a_den+a_num
        row.append(1 if covered(E8,th_num,a_den) else 0)
    per7 = (row[:7]==row[7:14])
    ok7 = ok7 and per7
    print(f"  a={a_num}/{a_den}: cover by j=0..6 {row[:7]}  7-periodic in j? {per7}")
print(f"  ==> cover EXACTLY 7-periodic in integer part j? {ok7}  (each band = Z/7 twist, HYP-2703)")

print()
print("="*78)
print("PART 1 (H1): Gauss-map / continued-fraction self-similarity of the cover set")
print("="*78)
# Work in a = frac(theta) on a single band. Take the band j with j mod 7 = 1 (step-1 AP, the
# 'consec-like' band).  Cover-in-a set: C(a) for AP_k.  Then induce by the Gauss map on the
# slope a near a->0 (where the mechanical word is long runs of the minority letter).
#
# Renormalization claim: zoom the cover into a in [1/(n+1), 1/n] and renormalize by a->{1/a}
# (Gauss).  Compare the induced cover to the un-zoomed cover.  Self-similar => the SHAPE repeats.
def cover_set_in_a(E, band_j, AcDen=None):
    """Return cover intervals (in a in [0,1)) on band theta=band_j+a, exact breakpoints."""
    # theta=band_j+a, a in [0,1).  floor(e*theta)=floor(e*band_j + e*a)= e*band_j + floor(e*a).
    # breakpoints of floor(e*a) in a: m/e for m=1..e-1, plus jumps. Use 1/e grid via e*a integer.
    E=sorted(set(E)); Enz=[e for e in E if e!=0]
    bps=set([F(0),F(1)])
    for e in Enz:
        for m in range(1,e): bps.add(F(m,e))
    bps=sorted(bps); ivs=[]
    for i in range(len(bps)-1):
        a0,a1=bps[i],bps[i+1]
        am=(a0+a1)/2
        res=set((e*band_j + int(e*am))%7 for e in E)
        if len(res)==7: ivs.append((a0,a1))
    merged=[]
    for iv in ivs:
        if merged and merged[-1][1]==iv[0]: merged[-1]=[merged[-1][0],iv[1]]
        else: merged.append([iv[0],iv[1]])
    return [tuple(m) for m in merged]

for k in [7,8,9]:
    E=tuple(range(k))
    for bj in [1,2,3]:
        ivs=cover_set_in_a(E,bj)
        tot=sum(b-a for a,b in ivs)
        print(f"  k={k} band j={bj}: cover-meas-in-a={float(tot):.5f}  #ivs={len(ivs)}")
# Self-similarity probe: for AP_k, look at the cover near a->0 vs near a->1 and at Gauss preimage.
print("\n  Gauss-induction probe (k=8, band j=1): cover near a in (0,1) and its G-renormalization")
E=tuple(range(8))
ivs=cover_set_in_a(E,1)
near0=[(float(a),float(b)) for a,b in ivs if a<F(1,4)]
print(f"   cover intervals with a<1/4: {near0}")
# Gauss-renormalize: a in (1/(n+1),1/n) -> 1/a - n in (0,1).  Check if the cover in successive
# Farey/CF windows (1/(n+1),1/n) is a SCALED copy.
print("   Gauss windows (1/(n+1),1/n): cover fraction within each window (should stabilize if self-similar):")
for n in range(1,8):
    lo,hi=F(1,n+1),F(1,n)
    inside=F(0)
    for a,b in ivs:
        aa=max(a,lo); bb=min(b,hi)
        if bb>aa: inside+=bb-aa
    frac = inside/(hi-lo)
    print(f"     window (1/{n+1},1/{n}) len={float(hi-lo):.4f}: cover-fraction={float(frac):.5f}")

print()
print("="*78)
print("PART 2 (H2): is consec a RENORMALIZATION fixed/periodic point?")
print("="*78)
# Renormalization map on offset sets: 'slope-renormalize' E.  Two natural candidate renorm maps:
#  (R1) DILATION (the inverse of compression): E -> spread.  consec is the unique FIXED set of
#       the COMPRESSION map (push all gaps to 1). Test: compression is a contraction whose unique
#       fixed point is consec, and measS7 increases monotonically along compression orbits.
#  (R2) Multiply offsets by a unit mod 7 and reduce (the Z/7 'twist' = cat-map-like linear map on
#       the torus).  Does measS7 have a symmetry/orbit under E -> (u*E mod something)?
def compress_step(E):
    """One compression step: reduce each gap>1 by 1 (toward consec). Fixed point = consec."""
    E=sorted(set(E)); out=[E[0]]
    for i in range(1,len(E)):
        g=E[i]-E[i-1]
        out.append(out[-1] + (g-1 if g>1 else 1))
    # renormalize min to 0
    m=out[0]
    return tuple(x-m for x in out)

print("  (R1) compression orbit toward consec; track measS7 (should rise to the consec fixed pt)")
for seed in [(0,1,2,3,4,5,6,9), (0,2,4,6,8,10,12,14), (0,1,2,3,4,5,6,7,12)]:
    k=len(seed)
    cur=tuple(seed); seen=set()
    chain=[]
    while cur not in seen:
        seen.add(cur)
        chain.append((cur, measS7(cur)))
        nxt=compress_step(cur)
        if nxt==cur: break
        cur=nxt
    consec=tuple(range(k))
    fixed = (chain[-1][0]==consec)
    mono = all(chain[i][1] <= chain[i+1][1] for i in range(len(chain)-1))
    print(f"   seed {seed}: orbit len={len(chain)}, ends at consec? {fixed}, measS7 monotone-up? {mono}")
    print(f"      measS7 chain: {[float(v) for _,v in chain]}")

# (R2): the Z/7 'cat-map' twist E -> u*E mod (something). Test scale-invariance orbit (THM-531).
print("\n  (R2) scale orbit E -> c*E (the dilation that PRESERVES measS7 -- scale invariance):")
for E in [(0,1,2,3,4,5,6,7),(0,2,4,6)]:
    vals=set()
    for c in range(1,6):
        Ec=tuple(c*e for e in E)
        vals.add(measS7(Ec))
    print(f"   E={E}: measS7 over scalings c=1..5: {sorted(float(v) for v in vals)}  invariant? {len(vals)==1}")

print()
print("="*78)
print("PART 3 (H3): Lyapunov/contraction -- does spread -> iid surjection at the Gauss rate?")
print("="*78)
# iid baseline: a random map e->uniform Z/7 surjects with prob p_iid(k) = 7! S(k,7)/7^k.
def stirling2(n,k):
    # exact
    return sum((-1)**(k-j)*factorial(k)//(factorial(j)*factorial(k-j))*(j**n) for j in range(k+1))//factorial(k)
def p_iid(k):
    return F(factorial(7)*stirling2(k,7), 7**k)
print("  iid surjection baseline p_iid(k)=7!*S(k,7)/7^k vs measS7(consec) and a WIDE spread:")
lam_gauss = pi*pi/(6*log(2))   # Lyapunov exponent of the Gauss map (Levy)
print(f"  (Gauss-map Lyapunov exponent lambda = pi^2/(6 ln2) = {lam_gauss:.5f})")
for k in [7,8,9,10]:
    consec=tuple(range(k))
    mc=measS7(consec)
    pid=p_iid(k)
    # a 'wide' spread: stretch by a factor to push slope orbit toward equidistribution
    wide=tuple(range(0, 3*k, 3))  # gap-3 spread, same k, span 3(k-1)
    mw=measS7(wide)
    # a primitive wide set spanning ~ a lot
    print(f"  k={k}: measS7(consec)={float(mc):.5f}  measS7(gap3 wide)={float(mw):.5f}  p_iid={float(pid):.5f}")
    print(f"        consec/p_iid ratio={float(mc/pid):.4f}  wide/p_iid={float(mw/pid):.4f}  (wide should -> 1)")

print("\n  Decay of |measS7(spread_t) - p_iid| as we widen (slope orbit equidistributes):")
k=8; pid=p_iid(k)
print(f"  k={k}, p_iid={float(pid):.5f}")
prev=None
for t in range(1,9):
    Et=tuple(range(0, t*k, t))  # gap-t spread
    if not primitive(Et):
        # still fine for measS7 magnitude; report
        pass
    m=measS7(Et)
    d=abs(float(m)-float(pid))
    ratio = (d/prev) if prev else float('nan')
    print(f"   gap={t} span={span(Et)}: measS7={float(m):.5f}  |meas-p_iid|={d:.5f}  ratio_to_prev={ratio:.4f}")
    prev=d

print("\nDONE.")
