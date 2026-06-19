#!/usr/bin/env python3
"""
lrc14_angleC_threegap_extremality_macmini_0618s7.py   (mac-mini-2026-06-18-S7)

ANGLE C — THREE-GAP / DILATION-ORBIT EXTREMALITY deployment test.

GOAL: pin down WHETHER the AP-extremality of meas(S7(E)) (the dangerous-row crux)
follows from a clean three-gap / gap-distribution statement, and produce the exact
three-gap structure of meas(S7(AP_k)).

Setup recap.  E = {0=e_1<e_2<...<e_k}, primitive (gcd=1). For x in [0,1) the k points
P(x) = {frac(e_i x)} sit on the circle. S7(E) = { x : every sector [j/7,(j+1)/7) (j=0..6)
contains at least one point of P(x) }.  meas(S7(E)) is the crux carrier; the claim is
  consecutive AP {0,1,...,k-1} MAXIMIZES meas(S7(E))  for the dangerous rows k=8..11.

THREE-GAP CONNECTION (the literature deployment, Angle C):
  For the AP E={0,1,...,k-1}, the orbit P(x) = {0, frac(x), frac(2x), ..., frac((k-1)x)}
  is EXACTLY a Steinhaus/three-distance configuration: the first k multiples of x.
  By the THREE-GAP THEOREM (Sos-Suranyi-Swierczkowski), P(x) has at most THREE distinct
  gaps, and they are explicit in the continued fraction of x. So "all 7 sectors hit" is a
  pure three-gap event for the AP, and meas(S7(AP_k)) should be computable / boundable from
  the three-gap law alone. For a GENERAL E, the orbit is a union of dilated APs => Liang's
  3d-distance theorem (a union of d APs has <= 3d gaps): MORE gaps possible, orbit LESS
  uniform => HARDER to fill all 7 sectors => smaller meas(S7).  THIS is the mechanism
  behind AP-extremality: AP = single three-gap orbit = most uniform = best sector cover.

This script:
  (A) recompute meas(S7(AP_k)) exactly and verify the "all 7 sectors hit" event is a
      three-gap event: confirm the orbit P(x) for the AP has <=3 distinct gaps at sample x.
  (B) TEST the three-gap-uniformity mechanism: does "max-gap of P(x) < 1/7" (a sufficient
      condition for ALL sectors hit, since a 1/7-net hits every sector) coincide with a
      three-gap-only event for the AP, and is meas{maxgap<1/7} = meas(N(E)) maximized by AP?
  (C) The SHARP question for the certificate: a set E with a NON-AP structure has its orbit
      = union of >1 dilated APs => Liang 3d gaps. Test: does increasing the number of
      "AP-components" (additive-dimension proxy) strictly REDUCE meas(S7)?  This is the
      monotone mechanism. We measure meas(S7) vs the number of 3-term APs / additive energy.
  (D) Decisive AP-extremality search per dangerous row k=8..11 over a primitive box,
      reporting whether ANY E beats AP and, if a clean three-gap inequality would suffice.
"""
import sys, itertools
from fractions import Fraction as F
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)

SEV = F(1,7)

# ---------- meas(S7): all 7 sectors hit (exact, breakpoint sweep) ----------
def measS7(E):
    E=sorted(set(E))
    bps=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): bps.add(F(m,7*e))
    bps=sorted(bps)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        secs=set()
        for e in E:
            y=(e*xm)%1
            secs.add(int(y*7))
        if len(secs)==7:
            total+=x1-x0
    return total

# ---------- gap structure of orbit P(x) for given E at given x ----------
def orbit_gaps(E,x):
    pts=sorted(set((F(e)*x)%1 for e in E))
    g=[]
    for i in range(len(pts)):
        a=pts[i]; b=pts[(i+1)%len(pts)] if i+1<len(pts) else pts[0]+1
        g.append(b-a)
    return g

def n_distinct_gaps(E,x):
    return len(set(orbit_gaps(E,x)))

# ---------- additive structure proxies ----------
def num_3aps(E):
    """count 3-term APs a,b,c with b-a=c-b inside E (a proxy for low additive dimension)."""
    S=set(E); c=0
    El=sorted(E)
    for a in El:
        for b in El:
            if b<=a: continue
            cc=2*b-a
            if cc in S: c+=1
    return c

def additive_energy(E):
    """E(A) = #{(a,b,c,d): a+b=c+d}."""
    from collections import Counter
    sums=Counter()
    El=list(E)
    for a in El:
        for b in El:
            sums[a+b]+=1
    return sum(v*v for v in sums.values())

def normalize_primitive(E):
    E=sorted(set(E)); E=[e-E[0] for e in E]
    g=0
    for e in E: g=gcd(g,e)
    if g>1: E=[e//g for e in E]
    return tuple(E)

print("="*92)
print("ANGLE C: three-gap / dilation-orbit extremality of meas(S7)  [mac-mini-2026-06-18-S7]")
print("="*92)

# ---------- (A) AP orbit is a three-gap configuration ----------
print("\n[A] AP orbit {0,x,2x,...,(k-1)x} is a THREE-GAP configuration (<=3 distinct gaps).")
for k in [8,9,10,11,13]:
    AP=list(range(k))
    # sample many x, record max #distinct gaps observed
    maxg=0
    import random
    random.seed(1)
    for _ in range(400):
        num=random.randint(1,500); den=random.randint(1,500)
        x=F(num,den)%1
        if x==0: continue
        d=n_distinct_gaps(AP,x)
        maxg=max(maxg,d)
    print(f"  k={k:2d}: AP orbit max #distinct gaps over 400 random rationals x = {maxg}  (<=3 expected)")

# ---------- (B) meas(S7(AP)) exact + the 1/7-net sufficient condition ----------
print("\n[B] meas(S7(AP_k)) exact, and meas(N)=meas{maxgap<=1/7} (the three-gap net event).")
print("    A 1/7-net (maxgap<=1/7) hits every sector => N subset S7. For the AP, maxgap is")
print("    a pure three-gap quantity (largest of the <=3 gaps).")
def measN_threegap(E):
    """measure of x in [0,1) where the orbit P(x) has max circular gap <= 1/7. exact."""
    E=sorted(set(E))
    diffs=set()
    for a in range(len(E)):
        for b in range(a+1,len(E)):
            diffs.add(E[b]-E[a])
    bps=set([F(0),F(1)])
    for d in diffs:
        for m in range(0,d+1): bps.add(F(m,d))
    bps=sorted(x for x in bps if 0<=x<=1)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        g=orbit_gaps(E,xm)
        if max(g)<=SEV: total+=x1-x0
    return total
for k in [8,9,10,11,12,13]:
    AP=list(range(k))
    s7=measS7(AP); nn=measN_threegap(AP)
    print(f"  k={k:2d}: meas(S7(AP))={float(s7):.6f}   meas(N=maxgap<=1/7)(AP)={float(nn):.6f}   "
          f"S7-N gap={float(s7-nn):.6f}")

# ---------- (C) monotone mechanism: more AP-components / higher energy => more cover? ----------
print("\n[C] Mechanism test: orbit of E = union of dilated APs (Liang 3d-distance).")
print("    Hypothesis: AP (1 component, 3 gaps) gives MAX cover; splitting into more")
print("    components (=> up to 3d gaps, less uniform) REDUCES meas(S7).")
print("    Probe families at k=8: AP, then progressively dissociated; report meas(S7), #3APs, energy.")
fam8 = {
  "AP {0..7}":            [0,1,2,3,4,5,6,7],
  "1 hole {0..6,8}":      [0,1,2,3,4,5,6,8],
  "2 blocks":             [0,1,2,3,40,41,42,43],
  "4 blocks":             [0,1,20,21,40,41,60,61],
  "Sidon-ish":            [0,1,3,7,12,20,30,44],
  "dissociated 2^i":      [0,1,2,4,8,16,32,64],
}
for name,E in fam8.items():
    En=normalize_primitive(E)
    s7=measS7(list(En))
    print(f"  {name:<22} meas(S7)={float(s7):.5f}  #3APs={num_3aps(En):3d}  energy={additive_energy(En):5d}  E={En}")

# ---------- (D) decisive AP-extremality over a primitive box, per dangerous row ----------
print("\n[D] AP-extremality (decisive small-box search) for dangerous rows k=8..11.")
print("    Search primitive normalized k-sets with max element <= B; does ANY beat AP?")
def search_row(k, B):
    AP=tuple(range(k))
    apval=measS7(list(AP))
    best=apval; bestE=AP; beaters=0; total=0
    # iterate subsets of {1..B} of size k-1 (0 always in), normalized primitive
    seen=set()
    for combo in itertools.combinations(range(1,B+1), k-1):
        E=(0,)+combo
        En=normalize_primitive(E)
        if En in seen: continue
        seen.add(En); total+=1
        v=measS7(list(En))
        if v>apval: beaters+=1
        if v>best: best=v; bestE=En
    return apval,best,bestE,beaters,total
for k,B in [(8,12),(9,12),(10,11),(11,11)]:
    apval,best,bestE,beaters,total=search_row(k,B)
    flag = "AP IS MAX" if best==apval else f"AP BEATEN by {bestE}"
    print(f"  k={k:2d} (box maxE<={B}, {total} primitive shapes): "
          f"meas(S7(AP))={float(apval):.6f}  best={float(best):.6f}  beaters={beaters}  => {flag}")

print("\n[E] Three-gap closed form attempt for meas(S7(AP_k)): is it a rational with")
print("    denominator dividing 7*lcm? Report exact fraction + denominator factorization.")
for k in [8,9,10,11,12,13]:
    AP=list(range(k))
    s7=measS7(AP)
    print(f"  k={k:2d}: meas(S7(AP)) = {s7}   (= {float(s7):.6f})  denom={s7.denominator}")

print("\nDONE.")
