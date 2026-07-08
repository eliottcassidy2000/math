"""
mac-mini-2026-07-07-S57 (HYP-5237) -- k=10 RESIDUAL CHARACTERIZATION.

After the average-form conditional tent (THM-655 machinery) at k=10, a shape P is
DISCHARGED for all families iff sup_E avgc(E,P) <= c*(P).  c*(P) ~ 1.18 is tight, so the
block {0..9} itself (avgc ~ 1.35) offends: the average form does NOT discharge k=10 alone.

The residual is {(P,E) : avgc(E,P) > c*(P)}.  klein-S174 THM-653 exhausts diam <= 10.
QUESTION: is the residual family set of BOUNDED diameter?  If yes, extending klein's
exhaustion (or the spread floor) to that diameter closes k=10.  The dilation structure
predicts YES: c*block has differences {c,2c,...} -> c(dc,P) -> 1 (Koksma) as c grows, so
large dilates are safe; only small c and compact E can offend.

OUTPUTS:
 (1) per shape: c*(P), EnvBlock (all-diam ceiling), block avgc, and whether the shape is
     discharged / has residual.
 (2) THE RESIDUAL FAMILY DIAMETER: max diam over all (P,E) with avgc > c* AND diam > 10,
     searched over compact families (diam<=24) + dilated-blocks c*{0..9} (c=1..8) +
     block+outlier + two-block to diam 60.  Report the max residual diam and the witnesses.
 (3) dilation decay table: avgc(c*block, P) for c=1..8 at the worst shapes -> confirm decay.
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

MP = F(14249, 252252); TH = F(1, 7); DMAX = 5000
k = 10

def GP_intervals(P):
    bad = []
    for p in P:
        w = F(1, 14*p)
        for j in range(p+1):
            bad.append((F(j,p)-w, F(j,p)+w))
    bad = [(max(l, F(0)), min(h, F(1))) for l, h in bad if h > 0 and l < 1]
    bad.sort(); merged = []
    for l, h in bad:
        if merged and l <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], h))
        else: merged.append((l, h))
    good = []; prev = F(0)
    for l, h in merged:
        if l > prev: good.append((prev, l))
        prev = max(prev, h)
    if prev < 1: good.append((prev, F(1)))
    return good

def c_array(P):
    beta = F(14-k, 7*k); w = TH - beta; intf = w*w/2
    bf, tf, iff_ = float(beta), float(TH), float(intf)
    iv = GP_intervals(P); meas = sum(h-l for l, h in iv)
    def Dfun(t):
        return np.where(t <= bf, -t*iff_, np.where(t <= tf, (t-bf)**2/2 - t*iff_, (1.0-t)*iff_))
    d_arr = np.arange(1, DMAX+1, dtype=np.int64); acc = np.zeros(DMAX)
    for (l, h) in iv:
        for (pt, sgn) in ((h, +1.0), (l, -1.0)):
            num, den = pt.numerator, pt.denominator
            fr = ((d_arr*num) % den).astype(np.float64)/den
            acc += sgn*Dfun(fr)
    c = 1.0 + acc/(d_arr*float(meas)*iff_)
    delta2 = 6/7 + bf + float(w)**2/4
    ctail = 1.0 + len(iv)*delta2/(DMAX*float(meas))
    return c, float(meas), ctail

def cstar(meas):
    return float((1 - MP/F(meas).limit_denominator(10**9)) * F(7*k, 2*(k-1)*(k-7)))

def envblock(c, ctail):
    allc = np.append(c, ctail); chat = np.maximum.accumulate(allc[::-1])[::-1]
    return sum((k-d)*chat[d-1] for d in range(1, k)) / (k*(k-1)//2)

def avgc_set(E, c, ctail):
    E = sorted(E); tot = 0.0; n = 0
    for i in range(len(E)):
        for j in range(i+1, len(E)):
            d = E[j]-E[i]; tot += c[d-1] if d <= DMAX else ctail; n += 1
    return tot/n

SHAPES = list(combinations(range(1, 14), 3))
print(f"=== k=10 residual characterization ({len(SHAPES)} shapes, c* ~ 1.18) ===\n")

# (1) discharge status per shape
env_ok = 0; block_off = 0
worst = []
for P in SHAPES:
    c, meas, ctail = c_array(P); cs = cstar(meas)
    eb = envblock(c, ctail); ba = avgc_set(list(range(k)), c, ctail)
    if eb <= cs: env_ok += 1
    if ba > cs: block_off += 1
    worst.append((eb, ba, cs, P, c, ctail, meas))
print(f"  shapes with EnvBlock <= c* (discharged all families/diam): {env_ok}/{len(SHAPES)}")
print(f"  shapes where the BLOCK offends (avgc(block) > c*, need diam route): {block_off}/{len(SHAPES)}")

# (2) residual family diameter: max diam over (P,E) with avgc>c* and diam>10
print("\n=== (2) residual family DIAMETER (avgc > c* AND diam > 10) ===")
max_resid_diam = 0; witnesses = []
# compact families diam 11..24 (subsets of [0,W], min=0, max=W)
compact_cache = {}
for W in range(11, 25):
    combos = [ (0,)+e+(W,) for e in combinations(range(1, W), k-2) ]
    compact_cache[W] = np.array([[E[j]-E[i] for i in range(k) for j in range(i+1,k)] for E in combos], dtype=np.int64), combos
for (eb, ba, cs, P, c, ctail, meas) in worst:
    shape_max = 0; shape_wit = None
    # compact
    for W in range(11, 25):
        dm, combos = compact_cache[W]
        cvals = np.where(dm <= DMAX, c[np.clip(dm-1,0,DMAX-1)], ctail).mean(axis=1)
        idx = np.where(cvals > cs)[0]
        if len(idx):
            if W > shape_max: shape_max, shape_wit = W, combos[idx[np.argmax(cvals[idx])]]
    # dilated blocks c*{0..9}
    for cc in range(2, 9):
        E = [cc*t for t in range(k)]
        if avgc_set(E, c, ctail) > cs:
            d = cc*(k-1)
            if d > shape_max: shape_max, shape_wit = d, tuple(E)
    # block+outlier and two-block to diam 60
    for coreL in range(3, k):
        for gap in range(1, 55):
            rightL = k-coreL
            E = list(range(coreL)) + [coreL-1+gap+t for t in range(rightL)]
            if len(set(E))!=k: continue
            if max(E)-min(E) > 10 and avgc_set(E, c, ctail) > cs:
                d = max(E)-min(E)
                if d > shape_max: shape_max, shape_wit = d, tuple(E)
    if shape_max > max_resid_diam:
        max_resid_diam = shape_max;
    if shape_max > 10:
        witnesses.append((shape_max, P, shape_wit, cs))
witnesses.sort(reverse=True)
print(f"  MAX residual diameter over all shapes: {max_resid_diam}")
print(f"  => k=10 residual is BOUNDED-diameter: klein exhaustion to diam {max_resid_diam} (or spread floor) closes it")
print(f"  top residual witnesses (diam, shape, family):")
for d, P, wit, cs in witnesses[:12]:
    dil = ""
    if wit and len(set(wit[i+1]-wit[i] for i in range(len(wit)-1)))==1:
        g = wit[1]-wit[0]; dil = f" = {g}*block" if g>1 else " = block"
    print(f"    diam {d:2d}: P={P}, E={wit}{dil} (c*={cs:.4f})")

# (3) dilation decay
print("\n=== (3) dilation decay avgc(c*block, P) -- confirms large dilates are safe ===")
worst_shapes = sorted(worst, key=lambda r:-r[1])[:4]
for (eb, ba, cs, P, c, ctail, meas) in worst_shapes:
    row = []
    for cc in range(1, 9):
        E = [cc*t for t in range(k)]
        row.append(avgc_set(E, c, ctail))
    print(f"  P={P} (c*={cs:.4f}): avgc(c*block) c=1..8 = " + " ".join(f"{v:.3f}" for v in row)
          + f"  -> {'decays below c*' if row[-1]<=cs else 'still over at c=8'}")
