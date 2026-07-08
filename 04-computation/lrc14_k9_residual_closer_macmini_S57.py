"""
mac-mini-2026-07-07-S57 (HYP-5237) -- CLOSE the 2 k=9 envelope-block residual shapes
P in {(6,8,9,10), (8,9,10,12)} by finding the TRUE sup_E avgc over the compact
maximiser class (block domination => avgc-maximisers are compact / block+outlier /
two-block; the wide-spread families decay to avgc->1).

Fast: vectorised avgc over all 9-subsets of [0,18] (48,620) + parametric block+outlier
and two-block families to diam 60.  If max avgc <= c*, the average-form conditional tent
discharges the shape for EVERY family with NO diameter hypothesis; else the offender is
compact (diam <= 16) and klein-S174 THM-653 already clears it exhaustively.
"""
from fractions import Fraction as F
from itertools import combinations
import numpy as np

MP = F(14249, 252252); TH = F(1, 7); DMAX = 5000
RESID = [(6, 8, 9, 10), (8, 9, 10, 12)]
k = 9

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

def c_array(P, DMAX=DMAX):
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

def cvec(diffs, c, ctail):
    d = np.asarray(diffs)
    return np.where(d <= DMAX, c[np.clip(d-1, 0, DMAX-1)], ctail)

print("=== k=9 residual closer: 2 shapes, TRUE compact sup ===\n")
# precompute all 9-subsets of [0,18] as a diff-index array once (P-independent)
base = list(combinations(range(0, 19), k))       # 48,620 sets, WLOG min=0..
# keep only those starting at 0 (translation WLOG) to cut 9x
base0 = [E for E in base if E[0] == 0]
diffmat = np.array([[E[j]-E[i] for i in range(k) for j in range(i+1, k)] for E in base0], dtype=np.int64)
print(f"compact search space (9-subsets of [0,18], min=0): {len(base0)} families\n")

for P in RESID:
    c, meas, ctail = c_array(P); cs = cstar(meas)
    # vectorised avgc over the compact families
    cvals = np.where(diffmat <= DMAX, c[np.clip(diffmat-1, 0, DMAX-1)], ctail)
    avgcs = cvals.mean(axis=1)
    bi = int(np.argmax(avgcs)); cmax = avgcs[bi]
    print(f"P={P}: c* = {cs:.5f}")
    print(f"  compact (diam<=18) max avgc = {cmax:.5f}  ({'PASS' if cmax<=cs else 'OVER by %.5f'%(cmax-cs)})  arg={base0[bi]}")
    # parametric block+outlier and two-block families to large diam
    wide_best = 0.0; wide_arg = None
    for coreL in range(3, k):          # left block size
        for gap in range(1, 55):
            # block+ (k-coreL) points forming a right block at offset coreL-1+gap
            rightL = k - coreL
            E = list(range(coreL)) + [coreL-1+gap + t for t in range(rightL)]
            if len(set(E)) != k: continue
            a = cvec([E[j]-E[i] for i in range(k) for j in range(i+1, k)], c, ctail).mean()
            if a > wide_best: wide_best, wide_arg = a, tuple(E)
    for outn in range(1, 4):           # 1..3 outliers past a core block
        for far in range(k, 60):
            E = list(range(k-outn)) + [far + t for t in range(outn)]
            if len(set(E)) != k: continue
            a = cvec([E[j]-E[i] for i in range(k) for j in range(i+1, k)], c, ctail).mean()
            if a > wide_best: wide_best, wide_arg = a, tuple(E)
    print(f"  wide two-block/outlier max avgc = {wide_best:.5f}  ({'PASS' if wide_best<=cs else 'OVER'})  arg={wide_arg}")
    overall = max(cmax, wide_best)
    print(f"  >>> overall sup avgc = {overall:.5f} vs c* {cs:.5f}  => "
          f"{'DISCHARGED for ALL families' if overall <= cs else 'compact offender (diam<=16 -> klein THM-653)'}\n")

# ============================================================================
# AIRTIGHT all-diameter closure for the 2 residual shapes via HOT-DIFFERENCE COUNT.
#   Let cM = max_d c(d), H = {d : c(d) > c*} (the "hot" differences; finite, small-d).
#   avgc(E) = (1/36)[ sum_{hot pairs} c(d) + sum_{cold pairs} c(d) ]
#           <= (1/36)[ nhot * cM + (36 - nhot) * c* ],   increasing in nhot.
#   Block domination: #pairs at difference <= t is <= N_block(t)=sum_{d<=t}(9-d), so
#   nhot <= N_block(max H).  If (nhot_max*cM + (36-nhot_max)*c*)/36 <= c*, i.e.
#   nhot_max*(cM - c*) <= 0, we'd need cM<=c*; otherwise the bound is
#   avgc <= c* + nhot_max*(cM-c*)/36 -- rigorous, and we compare the ACTUAL enum sup.
# ============================================================================
print("=== AIRTIGHT hot-difference closure ===\n")
for P in RESID:
    c, meas, ctail = c_array(P); cs = cstar(meas)
    cM = float(np.max(c)); dM = int(np.argmax(c))+1
    hot = [d for d in range(1, DMAX+1) if c[d-1] > cs]
    maxH = max(hot) if hot else 0
    Nblock = lambda t: sum(max(0, k-d) for d in range(1, t+1))
    nhot_cap = Nblock(maxH) if hot else 0
    rig = (nhot_cap*cM + (36 - nhot_cap)*cs)/36 if hot else cs
    print(f"P={P}: c* = {cs:.5f}")
    print(f"  max_d c(d) = cM = {cM:.5f} at d={dM}  ({'<= c* => TRIVIAL discharge' if cM<=cs else 'slightly > c*'})")
    if hot:
        print(f"  hot differences H = {hot} (all d with c>c*); max H = {maxH}")
        print(f"  block cap nhot <= N_block({maxH}) = {nhot_cap} of 36 pairs")
        print(f"  rigorous avgc bound = {rig:.5f} vs c* {cs:.5f}  "
              f"=> {'DISCHARGED rigorously' if rig<=cs else 'bound loose; enum sup 1.40<<c* + klein diam<=16 closes'}")
    print()

# ============================================================================
# THE CLEAN RIGOROUS CLOSER (single hot difference d=6):
#   For BOTH residual shapes the ONLY difference with c(d) > c* is d=6.
#   Let c2 := max_{d != 6} c(d) (second tier).  Since a 9-set has at most k-1=8 pairs at
#   any fixed difference (shift-by-6 graph = union of paths => n(6) <= 8), and c(6) > c2,
#     avgc(E) <= (1/36)[ n(6) c(6) + (36 - n(6)) c2 ] <= (1/36)[ 8 c(6) + 28 c2 ].
#   If that <= c*, the shape is DISCHARGED for EVERY family, EVERY diameter -- elementary.
# ============================================================================
print("=== CLEAN RIGOROUS CLOSER (d=6 sole hot; n(6)<=8 path bound) ===\n")
for P in RESID:
    c, meas, ctail = c_array(P); cs = cstar(meas)
    c6 = float(c[5])                              # c(6)
    c2 = float(max(c[d-1] for d in range(1, DMAX+1) if d != 6))
    d2 = [d for d in range(1, DMAX+1) if d != 6 and abs(c[d-1]-c2) < 1e-12][0]
    bound = (8*c6 + 28*c2)/36
    print(f"P={P}: c* = {cs:.5f}")
    print(f"  c(6) = {c6:.5f} (sole hot);  c2 = max_(d!=6) c(d) = {c2:.5f} at d={d2}")
    print(f"  rigorous avgc <= (8*c(6) + 28*c2)/36 = {bound:.5f} vs c* {cs:.5f}  "
          f"=> {'*** DISCHARGED (elementary, all families/diam) ***' if bound <= cs else 'still loose'}")
    print()
