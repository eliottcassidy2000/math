"""
EXACT E[maxgap] for LRC(14) families, and whether the AP minimizes it.
boxeph-2026-07-07-S1, HYP-4760.

E[maxgap](E) = integral_0^1 (circular max-gap of {frac(e x): e in E}) dx, EXACT.
On each order-interval every gap is linear in x, so maxgap is piecewise-linear
with breakpoints where two gaps cross OR the sorted order changes.  We integrate
the upper envelope exactly by subdividing at gap-crossings (also rational).

Key question: kps-S58 named the crux 'AP-minimality of E[maxgap]'.  The sampled
data suggests GW={1..11,13,24} has E[maxgap] BELOW the AP.  Settle it exactly.
"""
from fractions import Fraction as F

def order_breakpoints(E):
    """x-values in (0,1) where the circular sorted order of {frac(e x)} changes:
    e_i x - e_j x in Z  => x = m/(e_i - e_j).  Plus wrap coincidences e_i x in Z."""
    bps = set()
    n = len(E)
    for i in range(n):
        if E[i] != 0:
            for m in range(0, abs(E[i])+1):
                bps.add(F(m, abs(E[i])))
        for j in range(i+1, n):
            d = abs(E[i]-E[j])
            if d == 0: continue
            for m in range(0, d+1):
                bps.add(F(m, d))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    return bps

def gaps_linear(E, mid, a, b):
    """On (a,b) with sample mid, return list of (c,b0) so gap_s(x)=c*x+b0 >=0,
    consecutive circular gaps of {frac(e x)}. Order fixed on (a,b)."""
    fl = {e: (e*mid).__floor__() for e in E}
    order = sorted(E, key=lambda e: e*mid - fl[e])
    m = len(order)
    gaps = []
    for s in range(m):
        e1 = order[s]; e2 = order[(s+1) % m]
        if s < m-1:
            c = F(e2-e1); b0 = F(-(fl[e2]-fl[e1]))
        else:
            c = F(order[0]-order[-1]); b0 = F(-(fl[order[0]]-fl[order[-1]]) + 1)
        gaps.append((c, b0))
    return gaps

def integral_max_linear(gaps, a, b):
    """integral_a^b max_s (c_s x + b0_s) dx, exact.  Subdivide at pairwise
    crossings of the linear functions inside (a,b), then on each sub-interval
    the max is a single line; sum the trapezoids."""
    # crossing x of line i and j: c_i x + b0_i = c_j x + b0_j
    xs = {a, b}
    n = len(gaps)
    for i in range(n):
        ci, bi = gaps[i]
        for j in range(i+1, n):
            cj, bj = gaps[j]
            if ci == cj: continue
            x = (bj - bi)/(ci - cj)
            if a < x < b: xs.add(x)
    xs = sorted(xs)
    total = F(0)
    for lo, hi in zip(xs, xs[1:]):
        if lo == hi: continue
        mid = (lo+hi)/2
        # value of max at mid
        best = max(c*mid + b0 for c, b0 in gaps)
        # which line attains it (unique generically); integrate that line lo..hi
        # find argmax line
        bc, bb = max(gaps, key=lambda cb: cb[0]*mid + cb[1])
        # integral of bc*x+bb over [lo,hi]
        total += bc*(hi*hi - lo*lo)/2 + bb*(hi-lo)
    return total

def E_maxgap_exact(E):
    E = list(E)
    bps = order_breakpoints(E)
    total = F(0)
    for a, b in zip(bps, bps[1:]):
        if a == b: continue
        mid = (a+b)/2
        gaps = gaps_linear(E, mid, a, b)
        total += integral_max_linear(gaps, a, b)
    return total

# ---------------- run ----------------
families = {
    "AP {1..13}":            list(range(1,14)),
    "GW {1..11,13,24}":      [1,2,3,4,5,6,7,8,9,10,11,13,24],
    "13-ladder{1..12,26}":   [1,2,3,4,5,6,7,8,9,10,11,12,26],  # {1..12,13k}, k=2
    "10-ladder{1..9,11,12,13,20}": [1,2,3,4,5,6,7,8,9,11,12,13,20],
    "2*AP {2..26}":          [2*j for j in range(1,14)],
    "1..6,20..26":           [1,2,3,4,5,6,20,21,22,23,24,25,26],
}
print("EXACT E[maxgap] (k=13 families):")
vals = {}
for nm, E in families.items():
    v = E_maxgap_exact(E)
    vals[nm] = v
    print(f"  {nm:32s} E[maxgap] = {str(v):>18} = {float(v):.6f}")

ap = vals["AP {1..13}"]
print(f"\nAP value = {float(ap):.6f} = 1/7 + {float(ap-F(1,7)):.6f}  (1/7={1/7:.6f})")
print("Below AP?  " + ", ".join(f"{nm}:{'YES' if v<ap else 'no'}" for nm,v in vals.items() if nm!="AP {1..13}"))

# ---- search near-AP one/two-swap ladders for the TRUE minimizer of E[maxgap] ----
print("\n--- one-swap ladder scan  {1..13}\\{j} u {j*k}  (E[maxgap] exact) ---")
best = (ap, "AP")
for j in range(1,14):
    for k in range(2,6):
        E = [x for x in range(1,14) if x!=j] + [j*k]
        if len(set(E))<13: continue
        v = E_maxgap_exact(sorted(E))
        if v < best[0]:
            best = (v, f"{{1..13}}\\{{{j}}} u {{{j*k}}}={sorted(E)}")
    # brief print for a few
print(f"  best one-swap E[maxgap] = {float(best[0]):.6f} at {best[1]}")

# two-swap around GW neighborhood (small)
print("\n--- some two-swap near GW ---")
cands = [
    [1,2,3,4,5,6,7,8,9,10,13,24,12],  # perturb
    [1,2,3,4,5,6,7,8,9,10,11,13,26],
    [1,2,3,4,5,6,7,8,9,11,13,22,24],
    [1,2,3,4,5,6,7,8,10,11,13,20,24],
]
for E in cands:
    E=sorted(set(E))
    if len(E)!=13: continue
    v=E_maxgap_exact(E)
    print(f"  {E}  E[maxgap]={float(v):.6f}{'  <BELOW AP' if v<ap else ''}{'  <BELOW GW' if v<vals['GW {1..11,13,24}'] else ''}")

print("\nSUMMARY: AP minimizes the TAIL mu_{1/7} (confirmed elsewhere) but does the")
print("MEAN E[maxgap]?  If GW<AP here, 'AP-minimality of E[maxgap]' (kps-S58) is FALSE.")
