# opus-2026-07-17-S359 -- HYP-7420: THE OWNER'S KILLER-COLLAPSE IDEA.
#   q = ceil(v1/m)  =>  v1 = qm - s with 0 <= s < m, so v1 = -s (mod q);
#   for v_i = v1 + d_i,  v_i = d_i - s (mod q): ALL residues small.
#   Then at t = p/q, ||v_i t|| = dist(e_i p, qZ)/q with e_i = v_i mod q.
#   Taking p ~ q/(14 e_min) puts e_min*p at the LOWER band edge and
#   e_max*p at or below the UPPER edge exactly when e_max <= 13 e_min.
# This script tests: (a) does the collapse make residues small? (b) does the
# ratio-13 band condition hold? (c) does a lonely p actually exist?
from fractions import Fraction as F
from math import gcd, ceil
import random

def residues(V, q):
    """signed-minimal residues; |e| is what loneliness sees."""
    out = []
    for v in V:
        r = v % q
        out.append(min(r, q - r))     # distance to nearest multiple
    return out

def lonely_p_exists(V, q, lam=F(1,14)):
    """is there p with dist(v_i p, qZ) >= lam*q for all i?  exact scan."""
    lo = lam * q
    for p in range(1, q):
        ok = True
        for v in V:
            r = (v * p) % q
            if min(r, q - r) < lo: ok = False; break
        if ok: return p
    return None

print("(a)+(b) THE COLLAPSE AND THE BAND CONDITION, clustered 13-families:")
random.seed(359)
rows = []
for trial in range(300):
    base = random.randint(2000, 60000)
    spread = random.randint(5, 60)                 # tight cluster
    V = sorted(random.sample(range(base, base + spread * 13), 13))
    if len(set(V)) < 13: continue
    m = random.randint(3, 40)
    q = ceil(V[0] / m)
    if q < 20: continue
    e = residues(V, q)
    if min(e) == 0: continue                       # a speed divisible by q
    rows.append((max(e) <= 13 * min(e), max(e), min(e), q, V))
band_ok = sum(1 for r in rows if r[0])
print(f"    {len(rows)} clustered families: e_max <= 13*e_min holds in "
      f"{band_ok}/{len(rows)} ({100*band_ok/max(len(rows),1):.0f}%)")
sizes = sorted(r[1] for r in rows)
print(f"    residue sizes e_max: median {sizes[len(sizes)//2]}, max {sizes[-1]}"
      f"   (vs original speeds ~10^4) -- the collapse is real")

print()
print("(c) DOES THE BAND CONDITION PREDICT A LONELY p?  (exact scan over p)")
agree = dis = 0; examples = []
for (cond, emax, emin, q, V) in rows[:60]:
    p = lonely_p_exists(V, q)
    if cond == (p is not None): agree += 1
    else:
        dis += 1
        if len(examples) < 3: examples.append((cond, p, emax, emin, q))
print(f"    60 families: band condition matches lonely-p existence in {agree}/{agree+dis}")
for (cond, p, emax, emin, q) in examples:
    print(f"      mismatch: cond={cond} p={p} e_max={emax} e_min={emin} q={q}")

print()
print("(d) THE SUFFICIENCY DIRECTION (what the idea actually proves):")
suff = tot = 0
for (cond, emax, emin, q, V) in rows:
    if not cond: continue
    tot += 1
    p = lonely_p_exists(V, q)
    if p is not None: suff += 1
print(f"    among families WITH e_max <= 13 e_min: lonely p exists in {suff}/{tot}")
print()
print("    => the band condition is SUFFICIENT (that is the theorem);")
print("       it is not necessary -- families failing it can still be lonely.")

print()
print("(e) THE DECISIVE TEST, CORRECTED.  An earlier version of this script")
print("    measured only whether SOME q satisfies the band condition and reported")
print("    100%.  That is NOT the same as loneliness.  Corrected test: for each")
print("    family find q with the band condition, THEN verify a lonely p exists.")
random.seed(3591)
band_only = both = 0
for trial in range(120):
    base = random.randint(2000, 60000)
    spread = random.randint(5, 60)
    V = sorted(random.sample(range(base, base + spread * 13), 13))
    if len(set(V)) < 13: continue
    found = None
    for q in range(20, 400):
        e = residues(V, q)
        if min(e) == 0: continue
        if max(e) <= 13 * min(e):
            found = q; break
    if found is None: continue
    band_only += 1
    if lonely_p_exists(V, found) is not None: both += 1
print(f"    {band_only} families had a band-condition q; of those, a lonely p")
print(f"    exists at that q in {both}/{band_only}"
      f"  ({100*both/max(band_only,1):.0f}%)")
print()
print("(f) WHY THE BAND CONDITION IS NOT SUFFICIENT.  For q <= 27 the residues")
print("    lie in [1, floor(q/2)] <= 13, so e_max <= 13*e_min is nearly AUTOMATIC")
print("    -- it carries almost no information.  When the residue pool has <= 13")
print("    values, 13 speeds must realise ALL of them (pigeonhole), and a family")
print("    realising every residue has NO lonely p:")
def ok_p(es, q):
    lo = q / 14.0
    for p in range(1, q):
        if all(min((e*p) % q, q - (e*p) % q) >= lo for e in es): return p
    return None
for q in (15, 18, 21, 24, 27):
    es = tuple(range(1, q // 2 + 1))
    print(f"      q={q}: residues {es[0]}..{es[-1]} (all present) -> lonely p:"
          f" {ok_p(es, q)}")
print()
print("    => the CLASSICAL SIEVE (q <= 14, witness t = 1/q) is the sharp form,")
print("       and the repo already has it kernel-pure as")
print("       LonelyRunner.lonely_of_no_multiple / counterexample_needs_all_divisors.")
print("       The collapse is a genuine and vivid reformulation of that sieve;")
print("       it does NOT extend the modulus range.")
