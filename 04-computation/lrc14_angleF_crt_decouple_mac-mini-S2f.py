#!/usr/bin/env python3
"""
ANGLE F part 2 — the CRT DECOUPLING and the precise solvability obstruction.

We make the CRT structure explicit. The level is 1/14 with 14 = 2 * 7.
A witness tau must satisfy ||v tau|| >= 1/14 for every v in S.

CRT DECOMPOSITION OF THE PROBLEM.
We model the witness as tau = a / q with q chosen. The crux: which q?
 - The "covering" of S is exactly that every modulus 2..14 divides some speed.
   This is what makes the residue problem HARD: it ties the speeds to the moduli
   that define the safe band.

THIS SCRIPT establishes the following, all EXACT:

(A) The 7-adic + 2-adic SPLIT of the safe band. Show DANGER_q for q with 7|q
    and 2|q is governed mod 7 and mod 2 separately ONLY in the limit; for finite q
    the band 14*min(r,q-r)>=q couples them.

(B) The DENOMINATOR LADDER. For each in-scope S, compute the EXACT optimal tau = a/D
    and its denominator D. Tabulate D over many S3 sets: is D bounded? what primes
    appear? This tells us the modulus the CRT witness needs.

(C) The CORE CRT THEOREM attempt: if we fix tau = a/(14*m) where m = Vmax-scale,
    decompose a mod 14 (controls small speeds dividing 14) and a mod m (controls cluster).
    Determine when the two conditions are JOINTLY solvable.

(D) THE SOLVABILITY OBSTRUCTION. Characterize the S for which NO rational witness
    with denominator <= some bound exists, isolating the "irreducibly Diophantine" core.

(E) THE COVERING-SYSTEM REFORMULATION. Each speed v forbids a set A_bad(v) of residues
    a mod q. Safe witness exists iff Union A_bad(v) != Z/q. Compute the union exactly,
    and check whether it's a COVERING SYSTEM (covers all residues) for the natural q.

All exact (Fractions / integer residues). Distinguish PROVED vs VERIFIED.
"""
from fractions import Fraction as F
from math import gcd
from itertools import combinations

def is_safe_residue(r, q):
    d = min(r % q, (q - r) % q)
    return 14 * d >= q

def circ_norm(x):
    x = x - int(x)
    if x < 0: x += 1
    return min(x, 1 - x)

def is_q_covering(S):
    return all(any(u % q == 0 for u in S) for q in range(2,15))

def is_primitive(S):
    g = 0
    for u in S: g = gcd(g,u)
    return g == 1

def M_exact_and_argmax(speeds):
    """Exact M and an argmax tau = a/D. Build candidate breakpoints exactly."""
    cands = set([F(0)])
    for v in speeds:
        for n in range(v+1):
            for off in (F(1,14), F(-1,14), F(0), F(1,2), F(13,14), F(-13,14)):
                t = (F(n)+off)/v
                t -= int(t)
                if t < 0: t += 1
                cands.add(F(t))
    best = F(0); arg = F(0)
    for t in cands:
        m = min(circ_norm(v*t) for v in speeds)
        if m > best:
            best = m; arg = t
    return best, arg

def header(t):
    print("\n"+"="*72); print(t); print("="*72)

# ----------------------------------------------------------------------
header("PART A — the 7-adic / 2-adic structure of the safe band")
print("""
Level 1/14, 14 = 2*7. ||v tau|| >= 1/14 means v tau is at least 1/14 from an integer.
For tau = a/q, danger = {r : 14*min(r,q-r) < q}. The 'width' of danger is 2*ceil(q/14)-1
residues centered at 0. The fraction of safe residues -> 13/14 = 1 - 1/14 - ... as q->inf,
but the band is a SINGLE interval mod q (around 0), NOT a product mod 2 x mod 7.
This is why naive CRT mod 2 x mod 7 does NOT capture the level: the level couples them.
""")
for q in (14, 70, 98, 7*13, 7*14, 2*7*13):
    danger = [r for r in range(q) if not is_safe_residue(r,q)]
    print(f"  q={q:4d}: danger band = {danger}  (width {len(danger)} = 2*ceil(q/14)-1={2*((q+13)//14)-1})")

# ----------------------------------------------------------------------
header("PART B — the DENOMINATOR LADDER: where does the optimum tau live?")
print("""
For each in-scope (primitive, q-covering, S3) set, the exact M is attained at a
rational tau = a/D. We collect D and its prime factorization. If the witness is
constructible by CRT with modulus q, then q must be a multiple of D (or D itself).
We sweep a family of S3 sets and tabulate.
""")

def prime_factors(n):
    f = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            f[d] = f.get(d,0)+1; n//=d
        d += 1
    if n>1: f[n] = f.get(n,0)+1
    return f

# Generate S3 sets: P subset of {1..13} (covering-ish) + a cluster L all >13, k>=3.
# We build a representative sweep.
sweep = []
base_P_options = [
    [1,2,3,4,5,6,7,8,9,10,11,12,13],     # full small part (too big alone for 13-set; we trim)
]
# Construct genuine 13-sets: choose P (size 13-k) and cluster L (size k).
import random
random.seed(1)
found_sets = []
attempts = 0
for k in (3,4,5,6,7):
    cnt = 0
    while cnt < 8 and attempts < 200000:
        attempts += 1
        psize = 13 - k
        if psize < 1: break
        P = sorted(random.sample(range(1,14), psize))
        # cluster: k speeds > 13, bounded spread (<=30) to stay near extremal regime
        base = random.randint(14, 60)
        spread = random.randint(k-1, 30)
        Lpool = list(range(base, base+spread+1))
        if len(Lpool) < k: continue
        L = sorted(random.sample(Lpool, k))
        S = sorted(set(P+L))
        if len(S) != 13: continue
        if not is_primitive(S): continue
        if not is_q_covering(S): continue
        found_sets.append((k,S))
        cnt += 1

print(f"  generated {len(found_sets)} primitive q-covering S3 sets (bounded-spread clusters)")
denom_primes = {}
maxD = 0
all_geq = True
for k,S in found_sets:
    M, arg = M_exact_and_argmax(S)
    D = arg.denominator
    maxD = max(maxD, D)
    for p in prime_factors(D):
        denom_primes[p] = denom_primes.get(p,0)+1
    if M < F(1,14):
        all_geq = False
        print(f"  *** COUNTEREXAMPLE? k={k} S={S} M={M} < 1/14 ***")
print(f"  ALL sampled M >= 1/14 : {all_geq}")
print(f"  max optimum-denominator D over sample = {maxD}")
print(f"  primes appearing in optimum denominators (count): "
      f"{dict(sorted(denom_primes.items()))}")
print("""  INTERPRETATION: the optimum tau lives at denominators built from the BINDING
  PAIR difference/sum (a Diophantine denominator); these are NOT bounded by a fixed
  small q -- D grows with the cluster scale. So a SINGLE fixed CRT modulus cannot
  serve all S3 sets. The modulus must scale with Vmax.""")

# ----------------------------------------------------------------------
header("PART C — does a fixed-modulus CRT witness exist? (q independent of S)")
print("""
TEST: is there a FIXED q (small) such that EVERY in-scope S has a safe rational
witness a/q? If yes -> constructive uniform witness (closes S3). If no -> obstruction.
We test the natural candidate moduli and report the FRACTION of sets covered.
""")
cand_q = [14,28,42,56,70,84,98,126,168,182,252,7*13,7*7,2*3*5*7]
for q in cand_q:
    nok = 0
    fails = []
    for k,S in found_sets:
        ok = any(all(is_safe_residue((v*a)%q,q) for v in S) for a in range(q))
        if ok: nok += 1
        else: fails.append((k,S))
    print(f"  q={q:4d}: {nok}/{len(found_sets)} sets have a safe rational witness a/q"
          + (f"   (sample fail: k={fails[0][0]} {fails[0][1]})" if fails else "  <-- COVERS ALL"))
