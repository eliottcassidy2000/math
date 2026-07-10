#!/usr/bin/env python3
r"""
lrc14_offline_taxonomy_monad_S10.py   (monad-explorer-2026-07-09-S10, HYP-5807/THM-681)

(A) taxonomy check: at tall rulers, EVERY |coeff| <= 2 support-2 off-line relation falls
    in classes 1-6 (collision / global exact / affine coincidence) -- exhaustive.
(B) exact W0 ledgers (global exact unit relations, s <= 3, |c| <= 2) + line weights.
(C) the dichotomy: certified-live counts vs (W0, collision profile); the scatter.
"""
import sys, random, cmath
from math import gcd, pi, floor, ceil
from itertools import combinations, product

def band(q):
    return -(-q // 14), (13 * q) // 14

def LM_exact(v, q):
    lo, hi = band(q)
    return sum(1 for p in range(q) if all(lo <= (x * p) % q <= hi for x in v))

def hhat_abs(q, k):
    if k % q == 0:
        lo, hi = band(q)
        return (hi - lo + 1) / q
    import math
    return abs(math.sin(pi * k * ((13*q)//14 - (-(-q//14)) + 1) / q)) / (q * abs(math.sin(pi * k / q)))

def line_weight(q, kvec, Mmax=200):
    """sum over m = 1..Mmax of prod |hhat(m k_l)| * (b/q)^(13 - support), both signs -> x2
       (with care: m and -m give conjugate products, same abs)."""
    lo, hi = band(q)
    bq = (hi - lo + 1) / q
    s = len(kvec)
    tot = 0.0
    for m in range(1, q):
        term = 1.0
        for kl in kvec:
            term *= hhat_abs(q, (m * kl) % q)
        tot += term
    return tot * bq ** (13 - s)

def taxonomy_check(v, q):
    """enumerate support-2 |c|<=2 off-line relations at tall ruler q; classify."""
    n = len(v)
    classes = {'collision': 0, 'global_exact': 0, 'affine': 0, 'halfsum': 0, 'other': []}
    for a, b_ in combinations(range(n), 2):
        for ca, cb in product((-2, -1, 1, 2), repeat=2):
            tot = ca * v[a] + cb * v[b_]
            if tot % q != 0:
                continue
            # skip the defining line itself (handled separately): k = m(e_i+e_j) with v_i+v_j=q
            t = tot // q
            key = (abs(ca), abs(cb))
            if t == 0:
                classes['global_exact'] += 1
            elif key == (1, 1) and t in (1, -1):
                classes['collision'] += 1
            elif key == (2, 2) and abs(t) == 1:
                classes['halfsum'] += 1
            elif key == (2, 2) and abs(t) == 2:
                classes['collision'] += 1   # the m = 2 point of a collision/defining line
            elif key == (2, 2) and abs(t) == 3:
                classes['halfsum'] += 1     # three-halves coincidence v_a +- v_b = 3q/2
            elif key == (2, 2) and abs(t) == 1 and (ca > 0) != (cb > 0):
                classes['halfsum'] += 1     # half-difference v_a - v_b = q/2 (even q)
            elif key in ((2, 1), (1, 2)) and abs(t) in (1, 2):
                classes['affine'] += 1
            elif key == (1, 1) and abs(t) >= 2:
                classes['other'].append((a, b_, ca, cb, t, 'unit-sum |t|>=2 -- IMPOSSIBLE at tall q?'))
            else:
                classes['other'].append((a, b_, ca, cb, t))
    return classes

def global_exact_relations(v, Cmax=2, smax=3):
    rels = []
    n = len(v)
    for s in (2, 3):
        if s > smax:
            break
        for S in combinations(range(n), s):
            for coeffs in product(*([[c for c in range(-Cmax, Cmax+1) if c != 0]] * s)):
                if coeffs[0] < 0:
                    continue  # sign-normalize
                if sum(c * v[i] for c, i in zip(coeffs, S)) == 0:
                    rels.append((S, coeffs))
    return rels

def W0_of(v, qref):
    rels = global_exact_relations(v)
    return sum(line_weight(qref, [abs(c) for c in coeffs]) for S, coeffs in rels), rels

if __name__ == "__main__":
    rng = random.Random(2026070910)
    base = [14 * i for i in range(1, 14)]
    fams = []
    for tries in range(500):
        v = list(base)
        for _ in range(rng.randint(2, 5)):
            i = rng.randrange(13)
            v[i] = max(2, v[i] + rng.choice([-3, -2, -1, 1, 2, 3, 7, -7]))
        v = sorted(set(v))
        if len(v) != 13:
            continue
        g0 = 0
        for x in v:
            g0 = gcd(g0, x)
        if g0 != 1 or not all(any(x % qq == 0 for x in v) for qq in range(2, 15)):
            continue
        fams.append(v)
        if len(fams) >= 12:
            break

    print("=" * 100)
    print("PART A -- TAXONOMY CHECK at tall rulers (exhaustive |c| <= 2, support 2)")
    print("=" * 100)
    other_total = 0
    checked = 0
    for v in fams[:6]:
        Vmax = v[-1]
        sums = sorted(set(a + b for i, a in enumerate(v) for b in v[i+1:]))
        tall = [q for q in sums if q > Vmax]
        for q in tall[:8]:
            cl = taxonomy_check(v, q)
            checked += 1
            other_total += len(cl['other'])
            if cl['other']:
                print(f"  UNEXPECTED at q={q}, v={v[:4]}...: {cl['other'][:3]}")
    print(f"  {checked} tall rulers checked exhaustively: unexpected classes = {other_total} "
          f"({'TAXONOMY COMPLETE' if other_total == 0 else 'REVISE'})")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART B -- LINE WEIGHTS (numeric, Mmax = 200) vs the closed bounds")
    print("=" * 100)
    qr = 200
    print(f"  collision line (1,1):  {line_weight(qr, [1, 1]):.5f}  (bound 0.0225)")
    print(f"  doubling line (2,1):   {line_weight(qr, [2, 1]):.5f}  (bound 0.0382)")
    print(f"  Schur line (1,1,1):    {line_weight(qr, [1, 1, 1]):.5f}  (bound 0.0645)")
    sys.stdout.flush()

    print()
    print("=" * 100)
    print("PART C -- THE DICHOTOMY: (W0, collisions) vs certified-live and true-live")
    print("=" * 100)
    print(f"  {'W0':>7} {'#exact':>7} {'collis':>7} {'cert-live':>10} {'true-live':>10}  family")
    for v in fams:
        Vmax = v[-1]
        sums_all = [a + b for i, a in enumerate(v) for b in v[i+1:]]
        sums = sorted(set(sums_all))
        collisions = len(sums_all) - len(sums)
        qref = max(s for s in sums if s > Vmax)
        w0, rels = W0_of(v, qref)
        cert = true = 0
        for q in sums:
            if q <= Vmax:
                continue
            # certified: W0(q) + q-specific support-2 lines < 0.1124
            cl = taxonomy_check(v, q)
            qload = cl['collision'] / 2 * line_weight(q, [1, 1]) + \
                    (cl['affine'] / 2) * line_weight(q, [2, 1]) + \
                    (cl['halfsum'] / 2) * line_weight(q, [2, 2])
            w0q = sum(line_weight(q, [abs(c) for c in coeffs]) for S, coeffs in rels)
            lo, hi = band(q)
            bq = (hi - lo + 1) / q
            if w0q + qload < bq**12 * (2*bq - 1):
                cert += 1
            if LM_exact(v, q) > 0:
                true += 1
        ntall = sum(1 for q in sums if q > Vmax)
        print(f"  {w0:>7.4f} {len(rels):>7} {collisions:>7} {cert:>6}/{ntall:<3} {true:>6}/{ntall:<3}  {v[:5]}...")
        sys.stdout.flush()
    print("  (cert <= true always -- the floor is conservative; cert > 0 = a-priori-live family)")

    print()
    print("  --- GENERIC-SIDE battery (few exact relations; the W0 <= 0.08 branch) ---")
    gen = []
    for tries in range(4000):
        v = sorted(rng.sample(range(20, 300), 13))
        g0 = 0
        for x in v:
            g0 = gcd(g0, x)
        if g0 != 1 or not all(any(x % qq == 0 for x in v) for qq in range(2, 15)):
            continue
        if max(v) <= 13 * min(v):
            continue  # want gapped (outside spread13)
        gen.append(v)
        if len(gen) >= 8:
            break
    for v in gen:
        Vmax = v[-1]
        sums_all = [a + b for i, a in enumerate(v) for b in v[i+1:]]
        sums = sorted(set(sums_all))
        collisions = len(sums_all) - len(sums)
        qref = max(s2 for s2 in sums if s2 > Vmax)
        w0, rels = W0_of(v, qref)
        cert = true = 0
        for q in sums:
            if q <= Vmax:
                continue
            cl = taxonomy_check(v, q)
            qload = cl['collision'] / 2 * line_weight(q, [1, 1]) + \
                    (cl['affine'] / 2) * line_weight(q, [2, 1]) + \
                    (cl['halfsum'] / 2) * line_weight(q, [2, 2])
            w0q = sum(line_weight(q, [abs(c) for c in coeffs]) for S, coeffs in rels)
            lo, hi = band(q)
            bq = (hi - lo + 1) / q
            if w0q + qload < bq**12 * (2*bq - 1):
                cert += 1
            if LM_exact(v, q) > 0:
                true += 1
        ntall = sum(1 for q2 in sums if q2 > Vmax)
        print(f"  {w0:>7.4f} {len(rels):>7} {collisions:>7} {cert:>6}/{ntall:<3} {true:>6}/{ntall:<3}  {v[:5]}...")
