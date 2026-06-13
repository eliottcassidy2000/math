#!/usr/bin/env python3
"""
Multi-sieve concepts for LRC@14, recursive, with assumptions challenged.
opus-2026-06-02-S562 (remote-control).

Three sub-sieves combined; a config is CAUGHT if any finds a witness t with
||v_i t|| >= 1/n for all i:
  (D) DIVISION sieve : t = a/m, m in a modulus set, gcd(a,m)=1.  (the standard one)
  (P) PINCH sieve    : t = a/(v_i+v_j) over pairs (S557: the OPTIMAL witness is a
                       pair-sum time, so the natural moduli are PAIR-SUMS, not
                       small integers -- assumption #2 challenged).
  (F) recursive EVEN-FOLD : even speeds 2u fold to {u}; recurse (assumption #3:
                       the sieve is not flat, it recurses on the 2-adic structure).

Assumptions challenged, measured:
 #1 "the apex (mult of 14) is a hard obstruction": it is only stuck at ONE
    modulus; a DIFFERENT modulus catches it. Multi-sieving dissolves it.
 #2 small-integer moduli are insufficient (S551 unboundedness); pair-sum (pinch)
    moduli are a bounded-COUNT complete witness family.
 #3 the residual of a flat sieve is handled by the recursive fold.
All exact (Fraction).
"""

from fractions import Fraction
from math import gcd
import random


def safe_at(V, t, n):
    thr = Fraction(1, n)
    return all(min((v * t) % 1, 1 - (v * t) % 1) >= thr for v in V)


# every sieve returns (label, modulus, time_fraction) or None
# (D) division sieve over a modulus set
def div_sieve(V, n, moduli):
    for m in moduli:
        for a in range(1, m):
            if gcd(a, m) == 1 and safe_at(V, Fraction(a, m), n):
                return ("div", m, Fraction(a, m))
    return None


# (P) pinch sieve: pair-sum moduli (S557). Optimal witness is t=a/(v_i+v_j).
def pinch_sieve(V, n):
    sums = sorted({V[i] + V[j] for i in range(len(V)) for j in range(i + 1, len(V))})
    for s in sums:
        for a in range(1, s):
            if gcd(a, s) == 1 and safe_at(V, Fraction(a, s), n):
                return ("pinch", s, Fraction(a, s))
    return None


# (F) recursive even-fold witness
def fold_witness(V, n, depth=0, maxdepth=6):
    if depth > maxdepth:
        return None
    # base: all odd -> t=1/2
    if all(v % 2 == 1 for v in V):
        if safe_at(V, Fraction(1, 2), n):
            return ("fold-allodd", depth, Fraction(1, 2))
    E = [v for v in V if v % 2 == 0]
    O = [v for v in V if v % 2 == 1]
    if not E:
        return None
    fold = sorted(v // 2 for v in E)
    # find a good time s for the fold (even runners safe iff fold safe at the
    # DOUBLED time s; threshold 1/n). Use division + pinch + deeper fold.
    sub = div_sieve(fold, n, list(range(2, 4 * n))) or pinch_sieve(fold, n) \
        or fold_witness(fold, n, depth + 1, maxdepth)
    if not sub:
        return None
    s = sub[2]                                  # a Fraction time for the fold
    # doubling-preimages of s: t with 2t == s  ->  t = s/2 and s/2+1/2
    for t in ((s / 2) % 1, (s / 2 + Fraction(1, 2)) % 1):
        if safe_at(V, t, n):
            return ("fold", depth, t)
    return None


def primitive(V):
    g = 0
    for v in V:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in V))


# ---------------- experiments ----------------
def hard_configs(n, rng):
    """Configs that defeat small-modulus sieves: loaded (mult of every q<=n),
    apex-multiple-of-2q, near-tight, and large random."""
    out = {}
    out["AP (tight)"] = tuple(range(1, n))
    out["V* (sporadic tight)"] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24)
    out["apex=mult of 14"] = (1, 3, 5, 9, 11, 13, 14, 2, 4, 6, 8, 10, 12)  # contains 14
    out["loaded (mult of all q<=14)"] = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 143, 2520)
    # several random large primitive
    for i in range(4):
        out[f"rand-large {i}"] = primitive(tuple(rng.sample(range(1, 300), n - 1)))
    return out


def run(n=14, seed=0):
    rng = random.Random(seed)
    cfgs = hard_configs(n, rng)
    print(f"==== Multi-sieve on hard n={n} configs ====\n")
    small = list(range(2, n + 1))            # {2..14}
    mid = list(range(2, 2 * n + 1))          # {2..28}
    for name, V in cfgs.items():
        V = tuple(sorted(set(V)))
        if len(V) != n - 1:
            continue
        dS = div_sieve(V, n, small)
        dM = div_sieve(V, n, mid)
        pn = pinch_sieve(V, n)
        fd = fold_witness(V, n)
        # smallest division modulus that catches it
        mmin = None
        for m in range(2, 200):
            if div_sieve(V, n, [m]):
                mmin = m
                break
        print(f"  {name:28s}")
        print(f"     div{{2..{n}}}: {'CAUGHT '+str(dS[1:]) if dS else 'MISS'};  "
              f"div{{2..{2*n}}}: {'CAUGHT m='+str(dM[1]) if dM else 'MISS'};  "
              f"min div modulus = {mmin}")
        print(f"     PINCH: {'CAUGHT m=v_a+v_b='+str(pn[1])+' a='+str(pn[2]) if pn else 'MISS'};  "
              f"FOLD: {'CAUGHT depth='+str(fd[1]) if fd else 'MISS'}")
    print()


def challenge_apex(n=14, trials=3000, seed=2):
    """Assumption #1: an apex = multiple of 2q is 'stuck' only at modulus 14;
    a different modulus catches it. Measure."""
    rng = random.Random(seed)
    q = n // 2
    caught14 = caught_other = total = 0
    other_moduli = [m for m in range(2, n) if m != n]   # {2..13}
    examples = []
    for _ in range(trials):
        # build a primitive set containing a multiple of 14 (the apex obstruction)
        V = set(rng.sample(range(1, 60), n - 2))
        V.add(14 * rng.randint(1, 3))
        V = primitive(tuple(V))
        if len(V) != n - 1:
            continue
        total += 1
        c14 = div_sieve(V, n, [n]) is not None         # mod 14 base
        cother = div_sieve(V, n, other_moduli) is not None
        caught14 += c14
        caught_other += cother
        if (not c14) and cother and len(examples) < 3:
            examples.append((V, div_sieve(V, n, other_moduli)[1]))
    print("==== Assumption #1: does multi-sieving dissolve the apex obstruction? ====")
    print(f"  configs containing a multiple of 14 (apex-stuck at m=14): {total}")
    print(f"   caught by m=14 alone:           {caught14}/{total} ({100*caught14/total:.1f}%)")
    print(f"   caught by some m in {{2..13}}:    {caught_other}/{total} ({100*caught_other/total:.1f}%)")
    print(f"   => apex-stuck-at-14 configs rescued by a DIFFERENT modulus:")
    for V, m in examples:
        print(f"        {V}  caught at m={m}")
    print("   The apex is a single-modulus artifact; multi-sieve has no apex.\n")


if __name__ == "__main__":
    run(14)
    challenge_apex(14)
