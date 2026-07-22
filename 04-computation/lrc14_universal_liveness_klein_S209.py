#!/usr/bin/env python3
"""
lrc14_universal_liveness_klein_S209.py

HYP-5731(e) capstone: the UNIVERSAL SUPRA-VMAX LIVENESS law.

For q > Vmax the residues v_l mod q ARE the speeds (no reduction), so ruler
liveness is the clean modular statement:

    LIVE(S, q) :⟺ ∃ p ∈ Z_q : ∀ l, v_l·p mod q ∈ [⌈q/14⌉, q − ⌈q/14⌉].

Any live (q,p) certifies M(S) ≥ 1/14 at t = p/q directly (band ⟹ ‖v_l t‖ ≥
1/14; kps LRCPairSumDispatch is the Lean consumer). THM-668 says the max is
ATTAINED at pair-sum q; the S209 data suggests the SUPPLY is much larger:
every supra-Vmax q may work.

Tests over the S209 census generator (k≥8 covering, mid-band-heavy) +
adversarial instances:
 (1) PAIR-SUM q ∈ (Vmax, 2Vmax]: fraction live (expect 100%?).
 (2) ARBITRARY q ∈ (Vmax, 2Vmax] (all integers, not just pair-sums): live?
 (3) The minimum live DENSITY over everything (uniformity constant), and
     where it occurs; minimum over q per instance of LM(q).
 (4) Sanity boundary: q ≤ Vmax rulers can be dead (divisibility) — confirm
     the law needs q > Vmax.
Exact integer computation throughout.
"""

import random
from math import gcd

random.seed(31415)
QS = list(range(2, 15))


def is_covering(S):
    return all(any(s % q == 0 for s in S) for q in QS)


def LM(S, q):
    def r_safe(r):
        return 14 * r >= q and 14 * (q - r) >= q
    bad = bytearray(q)
    bad[0] = 1
    seen = set()
    for v in S:
        r = v % q
        key = min(r, (q - r) % q)
        if key in seen:
            continue
        seen.add(key)
        if r == 0:
            return 0
        g = gcd(r, q)
        rr, qq = r // g, q // g
        inv = pow(rr, -1, qq)
        for m in range(qq):
            s = m * g
            if not r_safe(s):
                p0 = (m * inv) % qq
                for t in range(g):
                    bad[p0 + t * qq] = 1
    return q - sum(bad)


def gen_instance(V):
    P = random.choice([(8, 9, 10, 12), (7, 9, 10, 11, 12), (11, 12, 13),
                       (10, 11, 12, 13), (9, 11, 13)])
    k = 13 - len(P)
    if k < 8:
        return None
    L = {V}
    missed = [q for q in QS if not any(p % q == 0 for p in P)]
    for q in missed:
        if any(u % q == 0 for u in L):
            continue
        lo, hi = -(-14 // q), V // q
        if lo > hi:
            return None
        L.add(q * random.randint(lo, hi))
    for _ in range(3):
        if len(L) < k:
            L.add(random.randint(max(14, V // 14 + 1), max(16, 9 * V // 14 - 1)))
    while len(L) < k:
        L.add(random.randint(14, V))
    S = sorted(set(P) | L)
    if len(S) == 13 and is_covering(S):
        return S
    return None


def main():
    print("=" * 78)
    print("UNIVERSAL SUPRA-VMAX LIVENESS test (klein-S209)")
    print("=" * 78)

    adversarial = [
        [12, 33, 46, 47, 68, 73, 79, 81, 85, 87, 91, 112, 120],
        [31, 33, 45, 48, 73, 76, 82, 86, 98, 102, 103, 104, 120],
        [62, 66, 69, 102, 109, 118, 120, 126, 130, 136, 159, 185, 200],
        [9, 16, 24, 33, 40, 47, 54, 62, 65, 70, 77, 84, 91],
    ]
    instances = list(adversarial)
    while len(instances) < 40:
        S = gen_instance(random.choice([120, 160, 200, 260]))
        if S is not None:
            instances.append(S)

    tot_ps = live_ps = 0
    tot_all = live_all = 0
    min_dens = 2.0
    min_at = None
    min_LM_pairsum_per_inst = []
    dead_sub = 0
    tot_sub = 0
    for S in instances:
        V = max(S)
        pair_qs = sorted({a + b for i, a in enumerate(S) for b in S[i:] if a + b > V})
        all_qs = list(range(V + 1, 2 * V + 1))
        inst_min = None
        for q in pair_qs:
            lm = LM(S, q)
            tot_ps += 1
            live_ps += 1 if lm > 0 else 0
            d = lm / (q - 1)
            if d < min_dens:
                min_dens, min_at = d, (tuple(S), q, 'pair-sum')
            inst_min = lm if inst_min is None or lm < inst_min else inst_min
        min_LM_pairsum_per_inst.append(inst_min)
        # arbitrary q: sample 40 per instance (cost control)
        for q in random.sample(all_qs, min(40, len(all_qs))):
            lm = LM(S, q)
            tot_all += 1
            live_all += 1 if lm > 0 else 0
            d = lm / (q - 1)
            if d < min_dens:
                min_dens, min_at = d, (tuple(S), q, 'arbitrary')
        # sub-Vmax sanity: sample 15
        for q in random.sample(range(15, V + 1), 15):
            tot_sub += 1
            if LM(S, q) == 0:
                dead_sub += 1

    print(f"\ninstances: {len(instances)} (4 adversarial/hard + census)")
    print(f"[1] PAIR-SUM q > Vmax:   live {live_ps}/{tot_ps}"
          f"  ({100*live_ps/tot_ps:.2f}%)")
    print(f"[2] ARBITRARY q in (Vmax, 2Vmax] (sampled): live {live_all}/{tot_all}"
          f"  ({100*live_all/tot_all:.2f}%)")
    print(f"[3] minimum live density observed: {min_dens:.4f} at q={min_at[1]}"
          f" ({min_at[2]}) of S={list(min_at[0])}")
    print(f"    min over instances of (min over pair-sum q of LM): "
          f"{min(min_LM_pairsum_per_inst)}")
    print(f"[4] sub-Vmax sanity: dead rulers {dead_sub}/{tot_sub} sampled "
          f"(divisibility kills exist below Vmax — the law needs q > Vmax)")
    print("\nIf [1] and [2] are 100%: liveness is a property of the MODULUS RANGE")
    print("(q > Vmax), not of the pair-sum structure — pair-sums matter for")
    print("ATTAINMENT (THM-668 Part 1), while the certificate supply is universal.")
    print("The a-priori target then collapses to: ONE modulus class, e.g. prove")
    print("LIVE(S, q) for q = a canonical choice (2Vmax? Vmax+v_min? a prime in")
    print("range?) for every covering 13-set.")
    print("DONE.")


if __name__ == '__main__':
    main()


def exhaustive_small():
    """Honesty check at small Vmax: all primitive covering 13-subsets of [1,18]
    (the 966 family): is every supra-Vmax q in (Vmax, 2Vmax] live? Bands are
    coarse at q <= 36, so failures would concentrate here."""
    from itertools import combinations
    from math import gcd as _g
    from functools import reduce
    tot_sets = 0
    fail_sets = 0
    tot_rulers = 0
    dead = []
    for S in combinations(range(1, 19), 13):
        if not is_covering(S):
            continue
        if reduce(_g, S) != 1:
            continue
        tot_sets += 1
        V = S[-1]
        for q in range(V + 1, 2 * V + 1):
            tot_rulers += 1
            if LM(list(S), q) == 0:
                dead.append((S, q))
    fail_sets = len({s for s, q in dead})
    print(f"\n[5] EXHAUSTIVE small-Vmax check: primitive covering 13-subsets of [1,18]: {tot_sets}")
    print(f"    supra-Vmax moduli tested: {tot_rulers}; DEAD: {len(dead)} "
          f"(on {fail_sets} sets)")
    for s, q in dead[:8]:
        print(f"      dead: q={q} S={list(s)}")
    if not dead:
        print("    UNIVERSAL SUPRA-VMAX LIVENESS holds exhaustively at small Vmax.")


if __name__ == '__main__':
    exhaustive_small()
