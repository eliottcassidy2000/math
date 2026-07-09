"""
lrc14_nogoodperiod_dichotomy_opus_S164.py   (opus-2026-07-08-S164, part 3)

THE DICHOTOMY (LRC14 capstone safety check): does EVERY hard-region k-cluster have a good period,
EXCEPT the tight (full-residue AP) cases?  If the only no-good-period clusters are dilated complete
APs {0,d,..,(k-1)d} with Vmax | k*d-structure (the M=1/14 tight cases, cited via LRC<=13), then the
good-period lemma is SAFE: [good period j*=O(k)]  OR  [tight AP, cited].  A NON-AP with no good
period would be a genuine gap.

For each k and every hard-region (E, Vmax), classify: has good period?  If not, is E a dilated AP
(longest-AP = k)?  Report any NON-AP no-good-period cluster (would be a red flag).
"""
import sys, random
from math import gcd


def has_good_period(E, Vmax):
    k = len(E)
    for j in range(1, Vmax):
        ph = sorted((e * j) % Vmax for e in E)
        mg = Vmax - ph[-1] + ph[0]
        for i in range(k - 1):
            mg = max(mg, ph[i + 1] - ph[i])
        if mg * 7 > Vmax:
            return True
    return False


def longest_ap(E):
    S = set(E); E = sorted(E); best = 1
    for i, a in enumerate(E):
        for b in E[i + 1:]:
            d = b - a
            if a - d in S:
                continue
            L = 2; x = b + d
            while x in S:
                L += 1; x += d
            best = max(best, L)
    return best


def is_full_ap(E):
    """E is a complete AP {a, a+d, ..., a+(k-1)d}?"""
    E = sorted(E); k = len(E)
    if k < 2:
        return True
    d = E[1] - E[0]
    return all(E[i + 1] - E[i] == d for i in range(k - 1))


def main():
    print("=" * 96)
    print("NO-GOOD-PERIOD DICHOTOMY: is every hard-region no-good-period cluster a (tight) full AP?")
    print("=" * 96)
    r = random.Random(2024)
    for k in range(8, 14):
        n_no = 0; n_tot = 0; nonap_bad = []
        no_examples = []
        # systematic small clusters + random
        cands = []
        # full APs (the expected tight cases)
        for d in range(1, 8):
            cands.append([i * d for i in range(k)])
        # near-APs (AP with defects) + random
        for d in range(1, 6):
            base = [i * d for i in range(k)]
            for jd in range(1, k):
                for dl in (-1, 1):
                    E = sorted(set(x + (dl if idx == jd else 0) for idx, x in enumerate(base)))
                    if len(E) == k and min(E) == 0:
                        cands.append(E)
        for _ in range(20000):
            spread = r.randint(k, 45)
            E = sorted(set([0] + r.sample(range(1, spread), k - 2) + [spread]))
            if len(E) == k:
                cands.append(E)
        seen = set()
        for E in cands:
            E = sorted(set(E))
            if len(E) != k or E[0] != 0:
                continue
            spread = max(E)
            for Vmax in range(spread + 1, (7 * spread) // 6 + 1):
                key = (tuple(E), Vmax)
                if key in seen:
                    continue
                seen.add(key)
                n_tot += 1
                if not has_good_period(E, Vmax):
                    n_no += 1
                    if not is_full_ap(E):
                        nonap_bad.append((tuple(E), Vmax, longest_ap(E)))
                    elif len(no_examples) < 3:
                        no_examples.append((tuple(E), Vmax))
        print(f"  k={k:2d}: {n_tot} (E,Vmax) configs; {n_no} with NO good period; "
              f"NON-AP no-good-period: {len(nonap_bad)}")
        if no_examples:
            print(f"        AP no-good-period examples: {no_examples[:2]}")
        if nonap_bad:
            print(f"        !!! NON-AP RED FLAGS: {nonap_bad[:3]}")
    print()
    print("  READING: if NON-AP no-good-period = 0 at every k, the dichotomy holds -- every hard")
    print("  cluster has a good period EXCEPT full APs (the tight M=1/14 cases, cited via LRC<=13).")
    print("  => the good-period lemma is [j*=O(k)] OR [tight AP]; no genuine counterexample.")


if __name__ == "__main__":
    main()
