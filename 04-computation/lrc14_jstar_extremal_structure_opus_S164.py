"""
lrc14_jstar_extremal_structure_opus_S164.py   (opus-2026-07-08-S164, part 2)

The AP is NOT extremal for j* (part 1: random > AP at k=8,9,10).  Find the TRUE j*-maximizer
structure and test: (1) is max j* <= ceil(7(k-1)/6) universally?  (2) does j* correlate with
longest-AP (kps-S90 interlock: small j* / pigeonhole-good for longest-AP<=8)?  (3) what IS the
extremal cluster?  j* is DILATION-INVARIANT (j*(cE,cVmax)=j*(E,Vmax)) so this lives on the
longest-AP axis -- the same axis as the density floor (my S155-S163 theme).
"""
import sys, random
from math import gcd


def jstar(E, Vmax):
    k = len(E)
    for j in range(1, Vmax):
        ph = sorted((e * j) % Vmax for e in E)
        mg = Vmax - ph[-1] + ph[0]
        for i in range(k - 1):
            g = ph[i + 1] - ph[i]
            if g > mg:
                mg = g
        if mg * 7 > Vmax:
            return j
    return None


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


def hard_vmaxes(spread):
    return list(range(spread + 1, (7 * spread) // 6 + 1))


def worst_over_vmax(E):
    best = -1; bV = None
    for V in hard_vmaxes(max(E)):
        js = jstar(E, V)
        if js is None:
            return None, V
        if js > best:
            best = js; bV = V
    return best, bV


def main():
    print("=" * 96)
    print("j*-EXTREMAL STRUCTURE + longest-AP interlock  (LRC14 capstone: j*=O(k))")
    print("=" * 96)
    r = random.Random(77)
    for k in range(8, 14):
        apb = -(-7 * (k - 1) // 6)
        # thorough search: random + AP + near-AP + block-based
        cands = []
        for _ in range(12000):
            spread = r.randint(k, 55)
            mid = sorted(r.sample(range(1, spread), k - 2))
            cands.append([0] + mid + [spread])
        for d in range(1, 20):
            cands.append([i * d for i in range(k)])
            base = [i * d for i in range(k)]
            for jd in range(1, k):
                for dl in (-2, -1, 1, 2):
                    E = sorted(set([x + (dl if idx == jd else 0) for idx, x in enumerate(base)]))
                    if len(E) == k and min(E) == 0:
                        cands.append(E)
        best = (-1, None, None, None)   # (j*, E, V, longestAP)
        # stratify by longest-AP
        strat = {}
        for E in cands:
            E = sorted(set(E))
            if len(E) != k or E[0] != 0:
                continue
            j, V = worst_over_vmax(E)
            if j is None:
                print(f"  !! k={k}: NO good period in hard band for E={tuple(E)} V={V}")
                continue
            L = longest_ap(E)
            strat.setdefault(L, -1)
            strat[L] = max(strat[L], j)
            if j > best[0]:
                best = (j, tuple(E), V, L)
        smax = max((L for L in strat), default=0)
        strat_str = " ".join(f"L{L}:{strat[L]}" for L in sorted(strat))
        print(f"  k={k:2d}: MAX j*={best[0]:3d} (bound {apb})  argmax longest-AP={best[3]}  E={best[1]}")
        print(f"        j* by longest-AP: {strat_str}")
    print()
    print("  READING: max j* vs bound ceil(7(k-1)/6); the argmax's longest-AP; and whether")
    print("  small longest-AP => small j* (kps-S90 interlock).  j* is dilation-invariant.")


if __name__ == "__main__":
    main()
