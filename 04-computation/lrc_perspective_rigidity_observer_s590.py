#!/usr/bin/env python3
"""
claudebox-2026-06-03-S590 : perspectives = rigidity, NOT chirality — the observer-coupled vs
observer-blind reading of the human's perspective curiosity, and its LRC meaning.

The human's curiosity (T075 "perspective conjecture", T174/INV-083 A093934):
  perspectives(n) := Σ_T (# vertex-orbits of T under Aut(T))  =  # rooted/pointed tournaments
  coincides with structures(n+1) := A000568(n+1) ... but only for small n.

n=3: 2 structures, 4 perspectives (3 transitive + 1 cyclic).  4 = A000568(4).
n=4: 4 structures, 12 perspectives (4+4+2+2).               12 = A000568(5).
n=5: 12 structures, 48 perspectives.                         48 ≠ 56 = A000568(6).  <-- BREAKS.

The repo's prior reading (HYP-1824/1825) explained the 56-48 = 8 gap as a "chirality residue /
8 alpha stencils" (56 = 12 self-converse + 44 chiral). THIS FILE argues that is the wrong
symmetry: the perspective count is governed by AUTOMORPHISM RIGIDITY (vertex-orbits), not chirality
(T vs T^op). We show:
  (1) perspectives(n)/structures(n) -> n : almost every large tournament is RIGID (n perspectives);
      the deficit from n is the SYMMETRY content.
  (2) the gap T(n+1) - perspectives(n) is NOT the chiral count and does NOT track self-converse;
      it is a generic small-n artifact of both sequences briefly equalling 2(n-1)!.
  (3) the rigid<->symmetric ("observer-coupled" vs "observer-blind") split is the tournament face
      of LRC's single-corrector vs multi-sieve transition (the apex obstruction and its dissolution
      under pair-sum sieving, HYP-2063 / HYP-2075).
"""
from math import factorial, gcd
from fractions import Fraction
from collections import Counter
import itertools

# ---- fast cycle-type Burnside (A000568 etc.) ------------------------------- #
def odd_partitions(n):
    def gen(rem, maxp):
        if rem == 0:
            yield (); return
        p = maxp
        while p >= 1:
            if p % 2 == 1 and p <= rem:
                for rest in gen(rem - p, p):
                    yield (p,) + rest
            p -= 1
    return list(gen(n, n))

def perms_of_type(n, parts):
    c = Counter(parts); denom = 1
    for k, m in c.items():
        denom *= (k ** m) * factorial(m)
    return factorial(n) // denom

def pair_orbits(parts):
    orb = 0
    for i, L in enumerate(parts):
        orb += (L - 1) // 2
        for L2 in parts[i + 1:]:
            orb += gcd(L, L2)
    return orb

def Tcount(n):
    tot = Fraction(0)
    for parts in odd_partitions(n):
        tot += perms_of_type(n, parts) * 2 ** pair_orbits(parts)
    return int(tot / factorial(n))

def perspectives(n):
    tot = Fraction(0)
    for parts in odd_partitions(n):
        fixed = sum(1 for L in parts if L == 1)
        tot += perms_of_type(n, parts) * fixed * 2 ** pair_orbits(parts)
    return int(tot / factorial(n))

# ---- self-converse count (chirality), via converse-twisted Burnside -------- #
def self_converse(n):
    """# iso classes T with T ≅ T^op. = (1/n!) Σ_π #{T : π(T) = T^op}.
    π(T)=T^op fixable iff π has the right cycle parity; count 2^{(orbits where consistent)}.
    We compute directly for small n by canonical enumeration."""
    if n > 6:
        return None
    edges = list(itertools.combinations(range(n), 2))
    seen = {}
    def canon(beats):
        best = None
        for p in itertools.permutations(range(n)):
            key = tuple(beats[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
            if best is None or key < best:
                best = key
        return best
    sc = 0
    classes = {}
    for bits in itertools.product([0, 1], repeat=len(edges)):
        beats = [[0] * n for _ in range(n)]
        for (i, j), b in zip(edges, bits):
            beats[i][j] = b; beats[j][i] = 1 - b
        c = canon(beats)
        if c in classes:
            continue
        classes[c] = beats
        opp = [[beats[j][i] for j in range(n)] for i in range(n)]
        if canon(opp) == c:
            sc += 1
    return sc, len(classes)

def orbit_distribution(n):
    """exact vertex-orbit-size multiset over iso classes (small n)."""
    edges = list(itertools.combinations(range(n), 2))
    classes = {}
    for bits in itertools.product([0, 1], repeat=len(edges)):
        beats = [[0] * n for _ in range(n)]
        for (i, j), b in zip(edges, bits):
            beats[i][j] = b; beats[j][i] = 1 - b
        best = None
        for p in itertools.permutations(range(n)):
            key = tuple(beats[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
            if best is None or key < best:
                best = key
        classes.setdefault(best, beats)
    dist = Counter()
    perclass = []
    for beats in classes.values():
        auts = [p for p in itertools.permutations(range(n))
                if all(beats[i][j] == beats[p[i]][p[j]] for i in range(n) for j in range(n))]
        seen = set(); orbits = []
        for v in range(n):
            if v in seen: continue
            orb = set(p[v] for p in auts); orbits.append(len(orb)); seen |= orb
        perclass.append((len(orbits), len(auts)))
        dist[len(orbits)] += 1
    return dist, perclass

def main():
    print("=" * 78)
    print("S590  perspectives = rigidity (not chirality); the observer-coupled/blind LRC reading")
    print("=" * 78)

    print("\n[1] THE COINCIDENCE AND ITS BREAK (Burnside, exact)")
    print("  n :   T(n)   persp(n)   T(n+1)   2(n-1)!   gap=T(n+1)-persp   persp/T(n)")
    for n in range(2, 11):
        Tn, pn, Tn1, two = Tcount(n), perspectives(n), Tcount(n + 1), 2 * factorial(n - 1)
        print(f"  {n:2d}: {Tn:8d} {pn:9d} {Tn1:8d} {two:9d}   {Tn1 - pn:11d}      {pn / Tn:.4f}")
    print("  => persp(n)=T(n+1) ONLY n<=4 (both = 2(n-1)! briefly); gap explodes after.")
    print("     persp(n)/T(n) -> n: almost every large tournament is RIGID (n distinct perspectives).")

    print("\n[2] THE SYMMETRY DEFICIT = n*T(n) - persp(n)  (the 'observer-blind' mass)")
    print("  n : n*T(n)  persp(n)  deficit  deficit/T(n) = avg(n - #orbits) [symmetry per class]")
    for n in range(2, 11):
        Tn, pn = Tcount(n), perspectives(n)
        d = n * Tn - pn
        print(f"  {n:2d}: {n*Tn:9d} {pn:9d} {d:8d}   {d/Tn:.4f}")
    print("  => deficit/T(n) -> 0: the rigid (observer-coupled) regime dominates asymptotically;")
    print("     symmetric (observer-blind) classes are a vanishing fraction but set the small-n gap.")

    print("\n[3] orbit-size distribution per n (the human's 3+1, 4+4+2+2 decomposition)")
    for n in range(3, 6):
        dist, perclass = orbit_distribution(n)
        print(f"  n={n}: #orbits-per-class multiset {dict(sorted(Counter(o for o,_ in perclass).items()))}"
              f"   (sum = perspectives = {sum(o for o,_ in perclass)})")
        print(f"        per-class (orbits,|Aut|): {sorted(perclass)}")

    print("\n[4] CHIRALITY IS A DIFFERENT SYMMETRY: gap vs self-converse/chiral split")
    print("  n : structures  self-converse  chiral   gap=T(n+1)-persp(n)")
    for n in range(3, 7):
        sc = self_converse(n)
        if sc is None: continue
        scn, total = sc
        gap = Tcount(n + 1) - perspectives(n)
        print(f"  {n:2d}: {total:8d}   {scn:10d}   {total-scn:6d}   {gap:6d}")
    print("  => the perspective gap (0,0,8,160,...) does NOT equal the chiral count nor the SC count;")
    print("     HYP-1825's '8 = chirality/stencil residue' matched a number, not the mechanism.")
    print("     The mechanism is automorphism RIGIDITY (vertex-orbits), not T<->T^op chirality.")

if __name__ == "__main__":
    main()
