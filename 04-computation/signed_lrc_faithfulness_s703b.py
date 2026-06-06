"""SIGNED LRC -- general-config faithfulness of the cut-spectrum.  monad-explorer-S703b.

S703 established (AP_n, n<=16): folded sign-orbit < 2^{n-2} IFF C=2n-1 composite, and EVERY
collision is FOLD-ONLY (the UNFOLDED cut-spectrum of AP_n is faithful).  This sharply corrects
S699's T4 ("collisions are config automorphisms").

This script:
  (1) Push the AP "iff C composite" confirmation to n=20 (optimized, orbit sizes only).
  (2) GENERAL CONFIGS: scan ALL speed-sets and find the ones whose UNFOLDED cut-spectrum is NOT
      faithful -- these are the genuine 'spectral automorphism' configs (the real meaning S699
      reached for).  Tabulate them and check what additive structure they carry.
  (3) Correlate AP fold-collision COUNT with the factorization of C.
"""
from itertools import product, combinations
from math import gcd


def folded_orbit_size(V, C):
    n1 = len(V)
    pairs = [(i, j) for i in range(n1) for j in range(i + 1, n1)]
    seen = set()
    for bits in product((0, 1), repeat=n1 - 1):
        col = (0,) + bits
        sig = []
        for (i, j) in pairs:
            if col[i] != col[j]:
                f = (V[i] + V[j]) % C
                sig.append(min(f, C - f))
            else:
                sig.append(abs(V[i] - V[j]))   # diff < C/2, no fold
        seen.add(tuple(sorted(sig)))
    return len(seen)


def unfolded_orbit_size(V):
    n1 = len(V)
    pairs = [(i, j) for i in range(n1) for j in range(i + 1, n1)]
    seen = set()
    for bits in product((0, 1), repeat=n1 - 1):
        col = (0,) + bits
        sig = tuple(sorted((V[i] + V[j]) if col[i] != col[j] else abs(V[i] - V[j]) for (i, j) in pairs))
        seen.add(sig)
    return len(seen)


def factor(m):
    f = {}
    d = 2
    while d * d <= m:
        while m % d == 0:
            f[d] = f.get(d, 0) + 1
            m //= d
        d += 1
    if m > 1:
        f[m] = f.get(m, 0) + 1
    return f


def main():
    print("=" * 88)
    print("PART 1 (extended): AP_n folded sign-orbit vs 2^{n-2}, and factorization of C=2n-1")
    print("=" * 88)
    print(f"{'n':>3} {'C':>4} {'factor(C)':>12} {'2^(n-2)':>9} {'folded':>8} {'coll':>6} {'prime?':>7}")
    for n in range(3, 21):
        C = 2 * n - 1
        V = list(range(1, n))
        fo = folded_orbit_size(V, C)
        tot = 2 ** (n - 2)
        fac = factor(C)
        isprime = len(fac) == 1 and list(fac.values())[0] == 1
        facs = "*".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(fac.items()))
        print(f"{n:>3} {C:>4} {facs:>12} {tot:>9} {fo:>8} {tot-fo:>6} {str(isprime):>7}")

    print("\n" + "=" * 88)
    print("PART 2: GENERAL configs whose UNFOLDED cut-spectrum is NOT faithful (genuine collisions)")
    print("=" * 88)
    print("Searching all size-(n-1) speed sets with speeds in [1,B]; report unfolded-non-faithful ones.")
    for n in range(3, 7):
        B = 12 if n <= 5 else 11
        k = n - 1
        tot = 2 ** (n - 2)
        bad = []
        for V in combinations(range(1, B + 1), k):
            V = list(V)
            if unfolded_orbit_size(V) < tot:
                g = 0
                for x in V:
                    g = gcd(g, x)
                bad.append((tuple(x // g for x in V), V))
        # dedup by primitive
        prim_seen = {}
        for prim, V in bad:
            prim_seen.setdefault(prim, V)
        print(f"\n n={n} (k={k} runners, speeds in [1,{B}], 2^(n-2)={tot}): "
              f"{len(prim_seen)} primitive unfolded-non-faithful configs")
        for prim, V in sorted(prim_seen.items())[:25]:
            uo = unfolded_orbit_size(list(prim))
            # additive structure: number of coincident sums a+b=c+d (Sidon defect)
            sums = {}
            for a, b in combinations(prim, 2):
                sums[a + b] = sums.get(a + b, 0) + 1
            energy = sum(c * (c - 1) for c in sums.values())  # ordered colliding sum-pairs
            print(f"   {str(prim):<20} unfolded-orbit={uo}/{tot}  sum-collisions(2E')={energy}")

    print("\n" + "=" * 88)
    print("PART 3: is unfolded-faithfulness EXACTLY 'Sidon (B2) set'?  test the correspondence")
    print("=" * 88)
    for n in range(4, 7):
        B = 12
        k = n - 1
        tot = 2 ** (n - 2)
        nf_sidon = nf_nonsidon = f_sidon = f_nonsidon = 0
        for V in combinations(range(1, B + 1), k):
            V = list(V)
            sums = {}
            sidon = True
            for a, b in combinations(V, 2):
                if a + b in sums:
                    sidon = False
                sums[a + b] = 1
            faithful = unfolded_orbit_size(V) == tot
            if faithful and sidon:
                f_sidon += 1
            elif faithful and not sidon:
                f_nonsidon += 1
            elif not faithful and sidon:
                nf_sidon += 1
            else:
                nf_nonsidon += 1
        print(f" n={n}: faithful&Sidon={f_sidon}  faithful&non-Sidon={f_nonsidon}  "
              f"NONfaithful&Sidon={nf_sidon}  NONfaithful&non-Sidon={nf_nonsidon}")
    print("\n(If 'NONfaithful&Sidon'=0 always: Sidon => faithful. If 'faithful&non-Sidon'>0:")
    print(" non-Sidon configs can still be faithful, so faithfulness is STRICTLY WEAKER than Sidon.)")


if __name__ == "__main__":
    main()
