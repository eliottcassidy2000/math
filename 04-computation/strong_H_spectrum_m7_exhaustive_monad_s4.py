"""EXHAUSTIVE strong-tournament H-spectrum at m=7 (monad-compute-2026-06-03-S4).

HYP-2180 (opus-S599s) established the strong-tournament H-spectra EXHAUSTIVELY only for
m<=6, and asserted from a *near-transitive scan* that at m=7:
  - minH_strong(7) = 23  (formula m^2-5m+9 = 49-35+9 = 23),
  - 21 is still NOT a strong value,
  - the durable gaps 35=7*5 and 49=7^2 "fill in" via new strong primes,
  - 7 and 21 remain non-strong-values (hence permanently forbidden H).

This script makes m=7 EXHAUSTIVE. The space is 2^C(7,2) = 2^21 labeled tournaments.

KEY OPTIMIZATION (reversal halving): the converse tournament (reverse ALL arcs) has the
SAME Hamiltonian-path count H (reverse each path) and the SAME strong-connectivity. Under
the arc-bit encoding, the converse = complement of all E=C(m,2) arc bits. With E odd
(E=21 for m=7) there are NO self-converse fixed points, so iterating bits in [0, 2^(E-1))
visits exactly one representative of every reversal pair => HALF the work, SAME H-value set.

We also test is_strong (cheap) before Hcount (expensive) so non-strong tournaments are
skipped. Sanity: m<=6 spectra are recomputed with the same code and checked against the
canon values from HYP-2180.

NB on Hcount semantics (definitions.md): H(T) = number of Hamiltonian *paths* (directed),
i.e. Redei's theorem count (always odd). Verified by Hcount(transitive)=1.
"""
import sys, time
from itertools import combinations

def make_arcs(m):
    return list(combinations(range(m), 2))

def build(bits, arcs, m):
    adj = [0]*m
    for k, (i, j) in enumerate(arcs):
        if bits >> k & 1:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return adj

def Hcount(m, adj):
    size = 1 << m
    dp = [[0]*m for _ in range(size)]
    for v in range(m):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(m):
            c = row[v]
            if not c:
                continue
            av = adj[v]
            rem = ~mask
            w = 0
            while av:
                if av & 1 and (rem >> w & 1):
                    dp[mask | (1 << w)][w] += c
                av >>= 1
                w += 1
    full = size - 1
    return sum(dp[full][v] for v in range(m))

def is_strong(m, adj):
    FULL = (1 << m) - 1
    # forward reachability from 0
    seen = 1; frontier = 1
    while frontier:
        nf = 0; mm = frontier
        while mm:
            v = (mm & -mm).bit_length() - 1; mm &= mm - 1
            nf |= adj[v]
        nf &= ~seen
        if not nf:
            break
        seen |= nf; frontier = nf
    if seen != FULL:
        return False
    # reverse adjacency, reachability to 0
    radj = [0]*m
    for v in range(m):
        av = adj[v]; w = 0
        while av:
            if av & 1:
                radj[w] |= 1 << v
            av >>= 1; w += 1
    seen = 1; frontier = 1
    while frontier:
        nf = 0; mm = frontier
        while mm:
            v = (mm & -mm).bit_length() - 1; mm &= mm - 1
            nf |= radj[v]
        nf &= ~seen
        if not nf:
            break
        seen |= nf; frontier = nf
    return seen == FULL

def strong_spectrum(m, progress_every=0):
    """Set of H-values over all strongly-connected tournaments on m vertices,
    via reversal halving (iterate one rep per converse-pair)."""
    arcs = make_arcs(m)
    E = len(arcs)                 # C(m,2)
    half = 1 << (E - 1)           # bits with top arc-bit = 0; one per reversal pair
    sv = {}
    t0 = time.time()
    strong = 0
    for bits in range(half):
        adj = build(bits, arcs, m)
        if not is_strong(m, adj):
            continue
        strong += 1
        h = Hcount(m, adj)
        sv[h] = sv.get(h, 0) + 1
        if progress_every and (bits & (progress_every - 1)) == 0 and bits:
            el = time.time() - t0
            frac = bits / half
            print(f"    [m={m}] {bits}/{half} ({frac:6.2%})  strong={strong}  "
                  f"|spec|={len(sv)}  elapsed={el:6.1f}s  eta={el/frac-el:6.1f}s", flush=True)
    return sv, strong, half, time.time() - t0

def main():
    print("=" * 78)
    print("EXHAUSTIVE strong-tournament H-spectrum (monad-compute-S4)")
    print("Reversal-halving: process one rep per converse pair (E odd => no fixed pts).")
    print("=" * 78)

    CANON = {3: {3}, 4: {5}, 5: {9, 11, 13, 15},
             6: {15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45}}

    all_spectra = {}
    # ---- sanity: m=3..6 (fast) ----
    for m in range(3, 7):
        sv, strong, half, dt = strong_spectrum(m)
        spec = set(sv)
        all_spectra[m] = spec
        ok = (spec == CANON[m])
        print(f"  m={m}: strong H-values = {sorted(spec)}")
        print(f"        #strong(pairs)={strong} of {half} reps | dt={dt:.2f}s | "
              f"matches canon HYP-2180: {ok}")
        if not ok:
            print(f"        !! MISMATCH vs canon {sorted(CANON[m])}", flush=True)

    # ---- the target: m=7 exhaustive ----
    print("\n" + "-" * 78)
    print("m=7 EXHAUSTIVE (2^21 tournaments, 2^20 reps after halving) ...")
    sys.stdout.flush()
    sv7, strong7, half7, dt7 = strong_spectrum(7, progress_every=1 << 16)
    spec7 = set(sv7)
    all_spectra[7] = spec7
    print(f"\n  m=7: #strong(reps)={strong7} of {half7} | dt={dt7:.1f}s")
    print(f"  m=7 strong H-spectrum ({len(spec7)} values):")
    print(f"    {sorted(spec7)}")
    minh = min(spec7)
    print(f"\n  minH_strong(7) = {minh}  (formula m^2-5m+9 = {7*7-5*7+9} => "
          f"{'MATCH' if minh == 7*7-5*7+9 else 'MISMATCH'})")
    for q in (7, 21, 35, 49, 63):
        print(f"    is {q:3d} a strong value at m=7? {q in spec7}")

    # full minH_strong sequence
    seq = [min(all_spectra[m]) for m in range(3, 8)]
    print(f"\n  minH_strong(m) for m=3..7 = {seq}  (formula m^2-5m+9: "
          f"{[m*m-5*m+9 for m in range(3,8)]})")

    # ---- multiplicative semigroup generated by ALL strong values (m<=7) ----
    primes = set().union(*all_spectra.values())
    B = 300
    ach = {1}; changed = True
    while changed:
        changed = False
        for a in list(ach):
            for g in primes:
                ag = a * g
                if ag <= B and ag not in ach:
                    ach.add(ag); changed = True
    gaps = [h for h in range(1, B + 1, 2) if h not in ach]
    print(f"\n  Semigroup generated by strong values (m<=7), achievable odds<= {B}:")
    print(f"    forbidden (gap) odd H-values <= {B}: {gaps}")
    print(f"    => H=7 forbidden: {7 not in ach};  H=21 forbidden: {21 not in ach}")
    print(f"    => H=35 now achievable: {35 in ach};  H=49 now achievable: {49 in ach};  "
          f"H=63 achievable: {63 in ach}")

    print("\nDONE. Compare m=7 spectrum against the HYP-2180 'near-transitive scan' claims.")

if __name__ == "__main__":
    main()
