"""
merged_metagraph_exact_s561.py — monad-researcher-2026-06-02-S561

EXACT merged-metagraph vertex sequences, with an INDEPENDENT brute-force
verification of the self-converse closed form.

The merged metagraph G_n/Z_2 (CLAUDE.md: the PRIMARY object) factors out
complement symmetry T <-> T^op.  Its vertex count splits as:

    T(n)  = A000568(n)              total tournament iso classes
    SC(n) = A002785(n)             self-converse (= self-complementary) classes
    V_merged(n) = (T + SC)/2       vertices of G_n/Z_2
    NS(n)       = (T - SC)/2       non-self-converse merged vertices
                                    (CLAUDE.md: 0,1,2,22,184 for n=3..7)

WHY THIS SESSION:
  - The exact self-converse closed form (self_comp_tournaments_s90dl.py) is
    validated only against OEIS A002785 (n<=14).  Beyond that it is an
    UNVERIFIED extension.  Here we add a fully INDEPENDENT check: direct
    orbit enumeration of ALL tournaments for n=3..7, canonicalising under the
    full relabelling group, and counting self-converse classes by hand.
    This confirms the FORMULA structurally (both even- and odd-n branches),
    not merely by OEIS lookup.
  - The merged sequences V_merged and NS were never saved as standalone exact
    sequences; we extend them to n=50.

Two independent engines for T(n) and SC(n):
  ENGINE A (closed form): Burnside over cycle types.
      T(n)  sums 2^t over partitions of n into ODD parts
            (the only sigma fixing any tournament), t = orbits of pairs.
      SC(n) sums over anti-automorphism cycle types: all cycles length
            2*(odd) (i.e. == 2 mod 4) plus, for odd n, one fixed point.
  ENGINE B (brute force, n<=7): enumerate all 2^C(n,2) tournaments, partition
      into iso classes by the relabelling action, count classes (=T) and
      self-converse classes (canon(reverse(rep)) == rep) (=SC).

Both engines are big-integer exact.  ENGINE A is used for the n=50 extension;
ENGINE B is the independent witness at small n.
"""

from fractions import Fraction
from math import gcd, factorial
from collections import Counter
from itertools import permutations
import time

# ---------------------------------------------------------------------------
# ENGINE A — closed form (Burnside over cycle types)
# ---------------------------------------------------------------------------

def odd_partitions(n, max_part=None):
    """Partitions of n into odd parts (nonincreasing), as tuples."""
    if max_part is None:
        max_part = n
    if n == 0:
        yield ()
        return
    top = min(n, max_part)
    if top % 2 == 0:
        top -= 1
    p = top
    while p >= 1:
        for rest in odd_partitions(n - p, p):
            yield (p,) + rest
        p -= 2

def permcount(cycle_type):
    """Size of the conjugacy class in S_N with the given cycle type."""
    counts = Counter(cycle_type)
    denom = 1
    for length, mult in counts.items():
        denom *= (length ** mult) * factorial(mult)
    return factorial(sum(cycle_type)) // denom

def pair_orbits_auto(parts):
    """t(sigma): orbits of UNORDERED pairs under sigma with cycle type `parts`
    (all parts odd).  t = sum (c-1)/2 + sum_{i<j} gcd(c_i,c_j)."""
    s = sum((c - 1) // 2 for c in parts)
    for i in range(len(parts)):
        for j in range(i + 1, len(parts)):
            s += gcd(parts[i], parts[j])
    return s

def T_closed(n):
    """A000568(n) via exact Burnside."""
    if n <= 1:
        return 1
    total = Fraction(0)
    nf = factorial(n)
    for parts in odd_partitions(n):
        total += Fraction(permcount(parts) * (1 << pair_orbits_auto(parts)), nf)
    assert total.denominator == 1, f"T({n}) not integer!"
    return int(total)

def SC_closed(n):
    """A002785(n): self-converse tournament classes via anti-automorphism Burnside.

    An anti-automorphism g of a tournament has all cycles of length 2*(odd)
    (== 2 mod 4), plus exactly one fixed point when n is odd.  Counting
    tournaments fixed by g-then-reverse over each such cycle type gives:

        SC(n) = (1/n!) * sum_{p odd-partition of floor(n/2)}
                   permcount(2*p) * 2^e(p) * extra(n, p)

    with e(p) = sum p_i + 2*sum_{i<j} gcd(p_i,p_j)  and
    extra = 1 (n even) or n*2^len(p) (n odd).
    """
    if n <= 1:
        return 1
    m = n // 2
    total = Fraction(0)
    nf = factorial(n)
    for p in odd_partitions(m):
        doubled = tuple(2 * x for x in p)
        e = sum(p)
        for i in range(len(p)):
            for j in range(i + 1, len(p)):
                e += 2 * gcd(p[i], p[j])
        extra = (n * (1 << len(p))) if (n % 2 == 1) else 1
        total += Fraction(permcount(doubled) * (1 << e) * extra, nf)
    assert total.denominator == 1, f"SC({n}) not integer!"
    return int(total)

# ---------------------------------------------------------------------------
# ENGINE B — brute force orbit enumeration (independent), n <= 7
# ---------------------------------------------------------------------------

def brute_force_TSC(n):
    """Return (T, SC) by enumerating all tournaments on n vertices and
    grouping into isomorphism classes under the full relabelling group.

    Bit b for pair (i,j), i<j, encodes orientation: bit=1 -> i->j, 0 -> j->i.
    Relabelling by a permutation perm acts on bits (with a flip whenever the
    image pair (perm[i],perm[j]) is in decreasing order).  The converse
    tournament flips EVERY bit.  A class is self-converse iff its converse
    representative lies in the same orbit.
    """
    pairs = [(i, j) for i in range(n) for j in range(i, n) if i < j]
    E = len(pairs)
    pos = {}
    for b, (i, j) in enumerate(pairs):
        pos[(i, j)] = b
    full = (1 << E) - 1

    # Precompute each permutation as a list of (src_bit, dst_bit, flip).
    perm_maps = []
    for perm in permutations(range(n)):
        mp = []
        for b, (i, j) in enumerate(pairs):
            pi, pj = perm[i], perm[j]
            if pi < pj:
                dst, flip = pos[(pi, pj)], 0
            else:
                dst, flip = pos[(pj, pi)], 1
            mp.append((b, dst, flip))
        perm_maps.append(mp)

    def transform(mask, mp):
        r = 0
        for src, dst, flip in mp:
            bit = ((mask >> src) & 1) ^ flip
            if bit:
                r |= (1 << dst)
        return r

    visited = bytearray(1 << E)
    T_count = 0
    SC_count = 0
    for mask in range(1 << E):
        if visited[mask]:
            continue
        # new representative; build its orbit
        orbit = set()
        for mp in perm_maps:
            t = transform(mask, mp)
            if not visited[t]:
                visited[t] = 1
            orbit.add(t)
        T_count += 1
        if (mask ^ full) in orbit:        # converse isomorphic -> self-converse
            SC_count += 1
    return T_count, SC_count

# ---------------------------------------------------------------------------
# Known OEIS anchors
# ---------------------------------------------------------------------------
KNOWN_T = {
    1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880,
    9: 191536, 10: 9733056, 11: 903753248, 12: 154108311168,
    13: 48542114686912, 14: 28401423719122304,
    15: 31021002160355166848, 16: 63530415842308265100288,
    17: 244912778438520759443245824,
    18: 1783398846284777975419600287232,
    19: 24605641171260376770598003978281472,
    20: 645022068557873570931850526424042500096,
}
KNOWN_SC = {
    1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176,
    9: 2752, 10: 8784, 11: 279968, 12: 1492288,
    13: 95458560, 14: 872687552,
}
# CLAUDE.md: NS-merged node count = 0,1,2,22,184 for n=3..7
KNOWN_NS = {3: 0, 4: 1, 5: 2, 6: 22, 7: 184}


def main():
    NMAX = 50
    BRUTE_MAX = 7
    out = []
    def pr(s=""):
        print(s); out.append(s)

    pr("=" * 78)
    pr("MERGED-METAGRAPH EXACT SEQUENCES  —  monad-researcher-2026-06-02-S561")
    pr("=" * 78)

    # ---- Stage 1: independent brute-force verification (n=3..7) ----
    pr("\n[1] INDEPENDENT brute-force verification (ENGINE B vs ENGINE A)")
    pr("    direct orbit enumeration of ALL tournaments, n=3..%d" % BRUTE_MAX)
    pr(f"    {'n':>2} | {'T_brute':>8} {'T_form':>8} | {'SC_brute':>8} {'SC_form':>8} | {'status'}")
    pr("    " + "-" * 60)
    all_ok = True
    for n in range(3, BRUTE_MAX + 1):
        t0 = time.time()
        tb, scb = brute_force_TSC(n)
        tf, scf = T_closed(n), SC_closed(n)
        ok = (tb == tf) and (scb == scf)
        all_ok = all_ok and ok
        status = "MATCH" if ok else "*** MISMATCH ***"
        dt = time.time() - t0
        pr(f"    {n:>2} | {tb:>8} {tf:>8} | {scb:>8} {scf:>8} | {status} ({dt:.1f}s)")
    pr(f"    => independent verification {'PASSED' if all_ok else 'FAILED'} "
       f"(both even & odd n branches confirmed by direct enumeration)")

    # ---- Stage 2: cross-check closed form vs OEIS anchors ----
    pr("\n[2] Closed form vs OEIS anchors  (T: A000568 n<=20, SC: A002785 n<=14)")
    oeis_ok = True
    for n in range(1, 21):
        tf = T_closed(n)
        if n in KNOWN_T and tf != KNOWN_T[n]:
            pr(f"    T({n}) MISMATCH: {tf} != {KNOWN_T[n]}"); oeis_ok = False
    for n in range(1, 15):
        scf = SC_closed(n)
        if n in KNOWN_SC and scf != KNOWN_SC[n]:
            pr(f"    SC({n}) MISMATCH: {scf} != {KNOWN_SC[n]}"); oeis_ok = False
    pr(f"    => OEIS cross-check {'PASSED' if oeis_ok else 'FAILED'}")

    # ---- Stage 3: NS cross-check vs CLAUDE.md ----
    pr("\n[3] NS = (T-SC)/2 vs CLAUDE.md stated values (n=3..7)")
    ns_ok = True
    for n in range(3, 8):
        ns = (T_closed(n) - SC_closed(n)) // 2
        exp = KNOWN_NS[n]
        mark = "OK" if ns == exp else "MISMATCH"
        if ns != exp: ns_ok = False
        pr(f"    NS({n}) = {ns}  (CLAUDE.md: {exp})  {mark}")
    pr(f"    => NS cross-check {'PASSED' if ns_ok else 'FAILED'}")

    # ---- Stage 4: exact extension to n=NMAX ----
    pr(f"\n[4] EXACT merged-metagraph sequences, n=1..{NMAX}")
    pr(f"    (T must be even - SC for V_merged,NS integral; all big-int exact)")
    pr(f"    {'n':>3} | {'V_merged=(T+SC)/2':>26} | {'NS=(T-SC)/2':>26}")
    pr("    " + "-" * 72)
    Vvals, NSvals, Tvals, SCvals = {}, {}, {}, {}
    for n in range(1, NMAX + 1):
        t0 = time.time()
        tn = T_closed(n)
        scn = SC_closed(n)
        assert (tn - scn) % 2 == 0, f"T-SC odd at n={n}!"
        V = (tn + scn) // 2
        NS = (tn - scn) // 2
        Tvals[n], SCvals[n], Vvals[n], NSvals[n] = tn, scn, V, NS
        dt = time.time() - t0
        if n <= 12 or n % 5 == 0 or n == NMAX:
            pr(f"    {n:>3} | {V:>26} | {NS:>26}")

    # ---- Stage 5: write b-files for the merged sequences ----
    import os
    base = os.path.dirname(os.path.abspath(__file__))
    def write_bfile(name, vals):
        path = os.path.join(base, name)
        with open(path, 'w') as f:
            for n in sorted(vals):
                f.write(f"{n} {vals[n]}\n")
        return path
    pV = write_bfile('b_vmerged_s561.txt', Vvals)
    pN = write_bfile('b_nsmerged_s561.txt', NSvals)
    pr(f"\n[5] Wrote b-files:\n    {pV}\n    {pN}")

    pr("\n" + "=" * 78)
    pr("SUMMARY")
    pr(f"  independent brute-force (n<=7): {'PASS' if all_ok else 'FAIL'}")
    pr(f"  OEIS anchors:                   {'PASS' if oeis_ok else 'FAIL'}")
    pr(f"  CLAUDE.md NS values:            {'PASS' if ns_ok else 'FAIL'}")
    pr(f"  V_merged(3..10) = " + ", ".join(str(Vvals[n]) for n in range(3, 11)))
    pr(f"  NS(3..10)       = " + ", ".join(str(NSvals[n]) for n in range(3, 11)))
    pr("=" * 78)

    outpath = os.path.join(base, '..', '05-knowledge', 'results',
                           'merged_metagraph_exact_s561.out')
    with open(outpath, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print(f"\nSaved transcript to {outpath}")


if __name__ == '__main__':
    main()
