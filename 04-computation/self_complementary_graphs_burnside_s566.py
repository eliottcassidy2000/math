"""
self_complementary_graphs_burnside_s566.py — monad-researcher-2026-06-02-S566

GOAL: Extend the validated anti-automorphism Burnside machinery (THM-283;
S565) to UNDIRECTED SIMPLE GRAPHS and their self-complementary subfamily:

    * simple undirected graphs       (total = OEIS A000088)
    * self-complementary graphs      (OEIS A000171)

"Self-complementary" here means: the graph is isomorphic to its complement
(swap edge <-> nonedge).  This is the "edge<->nonedge involution" that S565's
handoff flagged as the parallel target to orientation-reversal.

------------------------------------------------------------------------------
KEY DIFFERENCE FROM S565 DIRECTED FAMILIES
------------------------------------------------------------------------------
In directed structures (tournament/oriented/digraph), relabelling by g maps
arc i->j to arc g(i)->g(j).  If g^L transposes the two vertices of a pair
orbit (swap=1), the arc direction is FLIPPED, applying the orientation-reversal
involution iota_or.  The S565 engine therefore includes iota_or^{swap} in the
orbit monodromy.

For UNDIRECTED graphs, {i,j} and {j,i} are the same edge with the same color
(edge or nonedge).  Transposing the two vertices does NOT flip the color.  The
effective monodromy for UNDIRECTED structures ignores swap:

    plain iso  : every orbit contributes C=2 choices  (no swap flip)
    self-comp  : orbit of length L contributes 2 if L even, 0 if L odd
                 (complement involution 0<->1 applied L times returns to start
                  iff L is even)

------------------------------------------------------------------------------
VERIFICATION ANCHORS
------------------------------------------------------------------------------
OEIS A000088 (total undirected simple graphs):
  n=1..14: 1, 2, 4, 11, 34, 156, 1044, 12346, 274668, 12005168,
           1018997864, 165091172592, 50502031367952, 29054155657235488

OEIS A000171 (self-complementary simple graphs):
  n: 1  2  3  4  5   6  7   8    9    10  11   12      13
  a: 1, 0, 0, 1, 2,  0, 0, 10,  36,   0,  0, 720, 5765760
  (nonzero only when n ≡ 0 or 1 mod 4)

Brute-force independent verification: enumerate all 2^C(n,2) labeled graphs,
find iso classes, check if class equals its complement class.
"""
from math import factorial
from collections import Counter
from itertools import permutations
import time, os

# ===========================================================================
# SHARED BURNSIDE INFRASTRUCTURE  (copied from S565, unchanged)
# ===========================================================================

def permcount(cycle_type):
    counts = Counter(cycle_type)
    denom = 1
    for length, mult in counts.items():
        denom *= (length ** mult) * factorial(mult)
    return factorial(sum(cycle_type)) // denom

def gen_cycle_types(n):
    out = []
    def rec(rem, mx, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1):
            cur.append(p); rec(rem - p, p, cur); cur.pop()
    rec(n, n, [])
    return out

def representative(cycle_type):
    perm = [0] * sum(cycle_type)
    base = 0
    for L in cycle_type:
        for k in range(L):
            perm[base + k] = base + (k + 1) % L
        base += L
    return perm

def pair_orbit_data(perm):
    """Yield (L, swap_parity) for each orbit of unordered pairs under perm."""
    n = len(perm)
    seen = set()
    data = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen:
                continue
            a, b = i, j
            L = 0
            while True:
                seen.add((a, b))
                a, b = perm[a], perm[b]
                if a > b:
                    a, b = b, a
                L += 1
                if (a, b) == (i, j):
                    break
            x = i
            for _ in range(L):
                x = perm[x]
            swap = 0 if x == i else 1
            assert x in (i, j)
            data.append((L, swap))
    return data

# ===========================================================================
# ENGINE A — undirected Burnside  (new for S566)
# ===========================================================================

def burnside_undirected(n):
    """
    Return (total_iso, self_complementary) for undirected simple graphs on n
    vertices.

    * total_iso          -> OEIS A000088
    * self_complementary -> OEIS A000171

    Orbit monodromy for undirected: swap is IGNORED (no directional iota).
    Complement involution 0<->1 has Cfix=0.
    """
    if n <= 1:
        return 1, 1   # single graph (empty or one vertex), trivially self-comp
    nf = factorial(n)
    tot_iso = 0
    tot_sc = 0
    for ct in gen_cycle_types(n):
        perm = representative(ct)
        pc = permcount(ct)
        fix_iso = 1
        fix_sc = 1
        for (L, _swap) in pair_orbit_data(perm):
            fix_iso *= 2                          # undirected: always 2 choices
            fix_sc  *= (2 if L % 2 == 0 else 0)  # self-comp: orbit L must be even
        tot_iso += pc * fix_iso
        tot_sc  += pc * fix_sc
    assert tot_iso % nf == 0 and tot_sc % nf == 0, \
        f"non-integer Burnside at n={n}"
    return tot_iso // nf, tot_sc // nf

# ===========================================================================
# ENGINE B — brute-force independent verification  (small n)
# ===========================================================================

def brute_force_undirected(n):
    """
    Enumerate all 2^C(n,2) labeled undirected simple graphs on {0..n-1},
    group into iso classes under vertex relabeling, count:
      * total iso classes
      * classes isomorphic to their complement
    """
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(pairs)
    pos = {p: b for b, p in enumerate(pairs)}

    perm_maps = []
    for perm in permutations(range(n)):
        mp = []
        for (i, j) in pairs:
            pi, pj = perm[i], perm[j]
            if pi < pj:
                mp.append(pos[(pi, pj)])
            else:
                mp.append(pos[(pj, pi)])
        perm_maps.append(mp)

    total_cfgs = 1 << E      # 2^E labeled graphs
    visited = bytearray(total_cfgs)

    iso_count = 0
    sc_count = 0
    for v in range(total_cfgs):
        if visited[v]:
            continue
        # build iso class
        orbit = set()
        for mp in perm_maps:
            w = 0
            for k in range(E):
                if (v >> k) & 1:
                    w |= (1 << mp[k])
            visited[w] = 1
            orbit.add(w)
        iso_count += 1
        # complement: flip all bits in the E-bit word
        comp_v = v ^ ((1 << E) - 1)
        if comp_v in orbit:
            sc_count += 1
    return iso_count, sc_count

# ===========================================================================
# OEIS ANCHORS
# ===========================================================================
A000088_known = {
    1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346, 9:274668,
    10:12005168, 11:1018997864, 12:165091172592,
    13:50502031367952, 14:29054155657235488,
}

A000171_known = {
    # OEIS A000171 b-file (verified via fetch), offset=1
    1:1, 2:0, 3:0, 4:1, 5:2, 6:0, 7:0, 8:10, 9:36,
    10:0, 11:0, 12:720, 13:5600, 14:0, 15:0,
    16:703760, 17:11220000, 18:0, 19:0, 20:9168331776,
    # nonzero only when n ≡ 0 or 1 mod 4
}

def main():
    out = []
    def pr(s=""):
        print(s); out.append(s)

    pr("=" * 78)
    pr("SELF-COMPLEMENTARY UNDIRECTED GRAPHS (A000088, A000171) — monad-S566")
    pr("=" * 78)
    pr("Parallel to S565 (orientation-reversal directed families).")
    pr("Key difference: swap in pair-orbit does NOT apply iota for undirected.")

    NMAX = 40
    t_start = time.time()

    # ---- [0] compute everything via Engine A ----
    iso_vals = {}
    sc_vals = {}
    for n in range(1, NMAX + 1):
        iso_vals[n], sc_vals[n] = burnside_undirected(n)

    pr(f"\n[0] Engine A computed n=1..{NMAX} in {time.time()-t_start:.2f}s")

    # ---- [1] OEIS cross-check ----
    pr("\n[1] OEIS cross-check")
    ok88 = all(iso_vals[n] == A000088_known[n] for n in A000088_known)
    ok171 = all(sc_vals[n] == A000171_known[n] for n in A000171_known)
    for n in sorted(A000088_known):
        match = iso_vals[n] == A000088_known[n]
        if not match:
            pr(f"    A000088 MISMATCH n={n}: got {iso_vals[n]}, expected {A000088_known[n]}")
    for n in sorted(A000171_known):
        match = sc_vals[n] == A000171_known[n]
        if not match:
            pr(f"    A000171 MISMATCH n={n}: got {sc_vals[n]}, expected {A000171_known[n]}")
    pr(f"    A000088 (total graphs):            {'PASS' if ok88 else 'FAIL'} (checked n=1..{max(A000088_known)})")
    pr(f"    A000171 (self-complementary graphs): {'PASS' if ok171 else 'FAIL'} (checked n=1..{max(A000171_known)})")

    # ---- [2] brute-force independent verification ----
    pr("\n[2] Independent brute-force (Engine B) vs Engine A (n<=6)")
    brute_ok = True
    for n in range(1, 7):
        t0 = time.time()
        ib, sb = brute_force_undirected(n)
        ia, sa = iso_vals[n], sc_vals[n]
        ok = (ib == ia) and (sb == sa)
        brute_ok &= ok
        pr(f"    n={n}: brute(iso={ib}, sc={sb})  formula(iso={ia}, sc={sa})"
           f"  [{'MATCH' if ok else '*** MISMATCH ***'}]  ({time.time()-t0:.2f}s)")
    pr(f"    => brute-force verification: {'PASSED' if brute_ok else 'FAILED'}")

    # ---- [3] structural sanity ----
    pr("\n[3] Structural sanity checks")
    sanity_ok = True
    for n in range(1, NMAX + 1):
        iso, sc = iso_vals[n], sc_vals[n]
        if sc > iso:
            pr(f"    n={n}: sc={sc} > iso={iso}!"); sanity_ok = False
        if (iso - sc) % 2 != 0:
            pr(f"    n={n}: (iso-sc) odd — complement pairing not integral!")
            sanity_ok = False
        # A000171 is 0 unless n ≡ 0 or 1 mod 4
        if n % 4 not in (0, 1) and sc != 0:
            pr(f"    n={n}: sc={sc} != 0 but n mod 4 = {n%4}!")
            sanity_ok = False
    pr(f"    sc <= total, (iso-sc) even, sc=0 for n≢0,1 (mod 4): "
       f"{'OK' if sanity_ok else 'FAIL'}")

    # ---- [4] display sequences ----
    pr("\n[4] Sequences (OEIS A000088 and A000171)")
    pr(f"\n  A000088 (total simple graphs), n=1..20:")
    pr("    " + ", ".join(str(iso_vals[n]) for n in range(1, 21)))
    pr(f"\n  A000171 (self-complementary graphs), n=1..20:")
    pr("    " + ", ".join(str(sc_vals[n]) for n in range(1, 21)))
    pr(f"\n  Nonzero A000171 values only (n ≡ 0,1 mod 4, n=1..{NMAX}):")
    nz = [(n, sc_vals[n]) for n in range(1, NMAX+1) if sc_vals[n] > 0]
    for n, v in nz:
        pr(f"    A000171({n:3d}) = {v}")

    # ---- [5] ratio A000171 / A000088 for n ≡ 0,1 mod 4 ----
    pr("\n[5] Self-complementary fraction A000171(n)/A000088(n):")
    for n, sc in nz:
        frac = sc / iso_vals[n]
        pr(f"    n={n:3d}: {sc}/{iso_vals[n]} ≈ {frac:.6f}")

    # ---- [6] write b-files ----
    base = os.path.dirname(os.path.abspath(__file__))
    def write_bfile(name, vals, nmax):
        path = os.path.join(base, name)
        with open(path, 'w') as f:
            for n in range(1, nmax + 1):
                f.write(f"{n} {vals[n]}\n")
        return path
    nmax_bfile = NMAX
    p1 = write_bfile('b_graphs_total_s566.txt', iso_vals, nmax_bfile)
    p2 = write_bfile('b_graphs_selfcomp_s566.txt', sc_vals, nmax_bfile)
    pr(f"\n[6] Wrote b-files:\n    {p1}\n    {p2}")

    pr("\n" + "=" * 78)
    pr("SUMMARY")
    pr(f"  OEIS A000088 cross-check: {'PASS' if ok88 else 'FAIL'}")
    pr(f"  OEIS A000171 cross-check: {'PASS' if ok171 else 'FAIL'}")
    pr(f"  Brute-force verification: {'PASS' if brute_ok else 'FAIL'}")
    pr(f"  Structural sanity:        {'PASS' if sanity_ok else 'FAIL'}")
    pr(f"  A000171({NMAX}) = {sc_vals[NMAX]}")
    pr(f"  A000171 nonzero values for n=1..60: {[n for n,v in nz]}")
    pr("=" * 78)

    outpath = os.path.join(base, '..', '05-knowledge', 'results',
                           'self_complementary_graphs_burnside_s566.out')
    with open(outpath, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print(f"\nSaved transcript to {outpath}")

if __name__ == '__main__':
    main()
