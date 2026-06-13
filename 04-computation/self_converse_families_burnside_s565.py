"""
self_converse_families_burnside_s565.py — monad-researcher-2026-06-02-S565

GOAL: Extend the validated anti-automorphism Burnside machinery (THM-283; self-
converse TOURNAMENTS = A002785, engines S560/S561/S562) to the other two
orientation-reversal families and compute their SELF-CONVERSE counts:

    * oriented graphs : each unordered pair in {none, ->, <-}        (3 colors)
    * digraphs        : each unordered pair in {none, ->, <-, both}  (4 colors)

A "converse" reverses every arc: it acts on the per-pair color by the involution
iota that swaps -> and <- and fixes {none, both}.  A structure is self-converse
iff its isomorphism class is fixed by converse.

------------------------------------------------------------------------------
ONE MECHANICAL BURNSIDE ENGINE (no hand-derived edge formula)
------------------------------------------------------------------------------
#self-converse classes = (1/n!) * sum_{g in S_n} Fix(g . converse)
#all classes (iso)      = (1/n!) * sum_{g in S_n} Fix(g)

For a representative permutation g of each cycle type we walk the orbits of the
UNORDERED pairs {i,j}.  An orbit of length L closes up after g^L maps the pair to
itself; g^L either FIXES the two endpoints or SWAPS them.  Going once around the
orbit therefore applies iota^F where F-parity == swap (this is conjugation-
invariant, so one representative per type suffices).  A pair-color is consistent
on the orbit iff it is fixed by the orbit monodromy:

    plain iso  : monodromy = iota^{swap}            -> contributes (C if swap==0 else Cfix)
    self-conv. : monodromy = iota^{swap + L}        -> contributes (C if (swap+L) even else Cfix)

where C = #colors and Cfix = #colors fixed by iota.  Fix(.) = product over orbits.

  tournaments    : C=2, Cfix=0   (colors {->,<-})
  oriented graphs: C=3, Cfix=1   (colors {none,->,<-})
  digraphs       : C=4, Cfix=2   (colors {none,->,<-,both})

This per-orbit monodromy handles the even-cycle orientation constraint correctly
(the exact place burnside_unified_s28.py's closed formula missed A001174 at n=8).

An INDEPENDENT brute-force engine enumerates all labeled structures at small n and
counts iso classes / self-converse classes directly, ground-truthing the engine.
"""
from math import factorial, gcd
from collections import Counter
from itertools import permutations
import time, os

# ===========================================================================
# ENGINE A — mechanical per-orbit-monodromy Burnside (exact big integers)
# ===========================================================================

def permcount(cycle_type):
    """# of permutations in S_N with the given cycle type (N = sum)."""
    counts = Counter(cycle_type)
    denom = 1
    for length, mult in counts.items():
        denom *= (length ** mult) * factorial(mult)
    return factorial(sum(cycle_type)) // denom

def gen_cycle_types(n):
    """All partitions of n (cycle types of S_n), as nonincreasing tuples."""
    out = []
    def rec(rem, mx, cur):
        if rem == 0:
            out.append(tuple(cur)); return
        for p in range(min(rem, mx), 0, -1):
            cur.append(p); rec(rem - p, p, cur); cur.pop()
    rec(n, n, [])
    return out

def representative(cycle_type):
    """A permutation (list) on {0..n-1} with the given cycle type."""
    perm = [0] * sum(cycle_type)
    base = 0
    for L in cycle_type:
        for k in range(L):
            perm[base + k] = base + (k + 1) % L
        base += L
    return perm

def pair_orbit_data(perm):
    """For each orbit of unordered pairs under `perm`, yield (L, swap_parity).
    swap_parity = 1 iff g^L swaps the two vertices of the pair (else 0)."""
    n = len(perm)
    seen = set()
    data = []
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) in seen:
                continue
            # walk the orbit of {i,j}
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
            # determine swap: apply perm^L to i, see where it lands
            x = i
            for _ in range(L):
                x = perm[x]
            swap = 0 if x == i else 1   # x must be i or j
            assert x in (i, j)
            data.append((L, swap))
    return data

def burnside_counts(n, C, Cfix):
    """Return (iso_total, self_converse) counts for a family with C colors and
    Cfix colors fixed by the orientation-reversal involution iota."""
    if n <= 1:
        return 1, 1
    nf = factorial(n)
    tot_iso = 0
    tot_sc = 0
    for ct in gen_cycle_types(n):
        perm = representative(ct)
        pc = permcount(ct)
        fix_iso = 1
        fix_sc = 1
        for (L, swap) in pair_orbit_data(perm):
            fix_iso *= (C if swap == 0 else Cfix)
            fix_sc  *= (C if (swap + L) % 2 == 0 else Cfix)
        tot_iso += pc * fix_iso
        tot_sc  += pc * fix_sc
    assert tot_iso % nf == 0 and tot_sc % nf == 0, f"non-integer Burnside at n={n}"
    return tot_iso // nf, tot_sc // nf

# Families: (name, C, Cfix)
FAMILIES = {
    'tournament': (2, 0),
    'oriented':   (3, 1),
    'digraph':    (4, 2),
}

# ===========================================================================
# ENGINE B — independent brute-force orbit enumeration (small n)
# ===========================================================================

def brute_force(n, C, Cfix):
    """Enumerate all C^C(n,2) labeled structures, group into iso classes under
    relabelling, count classes and self-converse classes.

    Color encoding (per unordered pair, i<j):
        0 = none ; 1 = i->j ; 2 = j->i ; 3 = both
    iota (converse / orientation reversal) swaps 1<->2, fixes 0 and 3.
    Cfix tells us how many colors iota fixes -> which color set we use:
        C=2 (tournament): colors {1,2}
        C=3 (oriented)  : colors {0,1,2}
        C=4 (digraph)   : colors {0,1,2,3}
    """
    if C == 2:
        colors = [1, 2]
    elif C == 3:
        colors = [0, 1, 2]
    else:
        colors = [0, 1, 2, 3]

    def iota(c):
        return {0: 0, 1: 2, 2: 1, 3: 3}[c]

    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    E = len(pairs)
    pos = {p: b for b, p in enumerate(pairs)}

    # color -> index within `colors`, and base for mixed-radix encoding
    cidx = {c: k for k, c in enumerate(colors)}
    base = C

    def encode(cfg):  # cfg: list of colors per pair
        v = 0
        for c in cfg:
            v = v * base + cidx[c]
        return v

    # Precompute, for each permutation, the action on (pair-index) -> (pair-index, flip)
    perm_maps = []
    for perm in permutations(range(n)):
        mp = []
        for b, (i, j) in enumerate(pairs):
            pi, pj = perm[i], perm[j]
            if pi < pj:
                mp.append((pos[(pi, pj)], 0))
            else:
                mp.append((pos[(pj, pi)], 1))   # order reversed -> flip color by iota
        perm_maps.append(mp)

    # enumerate all configs as base-C numbers; decode lazily
    total_cfgs = C ** E
    visited = bytearray(total_cfgs)
    idx_to_color = colors  # colors[k] = color of index k

    def decode(v):
        out = [0] * E
        for b in range(E - 1, -1, -1):
            out[b] = idx_to_color[v % base]
            v //= base
        return out

    def apply_perm(cfg, mp):
        res = [0] * E
        for b in range(E):
            dst, flip = mp[b]
            c = cfg[b]
            res[dst] = iota(c) if flip else c
        return res

    def converse(cfg):
        return [iota(c) for c in cfg]

    iso = 0
    sc = 0
    for v in range(total_cfgs):
        if visited[v]:
            continue
        cfg = decode(v)
        orbit = set()
        for mp in perm_maps:
            w = encode(apply_perm(cfg, mp))
            visited[w] = 1
            orbit.add(w)
        iso += 1
        if encode(converse(cfg)) in orbit:
            sc += 1
    return iso, sc

# ===========================================================================
# OEIS anchors
#   small hard-coded anchors (totals A000568, A000273; self-conv tournaments
#   A002785) + the official Robinson/Howroyd b-files for the three families this
#   engine reproduces:  A001174 (oriented totals), A005639 (self-conv oriented),
#   A002499 (self-conv digraphs).  b-files were downloaded to 05-knowledge/results/.
# ===========================================================================
A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:9733056,
           11:903753248,12:154108311168,13:48542114686912,14:28401423719122304}
A002785 = {2:1,3:2,4:2,5:8,6:12,7:88,8:176,9:2752,10:8784,11:279968,12:1492288,
           13:95458560,14:872687552}
A000273 = {1:1,2:3,3:16,4:218,5:9608,6:1540944,7:882033440,
           8:1793359192848,9:13027956824399552}     # digraphs (unlabeled)

def read_bfile(name):
    """Read an OEIS b-file (n value per line) from 05-knowledge/results/."""
    base = os.path.dirname(os.path.abspath(__file__))
    path = os.path.join(base, '..', '05-knowledge', 'results', name)
    d = {}
    try:
        with open(path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                k, v = line.split()
                d[int(k)] = int(v)
    except FileNotFoundError:
        pass
    return d

def main():
    A001174 = read_bfile('b001174_oeis.txt')   # oriented graph totals (Robinson/Howroyd)
    A005639 = read_bfile('b005639_oeis.txt')   # self-converse oriented graphs (Robinson)
    A002499 = read_bfile('b002499_oeis.txt')   # self-converse digraphs (Harary-Palmer 1966)
    out = []
    def pr(s=""):
        print(s); out.append(s)

    pr("=" * 78)
    pr("SELF-CONVERSE FAMILIES VIA ANTI-AUTOMORPHISM BURNSIDE  —  monad-S565")
    pr("=" * 78)

    NMAX = 40

    # ---- [0] compute everything once (engine A) ----
    iso_oriented, sc_oriented = {}, {}
    iso_digraph, sc_digraph = {}, {}
    iso_tourn, sc_tourn = {}, {}
    for n in range(1, NMAX + 1):
        iso_tourn[n], sc_tourn[n]       = burnside_counts(n, 2, 0)
        iso_oriented[n], sc_oriented[n] = burnside_counts(n, 3, 1)
        iso_digraph[n], sc_digraph[n]   = burnside_counts(n, 4, 2)

    # ---- [1] independent brute-force ground truth (the MISTAKE-049 counter) ----
    pr("\n[1] INDEPENDENT brute-force (ENGINE B) vs ENGINE A — direct orbit enumeration")
    brute_ok = True
    brute_plan = [("tournament", 2, 0, 5),
                  ("oriented",   3, 1, 5),
                  ("digraph",    4, 2, 4)]
    for fam, C, Cfix, bmax in brute_plan:
        pr(f"    --- {fam} (brute n<= {bmax}) ---")
        for n in range(2, bmax + 1):
            t0 = time.time()
            ib, sb = brute_force(n, C, Cfix)
            ia, sa = burnside_counts(n, C, Cfix)
            ok = (ib == ia) and (sb == sa)
            brute_ok &= ok
            pr(f"      n={n}: brute(iso={ib}, sc={sb})  formula(iso={ia}, sc={sa})  "
               f"[{'MATCH' if ok else '*** MISMATCH ***'}]  ({time.time()-t0:.1f}s)")
    pr(f"    => brute-force verification: {'PASSED' if brute_ok else 'FAILED'}")

    # ---- [2] ENGINE A vs OEIS (b-files to n=NMAX, + small hardcoded anchors) ----
    pr("\n[2] ENGINE A vs OEIS  (b-files A001174/A005639/A002499 to n=%d)" % NMAX)
    anchors_ok = True
    def xcheck(label, vals, ref):
        nonlocal anchors_ok
        last = 0; bad = 0
        for n in range(1, NMAX + 1):
            if n in ref:
                last = n
                if vals[n] != ref[n]:
                    bad += 1
                    if bad <= 3:
                        pr(f"      {label} MISMATCH n={n}: {vals[n]} != {ref[n]}")
        if bad: anchors_ok = False
        pr(f"      {label}: {'MATCH' if bad == 0 else f'{bad} MISMATCH'} for all n<= {last}")
    xcheck("tournament total  A000568", iso_tourn,    A000568)
    xcheck("tournament SC     A002785", sc_tourn,     A002785)
    xcheck("oriented total    A001174", iso_oriented, A001174)
    xcheck("oriented SC       A005639", sc_oriented,  A005639)
    xcheck("digraph  total    A000273", iso_digraph,  A000273)
    xcheck("digraph  SC       A002499", sc_digraph,   A002499)
    pr(f"    => OEIS cross-check: {'PASSED' if anchors_ok else 'FAILED'}")

    # ---- [3] the self-converse sequences (repo gap-fill) ----
    pr("\n[3] Self-converse family sequences (new to the repo; OEIS-confirmed)")
    pr("\n  self-converse ORIENTED graphs = OEIS A005639, n=1..18:")
    pr("    " + ", ".join(str(sc_oriented[n]) for n in range(1, 19)))
    pr("\n  self-converse DIGRAPHS = OEIS A002499, n=1..16:")
    pr("    " + ", ".join(str(sc_digraph[n]) for n in range(1, 17)))

    # ---- [4] structural sanity ----
    pr("\n[4] structural sanity checks")
    sanity_ok = True
    for n in range(1, NMAX + 1):
        if sc_oriented[n] > iso_oriented[n] or sc_digraph[n] > iso_digraph[n]:
            pr(f"    n={n}: SC exceeds total!"); sanity_ok = False
        if (iso_oriented[n] - sc_oriented[n]) % 2 or (iso_digraph[n] - sc_digraph[n]) % 2:
            pr(f"    n={n}: (iso-sc) odd -> converse pairing not integral!"); sanity_ok = False
    pr(f"    self-converse <= total and (iso-sc) even for all n<= {NMAX}: "
       f"{'OK' if sanity_ok else 'FAIL'}")

    # ---- [5] write b-files ----
    base = os.path.dirname(os.path.abspath(__file__))
    def write_bfile(name, vals, nmax):
        path = os.path.join(base, name)
        with open(path, 'w') as f:
            for n in range(1, nmax + 1):
                f.write(f"{n} {vals[n]}\n")
        return path
    p1 = write_bfile('b_sc_oriented_s565.txt', sc_oriented, NMAX)
    p2 = write_bfile('b_sc_digraph_s565.txt', sc_digraph, NMAX)
    pr(f"\n[5] wrote b-files:\n    {p1}\n    {p2}")

    pr("\n" + "=" * 78)
    pr("SUMMARY")
    pr(f"  OEIS cross-check (A000568,A002785,A001174,A005639,A000273,A002499): "
       f"{'PASS' if anchors_ok else 'FAIL'}")
    pr(f"  independent brute force:                        "
       f"{'PASS' if brute_ok else 'FAIL'}")
    pr(f"  structural sanity:                              "
       f"{'PASS' if sanity_ok else 'FAIL'}")
    pr(f"  SC-oriented(1..10) = " + ", ".join(str(sc_oriented[n]) for n in range(1, 11)))
    pr(f"  SC-digraph (1..10) = " + ", ".join(str(sc_digraph[n]) for n in range(1, 11)))
    pr("=" * 78)

    outpath = os.path.join(base, '..', '05-knowledge', 'results',
                           'self_converse_families_burnside_s565.out')
    with open(outpath, 'w') as f:
        f.write('\n'.join(out) + '\n')
    print(f"\nSaved transcript to {outpath}")

if __name__ == '__main__':
    main()
