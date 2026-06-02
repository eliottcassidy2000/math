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
# ===========================================================================
A000568 = {1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880,9:191536,10:9733056,
           11:903753248,12:154108311168,13:48542114686912,14:28401423719122304}
A002785 = {2:1,3:2,4:2,5:8,6:12,7:88,8:176,9:2752,10:8784,11:279968,12:1492288,
           13:95458560,14:872687552}
A001174 = {1:1,2:2,3:7,4:42,5:582,6:21480,7:2142288,8:575016219,
           9:388311934228,10:617839462647568}      # oriented graphs (unlabeled)
A000273 = {1:1,2:3,3:16,4:218,5:9608,6:1540944,7:882033440,
           8:1793359192848,9:13027956824399552}     # digraphs (unlabeled)

def main():
    out = []
    def pr(s=""):
        print(s); out.append(s)

    pr("=" * 78)
    pr("SELF-CONVERSE FAMILIES VIA ANTI-AUTOMORPHISM BURNSIDE  —  monad-S565")
    pr("=" * 78)

    # ---- [1] ENGINE A reproduces OEIS total + self-converse anchors ----
    pr("\n[1] ENGINE A vs OEIS anchors (totals A000568/A001174/A000273; SC tournaments A002785)")
    anchors_ok = True
    def check(name, C, Cfix, tot_ref, sc_ref, nmax):
        nonlocal anchors_ok
        pr(f"    --- {name} (C={C}, Cfix={Cfix}) ---")
        for n in range(1, nmax + 1):
            iso, sc = burnside_counts(n, C, Cfix)
            tref = tot_ref.get(n)
            sref = sc_ref.get(n) if sc_ref else None
            tmark = "" if (tref is None or iso == tref) else f"  <-- TOTAL MISMATCH ({tref})"
            smark = "" if (sref is None or sc == sref) else f"  <-- SC MISMATCH ({sref})"
            if tmark or smark:
                anchors_ok = False
            if n <= 10:
                pr(f"      n={n:2d}: iso={iso}{tmark}   self-conv={sc}{smark}")
    check("tournaments", 2, 0, A000568, A002785, 12)
    check("oriented graphs", 3, 1, A001174, None, 10)
    check("digraphs", 4, 2, A000273, None, 9)
    pr(f"    => OEIS anchor cross-check: {'PASSED' if anchors_ok else 'FAILED'}")

    # ---- [2] independent brute-force ground truth ----
    pr("\n[2] INDEPENDENT brute-force (ENGINE B) vs ENGINE A")
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

    # ---- [3] the two NEW self-converse sequences ----
    pr("\n[3] NEW self-converse sequences (ENGINE A)")
    NMAX = 30
    sc_oriented = {}
    sc_digraph = {}
    iso_oriented = {}
    iso_digraph = {}
    for n in range(1, NMAX + 1):
        io, so = burnside_counts(n, 3, 1)
        idg, sdg = burnside_counts(n, 4, 2)
        iso_oriented[n], sc_oriented[n] = io, so
        iso_digraph[n], sc_digraph[n] = idg, sdg
    pr("\n  self-converse ORIENTED graphs, n=1..18:")
    pr("    " + ", ".join(str(sc_oriented[n]) for n in range(1, 19)))
    pr("\n  self-converse DIGRAPHS, n=1..16:")
    pr("    " + ", ".join(str(sc_digraph[n]) for n in range(1, 17)))

    # ---- [4] sanity: self-converse <= total, parity, complement-pair integrality ----
    pr("\n[4] structural sanity checks")
    sanity_ok = True
    for n in range(1, NMAX + 1):
        # self-converse count must be <= total iso count
        if sc_oriented[n] > iso_oriented[n] or sc_digraph[n] > iso_digraph[n]:
            pr(f"    n={n}: SC exceeds total!"); sanity_ok = False
        # (iso + sc) even  ->  #converse-pairs = (iso - sc)/2 integral
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
    pr(f"  OEIS anchors (A000568,A001174,A000273,A002785): "
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
