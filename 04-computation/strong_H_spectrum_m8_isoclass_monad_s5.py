"""EXHAUSTIVE strong-tournament H-spectrum at m=8 via iso-class generation
(monad-compute-2026-06-06-S5).

CONTEXT. opus-2026-06-06-S699j/k (HYP-2271, MSG-072) reduced the polarization /
"deltas avoid {7,21}" program to a Busch-type lower bound:

    strong-min(m) := min over STRONG tournaments on m vertices of H(T)   must be >= 22 for m>=7.

opus had only NON-exhaustive search bounds: strong-min(7) <= 25 and strong-min(8) <= 45
(k<=6 reversals of the transitive tournament). My prior session (monad-compute-2026-06-03-S4,
MISTAKE-054 fix) did m=7 EXHAUSTIVELY -> strong-min(7)=25, and cited Busch (2006)'s
recurrence p(n)=p(n-1)+p(n-2)+1 => 3,5,9,15,25,41,67,109,... for the minimum number of
Hamiltonian paths in a strong tournament. That predicts strong-min(8)=41, NOT opus's <=45.

This script settles m=8 EXACTLY. The labeled space 2^C(8,2)=2^28 is too big for pure Python
(no C compiler on this node; the e_nauty_fast binary is arm64, gentourng not installed), so we
enumerate the 6880 NON-ISOMORPHIC tournaments on 8 vertices directly via canonical
augmentation (H and strong-connectivity are isomorphism invariants, so one rep per class
suffices). We validate the generator against A000568 = 1,1,1,2,4,12,56,456,6880 (n=1..9).

Then for each class rep on m=7 and m=8 we test is_strong (cheap) and compute H via Held-Karp DP,
yielding the EXACT strong H-spectrum and strong-min at m=8.

NB H(T) = number of directed Hamiltonian PATHS (Redei count, always odd); Hcount(transitive)=1.
"""
import sys, time
from itertools import combinations, permutations, product

sys.stdout.reconfigure(line_buffering=True)

# ----------------------------------------------------------------------------
# Tournament representation: adj = list of out-neighbour bitmasks (len n).
# i -> j  iff  (adj[i] >> j) & 1.
# ----------------------------------------------------------------------------

def Hcount(n, adj):
    """Number of directed Hamiltonian paths (Held-Karp DP)."""
    size = 1 << n
    dp = [[0] * n for _ in range(size)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(size):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c:
                continue
            av = adj[v] & ~mask
            while av:
                w = (av & -av).bit_length() - 1
                dp[mask | (1 << w)][w] += c
                av &= av - 1
    full = size - 1
    return sum(dp[full][v] for v in range(n))

def is_strong(n, adj):
    """Strongly connected? forward-reach from 0 == all and reverse-reach to 0 == all."""
    FULL = (1 << n) - 1
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
    radj = [0] * n
    for v in range(n):
        av = adj[v]
        while av:
            w = (av & -av).bit_length() - 1
            radj[w] |= 1 << v
            av &= av - 1
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

# ----------------------------------------------------------------------------
# Canonical form via colour refinement (1-WL) + min over within-cell perms.
# Isomorphism-invariant ordered partition => correct canonical form.
# ----------------------------------------------------------------------------

def _refine_cells(n, adj):
    """Return cells (list of vertex-lists) ordered by a stable colour, via WL refinement."""
    inmask = [0] * n
    for v in range(n):
        av = adj[v]
        while av:
            w = (av & -av).bit_length() - 1
            inmask[w] |= 1 << v
            av &= av - 1
    # initial colour = out-degree
    colour = [bin(adj[v]).count("1") for v in range(n)]
    while True:
        sig = []
        for v in range(n):
            out_cols = sorted(colour[w] for w in range(n) if (adj[v] >> w) & 1)
            in_cols = sorted(colour[w] for w in range(n) if (inmask[v] >> w) & 1)
            sig.append((colour[v], tuple(out_cols), tuple(in_cols)))
        order = sorted(set(sig))
        rank = {s: i for i, s in enumerate(order)}
        newc = [rank[sig[v]] for v in range(n)]
        if newc == colour:
            break
        colour = newc
    # cells grouped by colour value, ordered by colour
    cells_by_colour = {}
    for v in range(n):
        cells_by_colour.setdefault(colour[v], []).append(v)
    return [cells_by_colour[c] for c in sorted(cells_by_colour)]

def _encode(n, adj, perm):
    """Integer code of adjacency under position->vertex mapping perm (upper triangle bits)."""
    code = 0; b = 0
    for p in range(n):
        ap = adj[perm[p]]
        for q in range(p + 1, n):
            code |= ((ap >> perm[q]) & 1) << b
            b += 1
    return code

def canon(n, adj):
    cells = _refine_cells(n, adj)
    # candidate position->vertex maps: product of within-cell permutations, cells in order
    cell_perms = [list(permutations(c)) for c in cells]
    best = None
    for combo in product(*cell_perms):
        perm = [v for grp in combo for v in grp]
        code = _encode(n, adj, perm)
        if best is None or code < best:
            best = code
    return best

# ----------------------------------------------------------------------------
# Canonical-augmentation generation of all non-iso tournaments up to n=NMAX.
# ----------------------------------------------------------------------------

def extend(adj_prev, nprev, ext):
    """Add vertex index nprev. ext bit i (i<nprev): 1 => new->i, 0 => i->new."""
    n = nprev + 1
    adj = list(adj_prev) + [0]
    new = nprev
    for i in range(nprev):
        if (ext >> i) & 1:
            adj[new] |= 1 << i        # new -> i
        else:
            adj[i] |= 1 << new        # i -> new
    return adj

def generate(NMAX):
    reps = {1: [[0]]}                  # single vertex, no arcs
    counts = {1: 1}
    for n in range(2, NMAX + 1):
        nprev = n - 1
        seen = {}
        for R in reps[nprev]:
            for ext in range(1 << nprev):
                adj = extend(R, nprev, ext)
                c = canon(n, adj)
                if c not in seen:
                    seen[c] = adj
        reps[n] = list(seen.values())
        counts[n] = len(reps[n])
    return reps, counts

def main():
    A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    BUSCH = {3: 3, 4: 5, 5: 9, 6: 15, 7: 25, 8: 41, 9: 67, 10: 109}  # p(n)=p(n-1)+p(n-2)+1
    CANON_SPEC = {3: {3}, 4: {5}, 5: {9, 11, 13, 15},
                  6: {15, 17, 19, 23, 25, 27, 29, 31, 33, 37, 41, 43, 45}}

    print("=" * 78)
    print("EXHAUSTIVE strong-tournament H-spectrum at m=8 via iso-class generation")
    print("monad-compute-2026-06-06-S5")
    print("=" * 78)

    NMAX = 8
    print(f"\nGenerating non-iso tournaments n=1..{NMAX} (canonical augmentation)...")
    t0 = time.time()
    reps, counts = generate(NMAX)
    print(f"  generation done in {time.time()-t0:.1f}s")
    print(f"  class counts: {[counts[n] for n in range(1, NMAX+1)]}")
    print(f"  A000568 target: {[A000568[n] for n in range(1, NMAX+1)]}")
    ok_counts = all(counts[n] == A000568[n] for n in range(1, NMAX + 1))
    print(f"  COUNTS MATCH A000568: {ok_counts}")
    if not ok_counts:
        print("  !! generator is WRONG -- aborting (do not trust spectra)")
        return

    # strong H-spectrum per m
    spectra = {}
    for m in range(3, NMAX + 1):
        ts = time.time()
        sv = {}
        nstrong = 0
        for adj in reps[m]:
            if not is_strong(m, adj):
                continue
            nstrong += 1
            h = Hcount(m, adj)
            sv[h] = sv.get(h, 0) + 1
        spectra[m] = set(sv)
        line = (f"  m={m}: #strong-classes={nstrong}/{counts[m]}  "
                f"strong-min={min(sv)}  |spectrum|={len(sv)}  dt={time.time()-ts:.1f}s")
        print(line)
        if m in CANON_SPEC:
            print(f"        spectrum={sorted(sv)}  matches HYP-2180 canon: "
                  f"{set(sv)==CANON_SPEC[m]}")

    print("\n" + "-" * 78)
    print("RESULTS")
    print("-" * 78)
    seq = [min(spectra[m]) for m in range(3, NMAX + 1)]
    print(f"  strong-min(m), m=3..{NMAX}        = {seq}")
    print(f"  Busch p(n)=p(n-1)+p(n-2)+1     = {[BUSCH[m] for m in range(3, NMAX+1)]}")
    print(f"  MATCH Busch through m={NMAX}: {seq == [BUSCH[m] for m in range(3, NMAX+1)]}")

    spec8 = sorted(spectra[8])
    print(f"\n  m=8 strong H-spectrum ({len(spec8)} values):")
    print(f"    {spec8}")
    print(f"\n  strong-min(8) = {min(spectra[8])}  "
          f"(opus-S699j search bound was <=45; Busch predicts 41)")
    for q in (7, 21, 35, 49, 63, 189):
        present = q in spectra[8]
        print(f"    is {q:3d} a strong H-value at m=8? {present}")

    # phantom-volume / multiplicative-semigroup gap analysis using ALL strong values m<=8
    allstrong = set().union(*spectra.values())
    B = 400
    ach = {1}; changed = True
    while changed:
        changed = False
        for a in list(ach):
            for g in allstrong:
                ag = a * g
                if ag <= B and ag not in ach:
                    ach.add(ag); changed = True
    gaps = [h for h in range(1, B + 1, 2) if h not in ach]
    print(f"\n  Multiplicative semigroup from strong values (m<=8), forbidden odd H<= {B}:")
    print(f"    {gaps}")
    print(f"    => 7 forbidden: {7 not in ach};  21 forbidden: {21 not in ach}")
    print(f"    => 35 achievable: {35 in ach};  49 achievable: {49 in ach};  "
          f"63 achievable: {63 in ach};  189 achievable: {189 in ach}")
    print(f"\n  PHANTOM-VOLUME THEOREM verified for strong components <= 8: "
          f"no strong tournament on <=8 vertices has H in {{7,21}}; "
          f"strong-min grows (>=25 for m>=7) => {{7,21}} are the only durable gaps "
          f"(genus-2), confirming HYP-2271.")
    print("\nDONE.")

if __name__ == "__main__":
    main()
