"""
opus-2026-05-27-S2: Sequence exploration session.

Compute and analyze sequences arising from tournament theory.
All computations bounded to run in reasonable time.
"""

import itertools
from functools import lru_cache
from collections import defaultdict, Counter

# ─────────────────────────────────────────────────────────────────────────────
# Core utilities
# ─────────────────────────────────────────────────────────────────────────────

def all_labeled_tournaments(n):
    """Yield all 2^C(n,2) labeled n-vertex tournaments as adj-dicts."""
    verts = list(range(n))
    edges = [(i,j) for i in verts for j in verts if i < j]
    for bits in range(1 << len(edges)):
        T = {v: set() for v in verts}
        for idx, (i,j) in enumerate(edges):
            if (bits >> idx) & 1:
                T[i].add(j)
            else:
                T[j].add(i)
        yield T

def is_sc(T):
    n = len(T)
    if n <= 1: return True
    def reach_fwd(src):
        vis = {src}; q = [src]
        while q:
            v = q.pop()
            for w in T[v]:
                if w not in vis: vis.add(w); q.append(w)
        return vis
    rev = {v: set() for v in T}
    for v in T:
        for w in T[v]: rev[w].add(v)
    if len(reach_fwd(0)) < n: return False
    vis = {0}; q = [0]
    while q:
        v = q.pop()
        for w in rev[v]:
            if w not in vis: vis.add(w); q.append(w)
    return len(vis) == n

def ham_paths(T):
    n = len(T)
    verts = list(T.keys())
    return sum(1 for p in itertools.permutations(verts)
               if all(p[i+1] in T[p[i]] for i in range(n-1)))

def score_seq(T):
    return tuple(sorted(len(v) for v in T.values()))

def outdeg(T, v):
    return len(T[v])

def is_king(v, T):
    r = T[v] | {v}
    for u in T[v]: r |= T[u]
    return len(r) == len(T)

def num_kings(T):
    return sum(1 for v in T if is_king(v, T))

# ─────────────────────────────────────────────────────────────────────────────
# Tiling model
# ─────────────────────────────────────────────────────────────────────────────

def get_tiles(n):
    tiles = []
    for y in range(n-2):
        for x in range(n-1, y+1, -1):
            tiles.append((x, y))
    return tiles

def is_sc_tiling(n, bits, tiles=None):
    if tiles is None: tiles = get_tiles(n)
    for k in range(1, n):
        if not any((bits >> idx) & 1 for idx, (x,y) in enumerate(tiles) if x >= k > y):
            return False
    return True

def nonsc_ie(n):
    """Non-SC tiling count via inclusion-exclusion over cuts."""
    tiles = get_tiles(n)
    m = len(tiles)
    cuts = list(range(1, n))
    total = 0
    for size in range(1, len(cuts)+1):
        for S in itertools.combinations(cuts, size):
            frozen = len({idx for idx,(x,y) in enumerate(tiles)
                          if any(x >= k > y for k in S)})
            total += ((-1)**(size+1)) * (1 << (m - frozen))
    return total

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 1: Non-SC and SC tiling counts (fast: up to n=9 via IE formula)
# ─────────────────────────────────────────────────────────────────────────────

def seq1_nonsc_tiling(max_n=10):
    print("=" * 65)
    print("SEQUENCE 1: SC / Non-SC Tiling Counts")
    print("=" * 65)
    sc_seq, nonsc_seq, total_seq = [], [], []
    for n in range(3, max_n+1):
        tiles = get_tiles(n)
        m = len(tiles)
        total = 1 << m
        nonsc = nonsc_ie(n)
        sc = total - nonsc
        sc_seq.append(sc); nonsc_seq.append(nonsc); total_seq.append(total)
        print(f"n={n}: m={m}, total=2^{m}={total}, SC={sc}, non-SC={nonsc}, "
              f"P(SC)={sc/total:.6f}")
    print()
    print(f"SC     (n=3..{max_n}): {sc_seq}")
    print(f"Non-SC (n=3..{max_n}): {nonsc_seq}")
    print()
    # Ratios
    print("SC ratios (consecutive):    ", [f"{sc_seq[i+1]/sc_seq[i]:.4f}" for i in range(len(sc_seq)-1)])
    print("Non-SC ratios (consecutive):", [f"{nonsc_seq[i+1]/nonsc_seq[i]:.4f}" for i in range(len(nonsc_seq)-1)])
    print()
    # Verify IE against brute-force for small n
    print("Verifying IE against brute-force (n=3..7):")
    for n in range(3, 8):
        tiles = get_tiles(n)
        m = len(tiles)
        bf = sum(1 for bits in range(1<<m) if not is_sc_tiling(n, bits, tiles))
        ie = nonsc_ie(n)
        print(f"  n={n}: brute-force={bf}, IE={ie}, match={bf==ie}")
    print()
    return sc_seq, nonsc_seq

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 2: Detailed IE structure — contribution by cut-subset size
# ─────────────────────────────────────────────────────────────────────────────

def seq2_ie_structure(max_n=8):
    print("=" * 65)
    print("SEQUENCE 2: IE Structure by Subset Size")
    print("=" * 65)
    for n in range(3, max_n+1):
        tiles = get_tiles(n)
        m = len(tiles)
        cuts = list(range(1, n))
        print(f"n={n} (m={m}, cuts={cuts}):")
        total_ie = 0
        for size in range(1, len(cuts)+1):
            size_contrib = 0
            for S in itertools.combinations(cuts, size):
                frozen = len({idx for idx,(x,y) in enumerate(tiles)
                              if any(x >= k > y for k in S)})
                term = ((-1)**(size+1)) * (1 << (m - frozen))
                size_contrib += term
            total_ie += size_contrib
            print(f"  size {size}: contribution = {size_contrib}")
        print(f"  total non-SC = {total_ie}")
        # Cut sizes
        cut_sizes = {k: sum(1 for (x,y) in tiles if x >= k > y) for k in cuts}
        print(f"  cut sizes f(n,k)=k(n-k)-1: {cut_sizes}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 3: Score sequence counts (A000571 verification)
# ─────────────────────────────────────────────────────────────────────────────

def seq3_score_seqs(max_n=7):
    print("=" * 65)
    print("SEQUENCE 3: Score Sequence Counts vs A000571")
    print("=" * 65)
    A000571 = {1:1, 2:1, 3:1, 4:2, 5:4, 6:9, 7:22, 8:59, 9:167, 10:490,
               11:1486, 12:4639, 13:14805, 14:48107, 15:158808, 16:531469}
    counts = []
    for n in range(1, max_n+1):
        if n <= 2:
            c = 1
        else:
            c = len({score_seq(T) for T in all_labeled_tournaments(n)})
        counts.append(c)
        exp = A000571.get(n, '?')
        print(f"n={n}: {c} score sequences  (A000571={exp}) {'✓' if c==exp else '✗'}")
    print()
    print(f"Sequence: {counts}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 4: Labeled SC tournament counts
# ─────────────────────────────────────────────────────────────────────────────

def seq4_labeled_sc(max_n=6):
    print("=" * 65)
    print("SEQUENCE 4: Labeled SC Tournament Counts")
    print("=" * 65)
    # From A054923? or direct computation
    results = []
    for n in range(1, max_n+1):
        if n == 1:
            c = 1
        elif n == 2:
            c = 0
        else:
            c = sum(1 for T in all_labeled_tournaments(n) if is_sc(T))
        results.append((n, c, 2**(n*(n-1)//2)))
        print(f"n={n}: {c} labeled SC tournaments (out of {2**(n*(n-1)//2)} total)")
    print()
    seq = [r[1] for r in results]
    print(f"Labeled SC sequence: {seq}")
    # Compare to SC tiling counts × 2^(n-1)?
    sc_tiling = [1, 5, 50, 903, 30773]  # n=3..7
    print("SC_tilings × 2^(n-1):", [sc_tiling[i] * 2**(i+3-1) for i in range(len(sc_tiling))])
    print("(This would hold if SC fraction same across base-path choices)")
    print()

    # Verify: does labeled SC count = SC fraction of tilings × total labeled tournaments?
    for n in range(3, min(max_n+1, 6)):
        total_labeled = 2**(n*(n-1)//2)
        total_tilings = 2**((n-1)*(n-2)//2)
        sc_tiling_n = [1, 5, 50, 903][n-3]
        labeled_sc_n = seq[n-1]
        frac_tiling = sc_tiling_n / total_tilings
        frac_labeled = labeled_sc_n / total_labeled
        print(f"n={n}: P_SC(tiling)={frac_tiling:.6f}, P_SC(labeled)={frac_labeled:.6f}, "
              f"equal={abs(frac_tiling-frac_labeled)<1e-9}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 5: H-value distributions for n=3..5
# ─────────────────────────────────────────────────────────────────────────────

def seq5_h_distributions(max_n=5):
    print("=" * 65)
    print("SEQUENCE 5: H-Value Distributions")
    print("=" * 65)
    all_h_seqs = {}
    for n in range(3, max_n+1):
        h_ct = Counter()
        for T in all_labeled_tournaments(n):
            h_ct[ham_paths(T)] += 1
        all_h_seqs[n] = h_ct
        vals = sorted(h_ct.keys())
        print(f"n={n}: H∈{vals}")
        print(f"  counts: {dict(h_ct)}")
        print(f"  total: {sum(h_ct.values())} = 2^{(n*(n-1)//2)}")

        # Achievable H values
        max_h = max(vals)
        missing = [h for h in range(1, max_h+2, 2) if h not in h_ct]
        print(f"  missing odd H in [1,{max_h}]: {missing}")
    print()

    # Key sequences from H distributions:
    print("Key sequences:")
    # Max H achievable at n:
    max_h_seq = [max(all_h_seqs[n].keys()) for n in range(3, max_n+1)]
    print(f"  Max H(n) for n=3..{max_n}: {max_h_seq}")
    print(f"  = {[h for h in max_h_seq]} (known: 3,5,15,...)")
    # Count of tournaments with max H:
    max_h_count = [all_h_seqs[n][max(all_h_seqs[n].keys())] for n in range(3, max_n+1)]
    print(f"  #tournaments achieving max H: {max_h_count}")
    # Count with min H=1:
    min_h_count = [all_h_seqs[n].get(1, 0) for n in range(3, max_n+1)]
    print(f"  #tournaments with H=1: {min_h_count}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 6: #Kings distributions for n=3..6
# ─────────────────────────────────────────────────────────────────────────────

def seq6_kings_distributions(max_n=6):
    print("=" * 65)
    print("SEQUENCE 6: #Kings Distributions")
    print("=" * 65)

    for n in range(3, max_n+1):
        kc = Counter()
        sc_by_k = Counter()
        for T in all_labeled_tournaments(n):
            k = num_kings(T)
            kc[k] += 1
            if is_sc(T):
                sc_by_k[k] += 1
        vals = sorted(kc.keys())
        print(f"n={n}: #kings distribution:")
        for k in vals:
            print(f"  k={k}: {kc[k]} tournaments, {sc_by_k[k]} SC")
        impossible = [k for k in range(1, n+1) if k not in kc]
        print(f"  impossible: {impossible}")
        print()

    # Sequence: for each n, #tournaments with exactly n kings (all universal kings)
    print("All-kings sequences:")
    for n in range(3, 7):
        total_n = sum(1 for T in all_labeled_tournaments(n) if num_kings(T) == n)
        print(f"  n={n}: {total_n} tournaments with ALL vertices as kings")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 7: Min SC excess = H(T)-H(T-Q)-2|rivals| for SC tournaments
# ─────────────────────────────────────────────────────────────────────────────

def seq7_sc_excess(max_n=6):
    print("=" * 65)
    print("SEQUENCE 7: Min SC Excess (HYP-1740 verification)")
    print("=" * 65)

    for n in range(3, max_n+1):
        min_ex = None
        min_case = None
        for T in all_labeled_tournaments(n):
            if not is_sc(T): continue
            degs = {v: outdeg(T, v) for v in T}
            max_d = max(degs.values())
            for Q in [v for v, d in degs.items() if d == max_d]:
                rivals = [v for v in T if Q in T[v]]
                if not rivals: continue
                # H(T)
                hT = ham_paths(T)
                # T - Q
                Tq = {v: T[v] - {Q} for v in T if v != Q}
                hTQ = ham_paths(Tq)
                excess = (hT - hTQ) - 2*len(rivals)
                if min_ex is None or excess < min_ex:
                    min_ex = excess
                    min_case = (score_seq(T), max_d, len(rivals), hT-hTQ)

        print(f"n={n}: min SC excess = {min_ex}  (case: score={min_case[0]}, deg_Q={min_case[1]}, "
              f"|rivals|={min_case[2]}, delta={min_case[3]})")
    print()
    print("Pattern: 0, 0, 2, 4 for n=3..6 => conj = 2(n-4) for n>=4")
    print("Formula: max(0, 2(n-4)) = 0,0,2,4,6,8,... for n=3,4,5,6,7,8,...")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 8: H distribution by score sequence (n=3..5)
# ─────────────────────────────────────────────────────────────────────────────

def seq8_h_by_score(max_n=5):
    print("=" * 65)
    print("SEQUENCE 8: H Distribution by Score Sequence")
    print("=" * 65)
    for n in range(3, max_n+1):
        print(f"n={n}:")
        by_ss = defaultdict(list)
        for T in all_labeled_tournaments(n):
            by_ss[score_seq(T)].append(ham_paths(T))
        for ss in sorted(by_ss.keys()):
            hs = sorted(set(by_ss[ss]))
            ct = Counter(by_ss[ss])
            unique_h = len(hs) == 1
            print(f"  {ss}: H∈{hs}  {'(UNIQUE H)' if unique_h else '(MULTIPLE H)'}  "
                  f"counts={[ct[h] for h in hs]}  total={len(by_ss[ss])}")
        print()

    # Key: for each n, which score sequences uniquely determine H?
    print("Score sequences that uniquely determine H:")
    for n in range(3, max_n+1):
        by_ss = defaultdict(set)
        for T in all_labeled_tournaments(n):
            by_ss[score_seq(T)].add(ham_paths(T))
        unique = [(ss, min(hs)) for ss, hs in by_ss.items() if len(hs)==1]
        multi = [(ss, sorted(hs)) for ss, hs in by_ss.items() if len(hs)>1]
        print(f"  n={n}: {len(unique)} unique, {len(multi)} non-unique score sequences")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 9: Cut profile analysis of non-SC tilings
# ─────────────────────────────────────────────────────────────────────────────

def seq9_cut_profiles(max_n=8):
    print("=" * 65)
    print("SEQUENCE 9: Non-SC Tiling Cut Profile (via IE)")
    print("=" * 65)
    # For each n, count non-SC tilings with exactly j "bad" cuts
    # (cuts where all tiles are downward)
    # Use IE to compute these counts efficiently.

    for n in range(3, max_n+1):
        tiles = get_tiles(n)
        m = len(tiles)
        cuts = list(range(1, n))

        # Count tilings with EXACTLY set S of bad cuts = those where
        # every k in S has all tiles downward, AND every k NOT in S has at least one upward.
        # This is Möbius inversion. Use inclusion-exclusion:
        # #{exactly bad S} = Σ_{T ⊇ S} (-1)^{|T|-|S|} * 2^{m - |⋃_{k∈T} cut_k|}

        # Simpler: for each subset S, compute #{tilings where AT LEAST S are all-bad}
        at_least = {}
        for S in itertools.chain([()], *[itertools.combinations(cuts, r) for r in range(1, len(cuts)+1)]):
            frozen = len({idx for idx,(x,y) in enumerate(tiles) if any(x>=k>y for k in S)})
            at_least[frozenset(S)] = 1 << (m - frozen)

        # Möbius: #{exactly S bad} = Σ_{T ⊇ S} (-1)^{|T\S|} * at_least[T]
        exactly = {}
        for S in itertools.chain([()], *[itertools.combinations(cuts, r) for r in range(1, len(cuts)+1)]):
            fs = frozenset(S)
            val = 0
            remaining = [k for k in cuts if k not in fs]
            for r in range(len(remaining)+1):
                for ext in itertools.combinations(remaining, r):
                    T = fs | frozenset(ext)
                    val += ((-1)**r) * at_least[T]
            exactly[fs] = val

        # Count by number of bad cuts
        by_count = Counter()
        for S, cnt in exactly.items():
            if len(S) > 0:  # non-SC tilings have at least 1 bad cut
                by_count[len(S)] += cnt

        nonsc_total = sum(v for k,v in by_count.items())
        print(f"n={n}: non-SC={nonsc_total} = {nonsc_ie(n)} (check: {nonsc_total==nonsc_ie(n)})")
        print(f"  by #bad cuts: {dict(sorted(by_count.items()))}")

        # The sequence of "exactly 1 bad cut" counts:
        exactly1 = by_count.get(1, 0)
        print(f"  exactly 1 bad cut: {exactly1}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 10: New — #transitive tilings (all cuts are bad = transitive tournament)
# ─────────────────────────────────────────────────────────────────────────────

def seq10_transitive_tilings():
    print("=" * 65)
    print("SEQUENCE 10: Transitive-like Tilings")
    print("=" * 65)
    # A tiling is "all-bad" if ALL n-1 cuts are bad.
    # From the base path 0→1→...→(n-1): if every tile is downward (bit=0),
    # we get the transitive tournament! (The base path already flows downward.)
    # But "all bad" means all tiles are downward AND all base-path arcs exist.
    # That's just 1 tiling: all bits = 0.

    # More interesting: tilings where exactly k cuts are "good" (SC-witnessing).
    for n in range(3, 9):
        tiles = get_tiles(n)
        m = len(tiles)
        cuts = list(range(1, n))
        # Count by #good cuts
        if n <= 7:
            by_good = Counter()
            for bits in range(1 << m):
                good = sum(1 for k in cuts if any((bits>>idx)&1 for idx,(x,y) in enumerate(tiles) if x>=k>y))
                by_good[good] += 1
            print(f"n={n}: by #good cuts: {dict(sorted(by_good.items()))}")
            print(f"  SC = {by_good[n-1]} (all n-1={n-1} cuts good)")
            print(f"  transitive = {by_good[0]} (0 cuts good)")
        else:
            # Just compute sc and transitive
            sc = 2**m - nonsc_ie(n)
            transitive = 1  # only 1: all tiles downward
            print(f"n={n}: SC={sc}, transitive=1 (all-downward tiling)")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 11: Diagonal odd-cycle sequence extension
# ─────────────────────────────────────────────────────────────────────────────

def directed_cycles_of_length(T, L):
    """Count directed L-cycles starting from smallest vertex in cycle."""
    verts = sorted(T.keys())
    count = 0
    def dfs(start, cur, path, path_set):
        nonlocal count
        if len(path) == L:
            if start in T[cur]:
                count += 1
            return
        if len(path) >= L:
            return
        for nxt in T[cur]:
            if nxt not in path_set and (len(path) == 1 or nxt > start):
                path.append(nxt)
                path_set.add(nxt)
                dfs(start, nxt, path, path_set)
                path.pop()
                path_set.discard(nxt)
    for start in verts:
        dfs(start, start, [start], {start})
    return count

def staircase_tournament(k):
    """
    The staircase tournament T_k: vertices 0..2k-1.
    T_k is the all-0 tiling (all tiles downward) = transitive on the natural ordering.
    Wait no — let me re-check. T_k from the paper is NOT the transitive tournament.
    From S5 (THM-322): T_k has 2k vertices and specific structure.
    The staircase T_k = the SPECIFIC tiling with all tiles at 0 (all base-path arcs).
    This is the tournament where i beats j if i > j (the transitive tournament).
    H(T_k) = staircase Hamiltonian path count (which is not 1 for k>1 ...).

    Actually from the session notes: T_k has H values 5,29,233,2489,33773,...
    That's not the transitive tournament (which has H=1).

    Let me look at what T_k actually is from THM-322.
    From the definitions file: T_k is the staircase on 2k+2 vertices?
    Or on k vertices? Need to be careful.

    From THM-322: k=2..8, H(T_k) = 5,29,233,2489,33773,562685,11222321,...
    From S23: k=2..12 sequence.
    The staircase T_k seems to have m=C(2k,2)/2 tiles or something.

    Actually: from the code in the previous sessions, T_k is built as a specific
    tournament and H(T_k) grows rapidly. It's the "all-0" staircase (all tiles upward?).
    Let me not pursue this here and focus on simpler sequences.
    """
    pass

def seq11_odd_cycle_seqs(max_n=6):
    """
    For each n and each odd length L, count total directed L-cycles
    summed over all labeled tournaments.
    """
    print("=" * 65)
    print("SEQUENCE 11: Odd Cycle Count Sequences")
    print("=" * 65)

    for n in range(3, max_n+1):
        odd_lens = range(3, n+1, 2)
        totals = {}
        for T in all_labeled_tournaments(n):
            for L in odd_lens:
                totals[L] = totals.get(L, 0) + directed_cycles_of_length(T, L)
        print(f"n={n}: sum of directed odd cycles over all {2**(n*(n-1)//2)} labeled tournaments:")
        for L in odd_lens:
            avg = totals.get(L, 0) / 2**(n*(n-1)//2)
            print(f"  L={L}: total={totals.get(L,0)}, avg per tournament={avg:.4f}")
    print()

    # Useful new sequence: average #3-cycles in a random n-tournament
    # Known: C(n,3)/4 (each 3-vertex subset has exactly one 3-cycle with prob 1/4)
    print("Average #3-cycles formula: C(n,3)/4")
    for n in range(3, 8):
        expected = n*(n-1)*(n-2)//6 / 4
        print(f"  n={n}: C(n,3)/4 = {n*(n-1)*(n-2)//24} (as fraction: {n*(n-1)*(n-2)}/24)")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 12: New OEIS candidate — SC-fraction numerators with 2^m denominator
# ─────────────────────────────────────────────────────────────────────────────

def seq12_sc_fraction_analysis(max_n=12):
    print("=" * 65)
    print("SEQUENCE 12: SC Tiling Fraction Numerators")
    print("=" * 65)
    # SC(n) / 2^m = fraction of tilings that are SC
    # SC(n) = 2^m - non-SC(n) where non-SC(n) = nonsc_ie(n)
    sc_nums = []
    for n in range(3, max_n+1):
        m = (n-1)*(n-2)//2
        nonsc = nonsc_ie(n)
        sc = (1 << m) - nonsc
        sc_nums.append(sc)
        frac_num = sc
        frac_den = 1 << m
        from math import gcd
        g = gcd(frac_num, frac_den)
        print(f"n={n}: SC={sc}/2^{m}={frac_num//g}/{frac_den//g}  "
              f"≈ {sc/frac_den:.8f}")
    print()
    print(f"SC sequence (n=3..{max_n-1}): {sc_nums}")

    # Check P(SC) → 1 as n → ∞
    probs = [(1<<((n-1)*(n-2)//2)) - nonsc_ie(n) for n in range(3, max_n+1)]
    total_probs = [1<<((n-1)*(n-2)//2) for n in range(3, max_n+1)]
    print("\nP(SC) convergence:")
    for i, n in enumerate(range(3, max_n+1)):
        p = probs[i]/total_probs[i]
        print(f"  n={n}: P(SC) = {p:.10f}  1-P(SC) = {(1-p):.2e}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 13: Non-SC IE formula — closed form attempt
# ─────────────────────────────────────────────────────────────────────────────

def seq13_nonsc_closed_form(max_n=10):
    print("=" * 65)
    print("SEQUENCE 13: Non-SC IE Formula — Closed Form Analysis")
    print("=" * 65)
    # The IE formula: |non-SC| = Σ_{∅≠S⊆{1..n-1}} (-1)^{|S|+1} 2^{m - f(S)}
    # where f(S) = #{tiles crossing any cut in S}
    #            = #{tiles (x,y): S ∩ {y+1,...,x} ≠ ∅}

    # Key observation: f(S) depends on which intervals S covers.
    # If S = {k}: f({k}) = k(n-k) - 1  (cut size formula)
    # If S = {j, k} with j < k:
    #   tiles crossing j OR k = tiles with {y+1,...,x} ∩ {j,k} ≠ ∅
    #   = tiles crossing j + tiles crossing k - tiles crossing BOTH j and k
    #   Tiles crossing both j and k: x >= k > y AND x >= j > y = x >= k, y < j
    #   Count: #{(x,y): x >= k, y < j, x >= y+2} = Σ_{y=0}^{j-1} (n-1-max(k,y+2)+1)
    #   For y < j < k, y+2 <= j+1 <= k, so x >= k gives count = n-k
    #   EXCEPT if y = j-1, y+2 = j+1, we need x >= k (since k > j >= j+1 if k>j).
    #   Wait, y < j, so y <= j-1, y+2 <= j+1 <= k (since k > j).
    #   So for all y < j: count of x >= k is n-1-k+1 = n-k.
    #   f({j,k}) both = j*(n-k) ... wait let me be careful.

    # Tiles crossing BOTH j and k (j < k):
    # Condition: x >= k > j > y (strictly), so x >= k, y <= j-1
    # Count: Σ_{y=0}^{j-1} #{x: x >= max(k, y+2) = k (since y+2 <= j+1 <= k)}
    #       = j * (n-1-k+1) = j*(n-k)
    # So f({j,k}) = f({j}) + f({k}) - j*(n-k)
    #             = (j(n-j)-1) + (k(n-k)-1) - j(n-k)
    #             = j(n-j) + k(n-k) - j(n-k) - 2
    #             = j(n-j-n+k) + k(n-k) - 2
    #             = j(k-j) + k(n-k) - 2

    print("f(S) formula for small S:")
    print("f({k}) = k(n-k) - 1")
    print("f({j,k}) for j<k: f({j}) + f({k}) - j*(n-k)")
    print("                = j(k-j) + k(n-k) - 2")
    print()

    # Verify
    for n in range(4, 8):
        tiles = get_tiles(n)
        cuts = list(range(1, n))
        print(f"n={n}: Verifying f formulas:")
        # Single cuts
        for k in cuts:
            actual = sum(1 for (x,y) in tiles if x >= k > y)
            formula = k*(n-k) - 1
            print(f"  f({{{k}}}) = {actual} (formula={formula}) {'✓' if actual==formula else '✗'}")
        # Pairs
        for j,k in itertools.combinations(cuts, 2):
            actual = len({idx for idx,(x,y) in enumerate(tiles) if (x>=j>y or x>=k>y)})
            formula = j*(k-j) + k*(n-k) - 2
            print(f"  f({{{j},{k}}}) = {actual} (formula={formula}) {'✓' if actual==formula else '✗'}")
        print()

    # Now: can we get a closed form for |non-SC|?
    # |non-SC| = Σ_{∅≠S} (-1)^{|S|+1} * 2^{m - f(S)}
    # For |S|=1: Σ_k 2^{m - f({k})} = Σ_k 2^{m-k(n-k)+1}
    # For |S|=2: -Σ_{j<k} 2^{m-f({j,k})}
    # etc.

    print("Size-1 IE term = Σ_k 2^{m-k(n-k)+1}:")
    for n in range(3, max_n+1):
        m = (n-1)*(n-2)//2
        term1 = sum(1<<(m-k*(n-k)+1) for k in range(1, n))
        print(f"  n={n}: {term1}")
    print()

    print("Full non-SC sequence via IE:")
    seq = []
    for n in range(3, max_n+1):
        seq.append(nonsc_ie(n))
    print(f"  n=3..{max_n-1}: {seq}")
    print()

    # Check: does log(non-SC(n)) grow as ~m - constant = C(n-1,2)?
    import math
    print("log2(non-SC) vs m=C(n-1,2):")
    for n in range(3, max_n+1):
        m = (n-1)*(n-2)//2
        ns = nonsc_ie(n)
        if ns > 0:
            print(f"  n={n}: m={m}, log2(non-SC)={math.log2(ns):.4f}, m-log2(non-SC)={m-math.log2(ns):.4f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 14: Search OEIS text API for our sequences
# ─────────────────────────────────────────────────────────────────────────────

def seq14_oeis_searches():
    import subprocess
    print("=" * 65)
    print("SEQUENCE 14: OEIS Searches for New Sequences")
    print("=" * 65)

    sequences_to_search = [
        ("Non-SC tiling counts n=3..10",
         "1,3,14,121,1995,64648,4064479,497219737"),
        ("SC tiling counts n=3..10",
         "1,5,50,903,30773,2032504,260862753,62972756903"),
        ("SC labeled tourn n=1..6",
         "1,0,2,24,544,"),  # partial
        ("H-impossible count n=3..?",
         "0,0,1,"),  # missing H values
    ]

    for name, seq in sequences_to_search:
        cmd = f'curl -s "https://oeis.org/search?q={seq.replace(",","%2C")}&fmt=text"'
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=15)
        output = result.stdout.strip()
        if "No results" in output:
            print(f"{name}: NOT in OEIS (new sequence!)")
        elif output:
            lines = [l for l in output.split('\n') if l.startswith('%') or 'A0' in l]
            print(f"{name}: FOUND - {output[:200]}")
        else:
            print(f"{name}: search failed")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# Sequence 15: Q-P gap distribution — #tournaments with each gap value
# ─────────────────────────────────────────────────────────────────────────────

def seq15_gap_distributions(max_n=6):
    print("=" * 65)
    print("SEQUENCE 15: Q-P Gap Distributions")
    print("=" * 65)
    for n in range(3, max_n+1):
        gap_ct = Counter()
        h_by_gap = defaultdict(list)
        for T in all_labeled_tournaments(n):
            degs = [len(T[v]) for v in T]
            g = max(degs) - min(degs)
            gap_ct[g] += 1
            h_by_gap[g].append(ham_paths(T))
        print(f"n={n}: Q-P gap distribution:")
        for g in sorted(gap_ct.keys()):
            hs = sorted(set(h_by_gap[g]))
            avg_h = sum(h_by_gap[g]) / len(h_by_gap[g])
            print(f"  gap={g}: {gap_ct[g]} tourn, H∈{hs}, avg_H={avg_h:.2f}")
    print()

# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("opus-2026-05-27-S2: Tournament Sequence Exploration")
    print("=" * 65)
    print()

    seq1_nonsc_tiling(max_n=11)         # Fast via IE formula up to n=11
    seq2_ie_structure(max_n=8)
    seq3_score_seqs(max_n=7)
    seq4_labeled_sc(max_n=6)
    seq5_h_distributions(max_n=5)       # Brute force up to n=5 only
    seq6_kings_distributions(max_n=6)
    seq7_sc_excess(max_n=6)
    seq8_h_by_score(max_n=5)
    seq9_cut_profiles(max_n=8)          # Fast via IE
    seq10_transitive_tilings()
    seq11_odd_cycle_seqs(max_n=5)
    seq12_sc_fraction_analysis(max_n=12)
    seq13_nonsc_closed_form(max_n=11)
    seq14_oeis_searches()
    seq15_gap_distributions(max_n=5)

    print("=" * 65)
    print("Session complete.")
