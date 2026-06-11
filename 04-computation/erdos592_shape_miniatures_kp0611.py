#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
THM-471 lab: the THM-460 B3 GENERAL-SHAPE recursive enumerator (POKE Task 1.2),
m=2 Schipperus-cutoff probes, and the first m=3 (open $1000 case) data.
kind-pasteur-2026-06-11-S1.

AMBIENT.  A(m; s, c) = functions [s]^m -> [c], ordered by largest disagreement,
where positions delta in [s]^m carry the exponent significance of the CNF
presentation of ordinals below w^(w^m) (THM-460 A): position delta = (d_0..d_{m-1})
has ordinal weight sum d_i * w^i, so d_{m-1} is most significant.  We linearize:
POS sorted by reversed-tuple lex DESCENDING, so POS[0] is the MOST significant
position; points are value-tuples over POS and tuple-lex order = ambient order
(preserving the m=1 invariant "list index order = ambient order").
At m=1 this reproduces ambient(s, C) of erdos592_chang_towers_macmini_s2 exactly.

SHAPES (tuple form of THM-460 B3; see THM-471 for the derivation from Lemma B2).
Exponent tuples delta in N^m (delta[i] = coefficient of w^i).  March convention
'j1' (faithful to THM-460 D's m=1 games: parts j = 1..M) or 'j0' (truncated march
j = 0..M-1, needed where j1 shapes outgrow the ambient — the m>=3 analogue of
THM-460 D's vacuousness guard).

  Bin(0) = single point.
  Finite peel  (least nonzero index 0):  two Bin(delta - e0) halves A < B, one
      common cross-split position p for ALL cross pairs, all internal splits
      strictly less significant than p.
  Limit peel  (least nonzero index i >= 1):  M order-separated parts
      Bin((delta - e_i) + j*(1,..,1) on positions 0..i-1),  j marching.
  BT(m, M) = M order-separated parts Bin((j,..,j)), j marching.

GAME.  Q^(m)(M; s, c): triangle-free graph on A(m; s, c) with no independent
BT(m, M).  THM-460 C1: SAT for all (s,c) at some fixed M => w^(w^m) -/-> (.,3)^2
strong; C2: positive relations (m=1 Chang, m=2 Schipperus) FORCE finite cutoffs.

VALIDATION (MISTAKE-067 discipline): the structural finder is cross-validated
against a brute-force enumerator on random graphs; at m=1 the checker is
cross-validated against is_binary_grid of the s2 script, and the M=2 j1 games
reproduce THM-460 D's SAT cells; M=1 reproduces the R(1,2)=3 Ramsey shadow.
"""

import itertools, time, random
from functools import lru_cache
from pysat.solvers import Glucose3

from erdos592_chang_towers_macmini_s2 import is_binary_grid


# ----------------------------------------------------------------- ambient

def positions(m, s):
    """[s]^m sorted MOST significant first (reversed-tuple lex, descending)."""
    return sorted(itertools.product(range(s), repeat=m),
                  key=lambda d: tuple(reversed(d)), reverse=True)


def ambient_points(m, s, c):
    """All functions [s]^m -> [c] as value-tuples over positions(m,s);
    tuple-lex order == largest-disagreement ambient order."""
    np = s ** m
    return list(itertools.product(range(c), repeat=np))


def split_pos(u, v):
    """Most significant disagreement = least index (positions are most-sig-first)."""
    for k in range(len(u)):
        if u[k] != v[k]:
            return k
    return None


# ----------------------------------------------------------------- shapes

def march_range(M, conv):
    return range(1, M + 1) if conv == 'j1' else range(0, M)


def peel(delta):
    """(kind, sub) where kind = 'fin' (i0 = 0) or ('lim', i0)."""
    i0 = next((i for i, d in enumerate(delta) if d > 0), None)
    if i0 is None:
        return ('leaf', None)
    sub = list(delta)
    sub[i0] -= 1
    if i0 == 0:
        return ('fin', tuple(sub))
    return (('lim', i0), tuple(sub))


def march_shape(sub, i0, j):
    out = list(sub)
    for i in range(i0):
        out[i] += j
    return tuple(out)


def shape_size(delta, M, conv):
    @lru_cache(maxsize=None)
    def sz(d):
        kind, sub = peel(d)
        if kind == 'leaf':
            return 1
        if kind == 'fin':
            return 2 * sz(sub)
        _, i0 = kind
        return sum(sz(march_shape(sub, i0, j)) for j in march_range(M, conv))
    return sz(delta)


def bt_part_list(m, M, conv):
    """Part exponents INCLUDING the j=0 single point under j0."""
    out = []
    for j in march_range(M, conv):
        out.append(tuple(j for _ in range(m)))
    return out


def bt_size(m, M, conv):
    return sum(shape_size(d, M, conv) for d in bt_part_list(m, M, conv))


# ------------------------------------------------------------- the checker

def is_bin_shape(delta, S, M, conv):
    """S = tuple of vectors sorted in ambient order. Checks the B3 grammar."""
    kind, sub = peel(delta)
    if kind == 'leaf':
        return len(S) == 1
    if kind == 'fin':
        n = len(S)
        if n % 2:
            return False
        A, B = S[:n // 2], S[n // 2:]
        p = split_pos(A[0], B[0])
        if p is None:
            return False
        for a in A:
            for b in B:
                if split_pos(a, b) != p:
                    return False
        for half in (A, B):
            for x, y in itertools.combinations(half, 2):
                sp = split_pos(x, y)
                if sp is None or sp <= p:
                    return False
        return is_bin_shape(sub, A, M, conv) and is_bin_shape(sub, B, M, conv)
    _, i0 = kind
    sizes = [shape_size(march_shape(sub, i0, j), M, conv) for j in march_range(M, conv)]
    if len(S) != sum(sizes):
        return False
    off = 0
    for j, szj in zip(march_range(M, conv), sizes):
        chunk = S[off:off + szj]
        off += szj
        if not is_bin_shape(march_shape(sub, i0, j), chunk, M, conv):
            return False
    return True
    # order-separation of chunks is automatic: S is sorted and chunks are consecutive


def is_bt(m, M, conv, S):
    """S sorted tuple of vectors: is it a BT(m, M)?"""
    sizes = [shape_size(d, M, conv) for d in bt_part_list(m, M, conv)]
    if len(S) != sum(sizes):
        return False
    off = 0
    for d, szj in zip(bt_part_list(m, M, conv), sizes):
        chunk = S[off:off + szj]
        off += szj
        if not is_bin_shape(d, chunk, M, conv):
            return False
    return True


# ------------------------------------------------------------ the finder

class BudgetExceeded(Exception):
    pass


class ShapeFinder:
    """Complete DFS for an independent BT(m,M) in a graph on the ambient.
    Structural recursion proposes, is_bin_shape disposes — TowerFinder
    architecture (THM-460 D), shape-generic, with two additions:
    (i) CONSTRUCTIVE cross-split pruning in the finite-peel case: choose the
        cross-split position p first (more significant than A's internal
        splits), then restrict the B-half's candidates to points agreeing with
        A[0] above p and differing at p.  Complete: a sorted Bin's second half
        lies in exactly that set (all of A and all of B agree at positions < p,
        and at p each half is constant with B's value above A's).
    (ii) an optional node budget: on explosion raises BudgetExceeded so callers
        report an honest TIMEOUT (never SAT) — MISTAKE-067 discipline."""

    def __init__(self, m, M, conv, L):
        self.m, self.M, self.conv = m, M, conv
        self.L = L
        self.N = len(L)
        self.npos = len(L[0])
        self.full = (1 << self.N) - 1
        # agree[k][i] = bitmask of points equal to point i at position k
        byval = [{} for _ in range(self.npos)]
        for i, v in enumerate(L):
            for k in range(self.npos):
                byval[k].setdefault(v[k], 0)
                byval[k][v[k]] |= 1 << i
        self.eqmask = [[byval[k][L[i][k]] for k in range(self.npos)]
                       for i in range(self.N)]

    def set_graph(self, nb, budget=None, anchor0=False):
        """anchor0: restrict the search to shapes whose MINIMUM is index 0.
        Complete ONLY for vertex-transitive graphs whose automorphisms act
        transitively while preserving the shape family — e.g. XOR-measurable
        rule graphs on {0,1}^npos (split positions are XOR-invariant and the
        fin-halves relabel under translation). Do NOT use on general graphs."""
        self.nb = nb  # list of int bitmasks
        self.budget = budget
        self.anchor0 = anchor0
        self.nodes = 0

    def _tick(self):
        self.nodes += 1
        if self.budget is not None and self.nodes > self.budget:
            raise BudgetExceeded

    def min_internal_split(self, A):
        """Most significant (= least index) split among pairs of A; npos if |A|<2."""
        best = self.npos
        for x, y in itertools.combinations(A, 2):
            sp = split_pos(self.L[x], self.L[y])
            if sp is not None and sp < best:
                best = sp
        return best

    def gen(self, delta, lo, taken_mask, allowed):
        """Yield sorted index-tuples forming an independent Bin(delta), indices
        > lo, none adjacent to taken_mask, all inside the `allowed` bitmask."""
        self._tick()
        kind, sub = peel(delta)
        if kind == 'leaf':
            if self.anchor0 and lo == -1 and taken_mask == 0:
                # the very first leaf of the whole search = the shape's minimum
                if allowed & 1:
                    yield (0,)
                return
            cand = allowed >> (lo + 1)
            i = lo + 1
            while cand:
                if cand & 1 and (self.nb[i] & taken_mask) == 0:
                    yield (i,)
                cand >>= 1
                i += 1
            return
        if kind == 'fin':
            for A in self.gen(sub, lo, taken_mask, allowed):
                amask = taken_mask
                for a in A:
                    amask |= 1 << a
                top = self.min_internal_split(A)
                a0 = A[0]
                # cross-split position p must be MORE significant (< top);
                # B agrees with A[0] at every position < p and differs at p.
                # Try p LEAST-significant-first: inner structures should consume
                # low significance, leaving high positions for outer levels
                # (most-significant-first starves outer levels and explodes).
                cands = []
                prefix = self.full
                for p in range(0, top):
                    cands.append(allowed & prefix & ~self.eqmask[a0][p])
                    prefix &= self.eqmask[a0][p]
                for p in range(top - 1, -1, -1):
                    b_allowed = cands[p]
                    if not b_allowed:
                        continue
                    for B in self.gen(sub, A[-1], amask, b_allowed):
                        S = A + B
                        if is_bin_shape(delta, tuple(self.L[i] for i in S),
                                        self.M, self.conv):
                            yield S
            return
        _, i0 = kind
        shapes = [march_shape(sub, i0, j) for j in march_range(self.M, self.conv)]

        def rec(k, lo2, mask):
            if k == len(shapes):
                yield ()
                return
            for P in self.gen(shapes[k], lo2, mask, allowed):
                pmask = mask
                for p in P:
                    pmask |= 1 << p
                for rest in rec(k + 1, P[-1], pmask):
                    yield P + rest
        yield from rec(0, lo, taken_mask)

    def find_bt(self):
        """Returns a tuple of indices (found), None (none exists — complete),
        or raises BudgetExceeded (budget tripped — NOT a certificate)."""
        parts = bt_part_list(self.m, self.M, self.conv)

        def rec(k, lo, mask):
            if k == len(parts):
                return ()
            for P in self.gen(parts[k], lo, mask, self.full):
                pmask = mask
                for p in P:
                    pmask |= 1 << p
                rest = rec(k + 1, P[-1], pmask)
                if rest is not None:
                    return P + rest
            return None
        return rec(0, -1, 0)


# ----------------------------------------------------- brute-force cross-check

def brute_bt(m, M, conv, L, nb):
    """Enumerate all index-combos of total BT size, filter independence +
    chunk-structure (completeness: parts are order-separated, so the sorted
    combo's consecutive-chunk partition is forced)."""
    total = bt_size(m, M, conv)
    N = len(L)
    if total > N:
        return None
    for combo in itertools.combinations(range(N), total):
        ok = True
        for a, b in itertools.combinations(combo, 2):
            if (nb[a] >> b) & 1:
                ok = False
                break
        if not ok:
            continue
        if is_bt(m, M, conv, tuple(L[i] for i in combo)):
            return combo
    return None


def crossvalidate():
    print("=== cross-validation: ShapeFinder vs brute force on random graphs ===", flush=True)
    rng = random.Random(23)
    cases = [(1, 2, 'j1', 2, 3), (1, 2, 'j1', 3, 2), (2, 1, 'j1', 2, 2),
             (2, 2, 'j0', 2, 2), (2, 1, 'j1', 2, 3), (3, 1, 'j1', 2, 2)]
    for (m, M, conv, s, c) in cases:
        L = ambient_points(m, s, c)
        N = len(L)
        if bt_size(m, M, conv) > N:
            print(f"   m={m} M={M} {conv} (s,c)=({s},{c}): shape {bt_size(m,M,conv)} > N={N}, skip", flush=True)
            continue
        trials = 20 if N <= 100 else 6
        finder = ShapeFinder(m, M, conv, L)
        bad = 0
        for tr in range(trials):
            p = 0.35 * rng.random()
            nb = [0] * N
            for i in range(N):
                for j in range(i + 1, N):
                    if rng.random() < p:
                        nb[i] |= 1 << j
                        nb[j] |= 1 << i
            finder.set_graph(nb)
            r1 = finder.find_bt()
            if N <= 100 or bt_size(m, M, conv) <= 4:
                r2 = brute_bt(m, M, conv, L, nb)
                if (r1 is None) != (r2 is None):
                    bad += 1
                    print(f"   DISAGREE m={m} M={M} {conv} ({s},{c}) tr={tr}: "
                          f"struct={r1 is not None} brute={r2 is not None}", flush=True)
            if r1 is not None:
                for a, b in itertools.combinations(r1, 2):
                    assert not (nb[a] >> b) & 1, "found shape not independent!"
                assert is_bt(m, M, conv, tuple(L[i] for i in r1)), "found shape fails checker!"
        print(f"   m={m} M={M} {conv} (s,c)=({s},{c}) N={N} |BT|={bt_size(m,M,conv)}: "
              f"{trials} trials {'OK' if bad == 0 else str(bad) + ' BAD'}", flush=True)


def validate_m1_against_old():
    print("\n=== m=1 checker agreement with is_binary_grid (s2 script) ===", flush=True)
    # at m=1, Bin((h)) = binary h-grid; compare on all 4-subsets of [3]^... use s=3,c=3
    L = ambient_points(1, 3, 3)   # functions [3]->[3] = triples, 27 points
    mism = 0
    for S in itertools.combinations(L, 4):
        a = is_bin_shape((2,), tuple(sorted(S)), 2, 'j1')
        b = is_binary_grid(tuple(sorted(S)))
        if a != b:
            mism += 1
    print(f"   all 4-subsets of [3]^3 as binary 2-grids: {mism} mismatches "
          f"({'OK' if mism == 0 else 'BAD'})", flush=True)


# ------------------------------------------------------------------ the game

def solve_shape_game(m, M, conv, s, c, tlimit=1200, verbose=True, batch=40,
                     node_budget=3_000_000):
    L = ambient_points(m, s, c)
    N = len(L)
    size = bt_size(m, M, conv)
    if size > N:
        print(f"   VACUOUS m={m} M={M} {conv} ({s},{c}): |BT|={size} > N={N}", flush=True)
        return 'vacuous'
    ev = {}
    cnt = 0
    for i in range(N):
        for j in range(i + 1, N):
            cnt += 1
            ev[(i, j)] = cnt
    sol = Glucose3()
    lazy_tri = N > 130
    if not lazy_tri:
        for i, j, k in itertools.combinations(range(N), 3):
            sol.add_clause([-ev[(i, j)], -ev[(i, k)], -ev[(j, k)]])
    finder = ShapeFinder(m, M, conv, L)
    t0 = time.time()
    shapes_blocked = tris_blocked = 0
    while True:
        if time.time() - t0 > tlimit:
            if verbose:
                print(f"   TIMEOUT m={m} M={M} {conv} ({s},{c}) "
                      f"(shapes={shapes_blocked}, tris={tris_blocked}, {time.time()-t0:.1f}s)", flush=True)
            return None
        if not sol.solve():
            if verbose:
                print(f"   UNSAT  m={m} M={M} {conv} ({s},{c}) N={N} "
                      f"(shapes={shapes_blocked}, tris={tris_blocked}, {time.time()-t0:.1f}s)", flush=True)
            return False
        model = set(l for l in sol.get_model() if l > 0)
        nb = [0] * N
        edges = []
        for (i, j), v in ev.items():
            if v in model:
                nb[i] |= 1 << j
                nb[j] |= 1 << i
                edges.append((i, j))
        if lazy_tri:
            # block triangles present in the model (batch)
            found_tri = 0
            for i, j in edges:
                common = nb[i] & nb[j]
                if common:
                    k = common.bit_length() - 1
                    a, b, cc = sorted((i, j, k))
                    sol.add_clause([-ev[(a, b)], -ev[(a, cc)], -ev[(b, cc)]])
                    tris_blocked += 1
                    found_tri += 1
                    if found_tri >= batch:
                        break
            if found_tri:
                continue
        finder.set_graph(nb, budget=node_budget)
        try:
            bad = finder.find_bt()
        except BudgetExceeded:
            if verbose:
                print(f"   TIMEOUT m={m} M={M} {conv} ({s},{c}) — finder node budget "
                      f"exceeded (shapes={shapes_blocked}, tris={tris_blocked}, "
                      f"{time.time()-t0:.1f}s) [NOT a SAT certificate]", flush=True)
            return None
        if bad is None:
            if verbose:
                print(f"   SAT    m={m} M={M} {conv} ({s},{c}) N={N} ({len(edges)} edges, "
                      f"shapes={shapes_blocked}, tris={tris_blocked}, {time.time()-t0:.1f}s)", flush=True)
            return True
        sol.add_clause([ev[(min(a, b), max(a, b))] for a, b in itertools.combinations(bad, 2)])
        shapes_blocked += 1


def main():
    t0 = time.time()
    print("=== shape sizes (the miniature-design table) ===", flush=True)
    for m in (1, 2, 3):
        for M in (1, 2):
            for conv in ('j1', 'j0'):
                print(f"   m={m} M={M} {conv}: |BT| = {bt_size(m, M, conv)}", flush=True)

    crossvalidate()
    validate_m1_against_old()

    print("\n=== m=1 calibration: M=1 = Ramsey shadow R(1,2)=3 ===", flush=True)
    for (s, c) in ((1, 2), (1, 3), (2, 2)):
        solve_shape_game(1, 1, 'j1', s, c, tlimit=120)

    print("\n=== m=1 calibration: M=2 j1 reproduces THM-460 D SAT cells ===", flush=True)
    for (s, c) in ((2, 3), (3, 2), (2, 4)):
        solve_shape_game(1, 2, 'j1', s, c, tlimit=600)

    print("\n=== m=2 probes (Schipperus positive => cutoffs FORCED, THM-460 C2) ===", flush=True)
    print("   --- M=1 j1 (shape = Bin((1,1)), %d leaves) ---" % bt_size(2, 1, 'j1'), flush=True)
    for (s, c) in ((2, 2), (2, 3), (2, 4), (2, 5), (3, 2), (2, 6)):
        r = solve_shape_game(2, 1, 'j1', s, c, tlimit=1500)
        if r is False:
            print(f"   *** FIRST m=2 M=1 CUTOFF at (s,c)=({s},{c}) ***", flush=True)
    print("   --- M=2 j0 (truncated march, %d leaves) ---" % bt_size(2, 2, 'j0'), flush=True)
    for (s, c) in ((2, 2), (2, 3), (2, 4), (3, 2)):
        r = solve_shape_game(2, 2, 'j0', s, c, tlimit=1500)
        if r is False:
            print(f"   *** m=2 M=2(j0) CUTOFF at (s,c)=({s},{c}) ***", flush=True)

    print("\n=== m=3 PROBES of the open case (first data ever) ===", flush=True)
    print("   --- M=1 j1 (shape = Bin((1,1,1)), %d leaves) ---" % bt_size(3, 1, 'j1'), flush=True)
    for (s, c) in ((2, 2),):
        solve_shape_game(3, 1, 'j1', s, c, tlimit=2400)
    print("   --- M=2 j0 (%d leaves) ---" % bt_size(3, 2, 'j0'), flush=True)
    for (s, c) in ((2, 2),):
        solve_shape_game(3, 2, 'j0', s, c, tlimit=2400)

    print(f"\n=== TOTAL {time.time()-t0:.1f}s ===", flush=True)


if __name__ == "__main__":
    main()
