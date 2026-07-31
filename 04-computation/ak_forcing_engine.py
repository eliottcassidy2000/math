#!/usr/bin/env python3
"""
ak_forcing_engine.py — klein-S691 (2026-07-28)

Exact engine for the *arithmetic Kakeya* forcing-graph benchmark
("constructible graph + forcing pair" certificates for AK(alpha)).

BENCHMARK (external, verbatim semantics implemented here):
  X ⊂ Z^2 finite with (0,0) ∈ X and a+b ≠ 0 for all nonzero (a,b) ∈ X.
  A constructible graph G is given by dims d_1..d_k and labels
  f_i : d_1×…×d_{i-1}×(d_i−1) → X.  Vertices V = d_1×…×d_k (1-indexed
  tuples), n(G) = ∏ d_i,  m(G) = Σ_i Σ_prefix 1[f_i(prefix)≠0]·d_{i+1}···d_k.
  An X-set is (R, T): T ⊆ V, R = list of singleton-support functions V→X.
  Forcing operations on (R,T):
   (1) edge rule — for e1 with i-prefix (a_1..a_i) and e2 with i-prefix
       (a_1..a_i+1), where x = f_i(a_1..a_i) ≠ 0: add [e1↦x, e2↦−x] to R.
   (2) if some f ∈ ⟨R⟩ has support ⊆ T ∪ {e} and f(e) = (a,−a), a≠0:
       add e to T.
   (3) closure of R under differences (⟹ under the generated Z-lattice).
  Forcing succeeds iff T reaches V.  Score = (m+r)/(n−t), with r = |R_initial|
  and t = |T_initial|.  A valid certificate of score s is claimed to prove
  AK(s); trivial baseline s = 2; the open target is s ≤ 1.675 (the current
  human record is the Katz–Tao 2002 exponent, the root 1.67513… of
  x^3 − 4x + 2 = 0; any certificate ≤ 1.675 would be a new theorem, IF the
  benchmark's soundness theorem for this certificate format is correct).

TWO READINGS OF RULE (1) — both implemented:
  * LOOSE (the literal spec text): e1, e2 need only MATCH ON THE FIRST i
    COORDINATES (prefix p and p+); their suffixes are unconstrained.
  * STRICT (the intuitive/glued-graph reading): e1, e2 must additionally
    share the SAME suffix — the classical "k+1 copies of H glued along
    corresponding vertices" picture, edges connect corresponding vertices.
  These are NOT equivalent: loose-mode derives within-layer transporters
    x(δ_{(p,s)} − δ_{(p,s')})  =  x(δ_{(p,s)} − δ_{(p+,s'')}) − x(δ_{(p,s')} − δ_{(p+,s'')})
  which strict mode cannot produce.  See ak_exploit_cert.py: loose mode
  admits a trivial family with score → 1, i.e. the literal spec text is
  UNSOUND as an AK certificate format (it would "prove" the full arithmetic
  Kakeya conjecture).  The strict game is the mathematically honest target.

FORCING = LATTICE-SECTION GAME (implementation form):
  Maintain the lattice L ⊆ (Z^2)^V spanned by all seed vectors x·δ_v and all
  rule-(1) edge vectors.  A vertex v ∉ T fires iff L contains a vector
  supported on T ∪ {v} whose v-column is a nonzero multiple of (1,−1).
  Scale-invariance makes the Q-span test equivalent (clear denominators),
  so exact rational elimination decides firing.  T-columns are wildcards:
  project them out and ask for vectors supported on the single pair of
  v-columns.  Completeness of the pivot-row test: any w in the row span of
  an RREF basis is w = Σ w[pivot_i]·row_i, so supp(w) ⊆ pair(v) forces every
  contributing row's pivot to lie in pair(v); it suffices to inspect rows
  whose pivot column belongs to pair(v) and which vanish elsewhere.

This file is the ENGINE + VERIFIER + baselines only.  Searches live in
ak_strict_search.py; the loose-mode exploit lives in ak_exploit_cert.py.

Status discipline: everything here is FINITE-EXACT computation; the AK(s)
consequence of a certificate is CITED from the external benchmark statement,
not re-proved here.
"""

from fractions import Fraction
from itertools import product
import ast, sys

Vec2 = tuple  # (a, b) integers


# ----------------------------------------------------------------------
# instance
# ----------------------------------------------------------------------

class AKInstance:
    def __init__(self, X, dims, fs, T0, R):
        """
        X    : list of (a,b) int pairs, must contain (0,0)
        dims : list of positive ints [d_1..d_k]
        fs   : list of k dicts; fs[i0] maps prefix tuples
               (e_1..e_{i0+1}) with 1<=e_j<=d_j (j<=i0), 1<=e_{i0+1}<=d_{i0+1}-1
               to elements of X.  Omitted entries mean (0,0).
        T0   : list of vertex tuples (initial T)
        R    : list of (vertex_tuple, x) seed pairs, x in X
        """
        self.X = [tuple(x) for x in X]
        self.dims = list(dims)
        self.fs = [dict(f) for f in fs]
        self.T0 = [tuple(v) for v in T0]
        self.R = [(tuple(v), tuple(x)) for (v, x) in R]
        self.k = len(dims)
        self.verts = list(product(*[range(1, d + 1) for d in dims]))
        self.vid = {v: i for i, v in enumerate(self.verts)}
        self.n = len(self.verts)

    # ---- validation against the benchmark's constraints ----
    def validate(self):
        errs = []
        if (0, 0) not in self.X:
            errs.append("X must contain (0,0)")
        for (a, b) in self.X:
            if (a, b) != (0, 0) and a + b == 0:
                errs.append(f"X element {(a,b)} has a+b=0")
        if any(d < 1 for d in self.dims):
            errs.append("dims must be >= 1")
        if len(self.fs) != self.k:
            errs.append("need one f_i per dimension")
        for i0, f in enumerate(self.fs):
            for pre, x in f.items():
                if len(pre) != i0 + 1:
                    errs.append(f"f_{i0+1} key {pre} has wrong length")
                    continue
                ok = all(1 <= pre[j] <= self.dims[j] for j in range(i0)) \
                     and 1 <= pre[i0] <= self.dims[i0] - 1
                if not ok:
                    errs.append(f"f_{i0+1} key {pre} out of domain")
                if tuple(x) not in self.X:
                    errs.append(f"f_{i0+1}[{pre}]={x} not in X")
        for v in self.T0:
            if v not in self.vid:
                errs.append(f"T0 vertex {v} not a vertex")
        for (v, x) in self.R:
            if v not in self.vid:
                errs.append(f"R seed vertex {v} not a vertex")
            if x not in self.X:
                errs.append(f"R seed value {x} not in X")
            if x == (0, 0):
                errs.append("R seed value (0,0) is not singleton-support")
        return errs

    # ---- cost parameters ----
    def m(self):
        tot = 0
        for i0, f in enumerate(self.fs):
            sfx = 1
            for d in self.dims[i0 + 1:]:
                sfx *= d
            tot += sum(sfx for x in f.values() if tuple(x) != (0, 0))
        return tot

    def r(self):
        return len(self.R)

    def t(self):
        return len(set(self.T0))

    def score(self):
        den = self.n - self.t()
        return Fraction(self.m() + self.r(), den) if den > 0 else None

    # ---- generator vectors of the lattice L ----
    def generators(self, mode):
        """Sparse vectors: dict vid -> (a,b).  mode in {'strict','loose'}."""
        gens = []
        for (v, x) in self.R:
            gens.append({self.vid[v]: (x[0], x[1])})
        for i0, f in enumerate(self.fs):
            suffixes = list(product(*[range(1, d + 1)
                                      for d in self.dims[i0 + 1:]]))
            for pre, x in f.items():
                x = tuple(x)
                if x == (0, 0):
                    continue
                prep = pre[:i0] + (pre[i0] + 1,)
                if mode == "strict":
                    for s in suffixes:
                        u, w = self.vid[pre + s], self.vid[prep + s]
                        gens.append({u: x, w: (-x[0], -x[1])})
                elif mode == "loose":
                    # span-equal basis for {x(δ_{(pre,s)}−δ_{(prep,s')})}:
                    s0 = suffixes[0]
                    u0, w0 = self.vid[pre + s0], self.vid[prep + s0]
                    gens.append({u0: x, w0: (-x[0], -x[1])})
                    for s in suffixes[1:]:
                        u = self.vid[pre + s]
                        gens.append({u: x, u0: (-x[0], -x[1])})
                        w = self.vid[prep + s]
                        gens.append({w: x, w0: (-x[0], -x[1])})
                else:
                    raise ValueError(mode)
        return gens


# ----------------------------------------------------------------------
# exact forcing closure
# ----------------------------------------------------------------------

def _rref(rows, ncols):
    """In-place fraction RREF; rows = list of lists of Fraction. Returns
    list of (pivot_col, row) for nonzero rows."""
    piv = []
    rr = 0
    for c in range(ncols):
        p = None
        for i in range(rr, len(rows)):
            if rows[i][c] != 0:
                p = i
                break
        if p is None:
            continue
        rows[rr], rows[p] = rows[p], rows[rr]
        inv = Fraction(1, 1) / rows[rr][c]
        rows[rr] = [e * inv for e in rows[rr]]
        for i in range(len(rows)):
            if i != rr and rows[i][c] != 0:
                fac = rows[i][c]
                rows[i] = [a - fac * b for a, b in zip(rows[i], rows[rr])]
        piv.append((c, rr))
        rr += 1
        if rr == len(rows):
            break
    return [(c, rows[i]) for (c, i) in piv]


def _pair_span_hits_11(vecs):
    """vecs: list of (a,b) Fractions.  Does span_Q contain (1,-1)?"""
    vecs = [v for v in vecs if v != (0, 0) and (v[0] != 0 or v[1] != 0)]
    if not vecs:
        return False
    # rank 2 -> everything; rank 1 -> need parallel to (1,-1)
    v0 = vecs[0]
    for v in vecs[1:]:
        if v0[0] * v[1] - v0[1] * v[0] != 0:
            return True
    return v0[0] * (-1) - v0[1] * 1 == 0  # v0 ∥ (1,−1)


def force(inst, mode, verbose=False, max_rounds=None):
    """Run the forcing closure.  Returns (success, T_final, trace)."""
    gens = inst.generators(mode)
    T = set(inst.vid[v] for v in inst.T0)
    n = inst.n
    trace = []
    rounds = 0
    while len(T) < n:
        rounds += 1
        if max_rounds and rounds > max_rounds:
            break
        free = [v for v in range(n) if v not in T]
        cols = []
        colof = {}
        for v in free:
            for c in (2 * v, 2 * v + 1):
                colof[c] = len(cols)
                cols.append(c)
        # dense projected rows (T-columns dropped = wildcards)
        rows = []
        for g in gens:
            row = [Fraction(0)] * len(cols)
            live = False
            for v, (a, b) in g.items():
                if v in T:
                    continue
                row[colof[2 * v]] = Fraction(a)
                row[colof[2 * v + 1]] = Fraction(b)
                live = live or a or b
            if live:
                rows.append(row)
        if not rows:
            break
        piv = _rref(rows, len(cols))
        fired = []
        # group pivot rows by vertex; a candidate vector supported on
        # pair(v) must be a combo of rows with pivot in pair(v) vanishing
        # elsewhere (completeness argued in module docstring)
        byv = {}
        for (c, row) in piv:
            v = cols[c] // 2
            byv.setdefault(v, []).append(row)
        ncols = len(cols)
        for v, rws in byv.items():
            j0, j1 = colof[2 * v], colof[2 * v + 1]
            ok_vecs = []
            if len(rws) == 1:
                r1 = rws[0]
                if all(r1[j] == 0 for j in range(ncols) if j != j0 and j != j1):
                    ok_vecs.append((r1[j0], r1[j1]))
            else:  # exactly two pivot rows can meet pair(v)
                r1, r2 = rws
                outside = [j for j in range(ncols)
                           if j != j0 and j != j1
                           and (r1[j] != 0 or r2[j] != 0)]
                if not outside:
                    ok_vecs = [(r1[j0], r1[j1]), (r2[j0], r2[j1])]
                else:
                    # nullspace of {(a,b): a·r1[j] + b·r2[j] = 0, j outside}
                    base, rank2 = None, False
                    for j in outside:
                        p, q = r1[j], r2[j]
                        if p != 0 or q != 0:
                            if base is None:
                                base = (p, q)
                            elif base[0] * q - base[1] * p != 0:
                                rank2 = True
                                break
                    if not rank2 and base is not None:
                        a, b = -base[1], base[0]
                        vec = (a * r1[j0] + b * r2[j0],
                               a * r1[j1] + b * r2[j1])
                        if vec != (0, 0):
                            ok_vecs.append(vec)
            if _pair_span_hits_11(ok_vecs):
                fired.append(v)
        if not fired:
            break
        for v in fired:
            T.add(v)
        trace.append(sorted(inst.verts[v] for v in fired))
        if verbose:
            print(f"  round {rounds}: fired {len(fired)} "
                  f"({[inst.verts[v] for v in fired][:6]}...)")
    return len(T) == n, T, trace


# ----------------------------------------------------------------------
# submission format (6 lines)
# ----------------------------------------------------------------------

def emit(inst, mode_note=""):
    s = inst.score()
    l1 = (f"score={s} m={inst.m()} r={inst.r()} n={inst.n} t={inst.t()}"
          + (f"  [{mode_note}]" if mode_note else ""))
    l2 = repr(sorted(set(inst.X)))
    l3 = repr(list(inst.dims))
    l4 = repr([{k: v for k, v in sorted(f.items())} for f in inst.fs])
    l5 = repr([tuple(v) for v in inst.T0])
    l6 = repr([{tuple(v): tuple(x)} for (v, x) in inst.R])
    return "\n".join([l1, l2, l3, l4, l5, l6])


def parse(text):
    lines = [ln for ln in text.strip().splitlines() if ln.strip()]
    X = ast.literal_eval(lines[1])
    dims = ast.literal_eval(lines[2])
    fs = ast.literal_eval(lines[3])
    T0 = ast.literal_eval(lines[4])
    Rd = ast.literal_eval(lines[5])
    R = []
    for d in Rd:
        (v, x), = d.items()
        R.append((v, x))
    return AKInstance(X, dims, fs, T0, R)


def verify(text_or_inst, mode, expect_score=None, name=""):
    inst = parse(text_or_inst) if isinstance(text_or_inst, str) else text_or_inst
    errs = inst.validate()
    if errs:
        print(f"[{name}] INVALID: {errs}")
        return False
    ok, T, trace = force(inst, mode)
    sc = inst.score()
    status = "FORCED" if ok else f"STUCK at |T|={len(T)}/{inst.n}"
    print(f"[{name}] mode={mode:6s} {status}  score={sc} "
          f"(m={inst.m()} r={inst.r()} n={inst.n} t={inst.t()})")
    if expect_score is not None and ok:
        assert sc == expect_score, (sc, expect_score)
    return ok


# ----------------------------------------------------------------------
# baselines
# ----------------------------------------------------------------------

def baseline_trivial():
    """Single vertex, two independent seeds: the trivial AK(2) proof."""
    X = [(0, 0), (1, 0), (0, 1)]
    return AKInstance(X, [1], [{}], [], [((1,), (1, 0)), ((1,), (0, 1))])


def baseline_rail_grid(L):
    """Strict-mode L×L rook-rail grid, score exactly 2 (validates strict
    telescoping dynamics: monochromatic rails deliver one direction each)."""
    X = [(0, 0), (1, 0), (0, 1)]
    f1 = {(i,): (0, 1) for i in range(1, L)}            # column rails (axis 1)
    f2 = {(i, j): (1, 0) for i in range(1, L + 1)       # row rails (axis 2)
          for j in range(1, L)}
    R = [((1, 1), (1, 0)), ((1, 1), (0, 1))]
    R += [((1, j), (0, 1)) for j in range(2, L + 1)]    # row 1 seeds
    R += [((i, 1), (1, 0)) for i in range(2, L + 1)]    # column 1 seeds
    return AKInstance(X, [L, L], [f1, f2], [], R)


if __name__ == "__main__":
    print("== ak_forcing_engine self-test ==")
    tv = baseline_trivial()
    verify(tv, "strict", Fraction(2), "trivial")
    verify(tv, "loose", Fraction(2), "trivial")
    for L in (3, 4):
        g = baseline_rail_grid(L)
        verify(g, "strict", Fraction(2), f"rail-grid L={L}")
    print("roundtrip:", parse(emit(tv)).score() == Fraction(2))
    print("OK")
