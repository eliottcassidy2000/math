#!/usr/bin/env python3
"""unknot1_decider.py -- three-valued decider core for "unknotting number 1".

Input: a 1-indexed planar diagram (PD) code, JSON list of 4-tuples
[a,b,c,d] per crossing (KnotTheory convention: a = incoming under-strand,
then b,c,d counterclockwise; c = outgoing under-strand; over-strand runs
b->d or d->b, reconstructed robustly from global in/out consistency of the
arc-successor structure, with the label-successor mod-2n rule as fallback).

Verdicts
  TRUE_CERTIFIED  -- an explicit crossing change is exhibited whose changed
                     diagram is reduced to the 0-crossing diagram by
                     face-local Reidemeister R1/R2 moves, AND the input is
                     independently certified nontrivial.  The latter condition
                     is essential: a crossing change from one unknot diagram
                     to another only proves u(K) <= 1, not u(K) = 1.
  FALSE_CERTIFIED -- a named obstruction proves u(K) >= 2:
                     (i)  Murasugi: |sigma(K)| <= 2 u(K), so |sigma| >= 4;
                     (ii) H1(Sigma(K)) must be cyclic if u = 1 (Montesinos);
                     (iii) Lickorish linking-form: u = 1 forces a generator g
                          of H1(Sigma(K)) = Z/det with lambda(g,g) = +-2/det.
                     Obstruction (iii) is DECISIVE only if the built-in gate
                     suite passes (trefoil/fig8 pass it, 7_4 must fail it);
                     otherwise it is report-only.
  UNKNOWN         -- neither side certified. This is expected sometimes: a
                     complete u=1 decision procedure is research-level (the
                     unknotting crossing may exist only in another diagram;
                     a full method needs bounded Reidemeister exploration
                     and deeper obstructions, e.g. Ozsvath-Szabo
                     d-invariants).

Method: checkerboard faces are rebuilt from the PD rotation system (at each
crossing the four arc-ends a,b,c,d are counterclockwise; faces = orbits of
the next-corner map).  Goeritz matrix over each color class, signature by
exact integer/Fraction congruence diagonalization, Gordon-Litherland
correction mu = sum of eta over type-II crossings.  The two sign
conventions (ETA_SIGN, TYPE_T0) are calibrated once against Knot Atlas
targets (3_1: sigma=-2 det 3; mirror: +2; 4_1: 0/5; 7_4: -2/15) and the
computation is cross-checked by requiring BOTH colorings to give the same
sigma and det on every call.

TODO (descoped this pass): Seifert-matrix / Alexander-polynomial route
(Seifert circles + band linking numbers; Delta(t) = det(V - tV^T)); would
sharpen the unknot screen (Delta == 1) and give genus heuristics.
"""

import json
import sys
from fractions import Fraction
from math import gcd

# Calibration constants, fixed by `--calibrate` against Knot Atlas targets.
ETA_SIGN = 1  # global sign of the Goeritz incidence eta
TYPE_T0 = 1   # parity bit of the type-II rule (see goeritz section)
# Calibrated 2026-07-28 via --calibrate: this combo gives 3_1 -> sigma -2,
# det 3 (Knot Atlas convention), mirror 3_1 -> +2, 4_1 -> 0/5, P(3,1,3)
# [= 7_4 up to mirror] -> |sigma| 2, det 15; it is the unique combo with
# trefoil sigma == -2.

VERDICT_TRUE = "TRUE_CERTIFIED"
VERDICT_FALSE = "FALSE_CERTIFIED"
VERDICT_UNKNOWN = "UNKNOWN"


def require(condition, message):
    """Optimization-stable internal validity check."""
    if not condition:
        raise ArithmeticError(message)


SCOPE_NOTES = (
    "Scope: UNKNOWN is possible; deciding u(K)=1 in general is "
    "research-level. TRUE needs (i) an independently certified nontrivial "
    "input, (ii) the unknotting change to be visible in THIS diagram, and "
    "(iii) the changed diagram to simplify monotonically under face-local "
    "R1/R2 (no R3 / no crossing-increasing exploration). A visible change "
    "without (i) certifies only u(K)<=1 and is reported as UNKNOWN. FALSE "
    "uses Murasugi |sigma|<=2u, cyclicity of H1 of the double branched "
    "cover, and the Lickorish linking-form test (gated by a 7_4 self-test). "
    "Deeper tools (Ozsvath-Szabo d-invariants, Alexander/Seifert route) are "
    "not implemented."
)

TREFOIL = [[1, 4, 2, 5], [3, 6, 4, 1], [5, 2, 6, 3]]
# figure-8 PD derived from the braid closure of (s1 s2^-1)^2 (the PD
# floating around in the task brief fails arc-count validation).
FIG8 = [[3, 8, 4, 1], [1, 7, 2, 6], [7, 4, 8, 5], [5, 3, 6, 2]]
# 7_4 as the pretzel P(3,1,3), traced by hand (the PD quoted in the task
# brief is not planar under the CCW convention).  A reduced alternating
# 7-crossing diagram with det 15 is 7_4, uniquely among 7-crossing knots;
# det is verified at runtime by the gate suite.
K7_4 = [[10, 2, 11, 1], [2, 10, 3, 9], [8, 4, 9, 3], [4, 12, 5, 11],
        [14, 6, 1, 5], [6, 14, 7, 13], [12, 8, 13, 7]]
EXAMPLE11 = [[4, 2, 5, 1], [10, 6, 11, 5], [8, 3, 9, 4], [2, 9, 3, 10],
             [11, 16, 12, 17], [7, 15, 8, 14], [15, 7, 16, 6],
             [13, 20, 14, 21], [17, 22, 18, 1], [21, 18, 22, 19],
             [19, 12, 20, 13]]
# An exact hostile control for the distinction u(K)=1 versus u(K)<=1.
# This is an unknot diagram obtained from the one-crossing unknot by the
# legal reverse-Reidemeister path R1+, R2+, R2+, R3, R3.  Greedy R1/R2
# stalls on the input, but changing crossing #4 makes greedy R1/R2 reach
# the empty diagram.  The unrepaired engine therefore returned a false
# TRUE_CERTIFIED for this u=0 knot.
HOSTILE_UNKNOT6 = [
    [1, 11, 2, 10], [6, 10, 7, 9], [3, 8, 4, 9],
    [11, 5, 12, 4], [7, 2, 8, 3], [5, 1, 6, 12],
]

# SECTION: linalg
# Exact linear algebra over Z / Q (Fractions only; no floats anywhere).


def _frac_copy(M):
    return [[Fraction(x) for x in row] for row in M]


def det_exact(M):
    """Signed determinant of an integer matrix, exact (returns int)."""
    m = len(M)
    if m == 0:
        return 1
    A = _frac_copy(M)
    det = Fraction(1)
    for k in range(m):
        piv = next((i for i in range(k, m) if A[i][k] != 0), None)
        if piv is None:
            return 0
        if piv != k:
            A[k], A[piv] = A[piv], A[k]
            det = -det
        det *= A[k][k]
        for i in range(k + 1, m):
            f = A[i][k] / A[k][k]
            if f:
                A[i] = [A[i][j] - f * A[k][j] for j in range(m)]
    require(det.denominator == 1, "exact determinant became nonintegral")
    return int(det)


def inverse_exact(M):
    """Exact inverse (Fraction entries) of a nonsingular matrix, or None."""
    m = len(M)
    A = _frac_copy(M)
    I = [[Fraction(int(i == j)) for j in range(m)] for i in range(m)]
    for k in range(m):
        piv = next((i for i in range(k, m) if A[i][k] != 0), None)
        if piv is None:
            return None
        A[k], A[piv] = A[piv], A[k]
        I[k], I[piv] = I[piv], I[k]
        p = A[k][k]
        A[k] = [x / p for x in A[k]]
        I[k] = [x / p for x in I[k]]
        for i in range(m):
            if i != k and A[i][k]:
                f = A[i][k]
                A[i] = [A[i][j] - f * A[k][j] for j in range(m)]
                I[i] = [I[i][j] - f * I[k][j] for j in range(m)]
    return I


def adjugate_int(M):
    """(adjugate as integer matrix, signed det).  M integer, det != 0.
    adj(M) = det(M) * M^{-1}; entries are provably integers."""
    m = len(M)
    d = det_exact(M)
    if m == 0:
        return [], 1
    if d == 0:
        return None, 0
    inv = inverse_exact(M)
    adj = []
    for i in range(m):
        row = []
        for j in range(m):
            v = inv[i][j] * d
            require(v.denominator == 1, "adjugate not integral -- bug")
            row.append(int(v))
        adj.append(row)
    return adj, d


def sym_signature_exact(M):
    """Signature (#pos - #neg pivots) of a symmetric integer matrix by
    exact congruence diagonalization over Q."""
    m = len(M)
    A = _frac_copy(M)
    pos = neg = 0
    for k in range(m):
        if A[k][k] == 0:
            j = next((i for i in range(k + 1, m) if A[i][i] != 0), None)
            if j is not None:  # symmetric swap of rows&cols k,j
                A[k], A[j] = A[j], A[k]
                for t in range(m):
                    A[t][k], A[t][j] = A[t][j], A[t][k]
            else:
                j = next((i for i in range(k + 1, m) if A[k][i] != 0), None)
                if j is None:
                    continue  # zero block: contributes 0 (degenerate)
                for t in range(m):  # add row/col j into k: pivot 2*A[k][j]
                    A[k][t] += A[j][t]
                for t in range(m):
                    A[t][k] += A[t][j]
        p = A[k][k]
        if p > 0:
            pos += 1
        else:
            neg += 1
        for i in range(k + 1, m):
            if A[i][k]:
                f = A[i][k] / p
                for t in range(m):
                    A[i][t] -= f * A[k][t]
                for t in range(m):
                    A[t][i] -= f * A[t][k]
    return pos - neg

# SECTION: diagram


class Diagram(object):
    """Internal combinatorial diagram.
    crossings: cid -> [eid,eid,eid,eid]  (slots 0..3, CCW; slot0 = under-in,
               slot2 = under-out, slots 1,3 = over strand)
    edges:     eid -> [(cid,slot),(cid,slot)]  (its two ends)
    heads:     (cid,slot) -> True iff the edge comes IN at this slot
    free_loops: number of closed 0-crossing circles produced by reduction."""

    def __init__(self):
        self.crossings = {}
        self.edges = {}
        self.heads = {}
        self.free_loops = 0

    @property
    def n(self):
        return len(self.crossings)

    def over_in_slot(self, c):
        return 1 if self.heads[(c, 1)] else 3

    def sign(self, c):
        # Calibrated: over b->d (over-in at slot 1) is a NEGATIVE crossing
        # (Knot Atlas 3_1 PD has writhe -3, 4_1 writhe 0 under this rule).
        return -1 if self.over_in_slot(c) == 1 else 1

    def copy(self):
        d = Diagram()
        d.crossings = {c: list(v) for c, v in self.crossings.items()}
        d.edges = {e: list(v) for e, v in self.edges.items()}
        d.heads = dict(self.heads)
        d.free_loops = self.free_loops
        return d


def parse_pd(pd):
    """Parse and validate a 1-indexed PD code; return a Diagram."""
    if not isinstance(pd, list) or any(
            not isinstance(t, list) or len(t) != 4 or
            any(not isinstance(x, int) or x < 1 for x in t) for t in pd):
        raise ValueError("PD must be a list of 4-lists of positive ints")
    n = len(pd)
    dg = Diagram()
    if n == 0:
        dg.free_loops = 1
        return dg
    occ = {}
    for c, t in enumerate(pd):
        for s, lab in enumerate(t):
            occ.setdefault(lab, []).append((c, s))
    if set(occ) != set(range(1, 2 * n + 1)):
        raise ValueError("arc labels must be exactly 1..2n (2n=%d)" % (2 * n))
    for lab, ps in occ.items():
        if len(ps) != 2:
            raise ValueError("arc label %d appears %d times (need exactly 2)"
                             % (lab, len(ps)))
    # roles: (c,slot) -> True (incoming) / False (outgoing)
    roles = {}
    for c in range(n):
        roles[(c, 0)] = True
        roles[(c, 2)] = False

    def other_occ(lab, here):
        a, b = occ[lab]
        return b if a == here else a

    # propagate over-strand orientations from global one-in/one-out per arc
    changed = True
    while changed:
        changed = False
        for c, t in enumerate(pd):
            if (c, 1) in roles:
                continue
            for s in (1, 3):
                oc = other_occ(t[s], (c, s))
                if oc in roles:
                    r = not roles[oc]
                    roles[(c, s)] = r
                    roles[(c, (s + 2) % 4)] = not r
                    changed = True
                    break
    for c, t in enumerate(pd):  # fallback: sequential-label successor rule
        if (c, 1) in roles:
            continue
        b, d = t[1], t[3]
        if d == b % (2 * n) + 1:
            roles[(c, 1)], roles[(c, 3)] = True, False
        elif b == d % (2 * n) + 1:
            roles[(c, 1)], roles[(c, 3)] = False, True
        else:
            raise ValueError("cannot orient over-strand at crossing %d" % c)
    for lab, ps in occ.items():
        rr = sorted(roles[p] for p in ps)
        if rr != [False, True]:
            raise ValueError("arc %d is not once-in/once-out (invalid PD)"
                             % lab)
    # successor permutation along the knot; must be a single 2n-cycle
    succ = {}
    for c, t in enumerate(pd):
        succ[t[0]] = t[2]
        if roles[(c, 1)]:
            succ[t[1]] = t[3]
        else:
            succ[t[3]] = t[1]
    seen, x = set(), 1
    for _ in range(2 * n):
        if x in seen:
            break
        seen.add(x)
        x = succ[x]
    if len(seen) != 2 * n or x != 1:
        raise ValueError("PD is not a knot (successor structure is not a "
                         "single 2n-cycle); links are out of scope")
    for c, t in enumerate(pd):
        dg.crossings[c] = list(t)
    for lab, ps in occ.items():
        dg.edges[lab] = list(ps)
    dg.heads = roles
    return dg

# SECTION: faces
# Corner (c,s) = the quadrant between arc-end slots s and s+1 (mod 4) at
# crossing c.  Faces are orbits of the next-corner map; checkerboard 2-color.


def _other_end(dg, eid, here):
    a, b = dg.edges[eid]
    if a == here:
        return b
    if b == here:
        return a
    raise AssertionError("edge/end inconsistency")


def build_faces(dg, strict=True):
    """Return (faces, cornerface): faces = list of corner-lists;
    cornerface[(c,s)] = face index.  Validates Euler count n+2 if strict."""
    corners = [(c, s) for c in dg.crossings for s in range(4)]
    cornerface = {}
    faces = []
    for c0 in corners:
        if c0 in cornerface:
            continue
        face = []
        cur = c0
        for _ in range(4 * dg.n + 1):
            face.append(cur)
            cornerface[cur] = len(faces)
            c, s = cur
            leave = (c, (s + 1) % 4)
            cur = _other_end(dg, dg.crossings[c][(s + 1) % 4], leave)
            if cur == c0:
                break
        else:
            raise ValueError("face walk failed to close (invalid diagram)")
        faces.append(face)
    if strict and dg.n > 0 and len(faces) != dg.n + 2:
        raise ValueError("face count %d != n+2 = %d: diagram not a planar "
                         "connected knot shadow" % (len(faces), dg.n + 2))
    return faces, cornerface


def checkerboard(dg, faces, cornerface):
    """2-color the faces (0/1) so faces adjacent across an edge differ."""
    nf = len(faces)
    color = [None] * nf
    adj = [set() for _ in range(nf)]
    for c in dg.crossings:
        for s in range(4):
            fa = cornerface[(c, s)]
            fb = cornerface[(c, (s - 1) % 4)]
            adj[fa].add(fb)
            adj[fb].add(fa)
    stack = [0]
    color[0] = 0
    while stack:
        f = stack.pop()
        for g in adj[f]:
            if color[g] is None:
                color[g] = 1 - color[f]
                stack.append(g)
            elif color[g] == color[f]:
                raise ValueError("diagram is not checkerboard colorable")
    if any(c is None for c in color):
        raise ValueError("face adjacency graph disconnected")
    return color

# SECTION: goeritz
# Goeritz matrix + Gordon-Litherland correction.
# At crossing c the corner pairs {0,2} and {1,3} get opposite colors.
# eta(c) = ETA_SIGN * (+1 if corners {1,3} are white else -1)   (unoriented)
# type II(c) = (corners13 white) XOR (over-in at slot 1) XOR TYPE_T0
# sigma(K) = sig(Goeritz_white) - sum_{type II} eta(c); computed for BOTH
# colorings and required to agree (structural cross-check of conventions).


def _goeritz_one_color(dg, faces, cornerface, color, white):
    wf = [i for i in range(len(faces)) if color[i] == white]
    idx = {f: i for i, f in enumerate(wf)}
    m = len(wf)
    G = [[0] * m for _ in range(m)]
    mu = 0
    for c in dg.crossings:
        p13white = (color[cornerface[(c, 1)]] == white)
        eta = ETA_SIGN * (1 if p13white else -1)
        bd = (dg.over_in_slot(c) == 1)
        if (p13white != bd) != bool(TYPE_T0):
            mu += eta
        ca, cb = ((c, 1), (c, 3)) if p13white else ((c, 0), (c, 2))
        fi, fj = idx[cornerface[ca]], idx[cornerface[cb]]
        if fi != fj:
            G[fi][fj] -= eta
            G[fj][fi] -= eta
    for i in range(m):
        G[i][i] = -sum(G[i][j] for j in range(m) if j != i)
    Gr = [row[1:] for row in G[1:]]  # delete face wf[0]
    dsigned = det_exact(Gr)
    sig = sym_signature_exact(Gr)
    return {"G": Gr, "det_signed": dsigned, "det": abs(dsigned),
            "sig": sig, "mu": mu, "sigma": sig - mu, "m": len(Gr)}


def goeritz_invariants(dg):
    """Return dict with sigma, det, and the (white) Goeritz matrix, after
    cross-checking that both checkerboard colorings agree."""
    if dg.n == 0:
        return {"sigma": 0, "det": 1, "G": [], "det_signed": 1,
                "dual_ok": True, "mu": 0, "sig": 0, "m": 0}
    faces, cornerface = build_faces(dg)
    color = checkerboard(dg, faces, cornerface)
    r0 = _goeritz_one_color(dg, faces, cornerface, color, 0)
    r1 = _goeritz_one_color(dg, faces, cornerface, color, 1)
    if r0["sigma"] != r1["sigma"] or r0["det"] != r1["det"]:
        raise ArithmeticError(
            "Gordon-Litherland dual-coloring mismatch (sigma %s vs %s, det "
            "%s vs %s): convention bug" %
            (r0["sigma"], r1["sigma"], r0["det"], r1["det"]))
    out = dict(r0)
    out["dual_ok"] = True
    return out


def writhe(dg):
    return sum(dg.sign(c) for c in dg.crossings)

# SECTION: lickorish
# H1(Sigma(K)) = Z^m / G Z^m, |H1| = det(K).  If u(K) = 1 then (Montesinos)
# Sigma(K) is +-det/2 surgery on a knot, so H1 is CYCLIC and (Lickorish,
# "The unknotting number of a classical knot") the linking form
# lambda(x,y) = +-(G^{-1})_{xy} mod 1 has a generator g with
# lambda(g,g) = +-2/det.  Testing both signs makes the result immune to the
# global sign convention of lambda.  Cyclic <=> gcd of adjugate entries = 1
# (the (m-1)-st determinant divisor D_{m-1}; d_{m-1} = D_{m-1}/D_{m-2} > 1
# iff D_{m-1} > 1 since all earlier invariant factors divide d_{m-1}).

LICKORISH_D_CAP = 10 ** 6


def linking_form_analysis(G):
    """Analyze the linking form presented by integer symmetric G.
    Returns dict: applicable, cyclic, d, q, obstructed, detail."""
    m = len(G)
    d0 = det_exact(G)
    d = abs(d0)
    base = {"applicable": False, "cyclic": None, "d": d, "q": None,
            "obstructed": False, "detail": ""}
    if d == 0:
        base["detail"] = "det 0 (not a knot Goeritz matrix)"
        return base
    if d == 1:
        base.update(cyclic=True, detail="det 1: H1(Sigma) trivial, "
                    "linking-form test vacuous")
        return base
    if d % 2 == 0:
        base["detail"] = "det even -- not a knot? test skipped"
        return base
    if d > LICKORISH_D_CAP:
        base["detail"] = "det too large for exhaustive generator loop"
        return base
    adj, _ = adjugate_int(G)
    gall = 0
    for row in adj:
        for v in row:
            gall = gcd(gall, abs(v))
    cyclic = (gcd(gall, d) == 1) if m > 1 else True
    if not cyclic:
        return {"applicable": True, "cyclic": False, "d": d, "q": None,
                "obstructed": True,
                "detail": "H1(Sigma(K)) is NON-CYCLIC (gcd of adjugate "
                          "entries %d > 1): u=1 impossible (Montesinos)"
                          % gcd(gall, d)}
    # find a generator: x with gcd(d, gcd(adj x)) = 1
    cands = [[int(i == k) for i in range(m)] for k in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            for sj in (1, -1):
                v = [0] * m
                v[i], v[j] = 1, sj
                cands.append(v)
    gen = None
    for x in cands:
        ax = [sum(adj[i][j] * x[j] for j in range(m)) for i in range(m)]
        c = 0
        for v in ax:
            c = gcd(c, abs(v))
        if gcd(c, d) == 1:
            gen = (x, ax)
            break
    if gen is None:
        base["detail"] = "cyclic, but no generator found in candidate set"
        base["cyclic"] = True
        return base
    x, ax = gen
    q = sum(x[i] * ax[i] for i in range(m)) % d
    targets = {2 % d, (-2) % d}
    ok = any((q * t * t) % d in targets
             for t in range(1, d) if gcd(t, d) == 1)
    return {"applicable": True, "cyclic": True, "d": d, "q": q,
            "obstructed": not ok,
            "detail": "linking form q/d = %d/%d; generator self-linkings "
                      "{q*t^2 mod d} %s contain +-2"
                      % (q, d, "DO" if ok else "do NOT")}


_GATE_CACHE = {}


def lickorish_gate_ok():
    """Decisive only if: trefoil & fig8 PASS the test, 7_4 FAILS it
    (Lickorish's classic u(7_4) >= 2 example, det 15, form 4/15)."""
    if "ok" in _GATE_CACHE:
        return _GATE_CACHE["ok"]
    ok = True
    try:
        for pd, d_expect, obstructed_expect in (
                (TREFOIL, 3, False), (FIG8, 5, False), (K7_4, 15, True)):
            inv = goeritz_invariants(parse_pd(pd))
            lk = linking_form_analysis(inv["G"])
            if (inv["det"] != d_expect or not lk["applicable"]
                    or lk["obstructed"] != obstructed_expect):
                ok = False
    except Exception:
        ok = False
    _GATE_CACHE["ok"] = ok
    return ok

# SECTION: reduction
# Greedy face-local Reidemeister simplification.  R1 = monogon face; R2 =
# bigon face on two distinct crossings whose common strand is over at both
# (equivalently under at both).  Moves are applied by "dissolving" the
# removed crossings: each removed crossing joins its slot pairs {0,2} and
# {1,3} into straight wires; edge chains through wires are merged, and
# chains that close up become free (0-crossing) unknotted circles.
# Face-locality makes both moves genuine planar Reidemeister moves.


def dissolve(dg, removed):
    removed = set(removed)
    wires = {}
    for c in removed:
        for s in range(4):
            wires[(c, s)] = (c, (s + 2) % 4)
    visited = set()

    def walk(eid, entry_end):
        """Follow the strand from edge eid entered at entry_end's far side;
        return (list of edge ids, final end) or (chain, None) if closed."""
        chain = [eid]
        cur_far = _other_end(dg, eid, entry_end)
        for _ in range(4 * len(dg.edges) + 4):
            if cur_far not in wires:
                return chain, cur_far
            nxt_port = wires[cur_far]
            nxt = dg.crossings[nxt_port[0]][nxt_port[1]]
            chain.append(nxt)
            cur_far = _other_end(dg, nxt, nxt_port)
            if nxt == chain[0] and len(chain) > 1:
                return chain[:-1], None  # closed loop
        raise AssertionError("dissolve walk did not terminate")

    # 1) chains anchored at surviving crossings
    for eid in list(dg.edges):
        if eid in visited:
            continue
        ends = dg.edges[eid]
        surv = [p for p in ends if p[0] not in removed]
        if len(surv) == 2:
            continue  # untouched edge
        if len(surv) == 1:
            start = surv[0]
            chain, final = walk(eid, start)
            require(final is not None and final[0] not in removed,
                    "dissolve chain did not end at a surviving crossing")
            visited.update(chain)
            keep = chain[0]
            for e in chain[1:]:
                if e in dg.edges:
                    del dg.edges[e]
            dg.edges[keep] = [start, final]
            dg.crossings[final[0]][final[1]] = keep
    # 2) leftover all-removed chains close into free loops
    for eid in list(dg.edges):
        if eid in visited or eid not in dg.edges:
            continue
        ends = dg.edges[eid]
        if ends[0][0] in removed and ends[1][0] in removed:
            chain, final = walk(eid, ends[0])
            for e in set(chain):
                if e in dg.edges:
                    del dg.edges[e]
            dg.free_loops += 1
    for c in removed:
        del dg.crossings[c]
        for s in range(4):
            dg.heads.pop((c, s), None)


def find_moves(dg):
    """Return ('R1', crossing) or ('R2', (X, Y)) or None."""
    faces, cornerface = build_faces(dg, strict=False)
    for f in faces:
        if len(f) == 1:
            return ("R1", f[0][0])
    for f in faces:
        if len(f) == 2:
            (X, sx), (Y, sy) = f
            if X == Y:
                continue
            e1 = dg.crossings[X][(sx + 1) % 4]
            e2 = dg.crossings[X][sx]
            if e1 == e2:
                continue
            # e1 over at X iff slot (sx+1) odd; over at Y iff slot sy odd
            if ((sx + 1) % 2) == (sy % 2):
                return ("R2", (X, Y))
    return None


def try_reduce(dg):
    """Greedy R1/R2 to exhaustion.  Returns (reduced_diagram, log, success)
    with success True iff the empty (0-crossing, 1-circle) diagram reached."""
    d = dg.copy()
    log = []
    while d.crossings:
        try:
            mv = find_moves(d)
        except Exception as ex:  # face build failed -> treat as stall
            log.append("stall(%s)" % ex)
            return d, log, False
        if mv is None:
            break
        if mv[0] == "R1":
            dissolve(d, [mv[1]])
            log.append("R1@%d" % mv[1])
        else:
            dissolve(d, list(mv[1]))
            log.append("R2@(%d,%d)" % mv[1])
    success = (not d.crossings) and d.free_loops == 1
    return d, log, success

# SECTION: search


def flip_crossing(dg, c):
    """Crossing change at c (swap over/under), keeping planar rotation:
    rotate the slot list so the old over-in slot becomes slot 0."""
    d = dg.copy()
    k = d.over_in_slot(c)  # 1 or 3
    old = d.crossings[c]
    d.crossings[c] = [old[(i + k) % 4] for i in range(4)]
    oldheads = {s: d.heads[(c, s)] for s in range(4)}
    for i in range(4):
        d.heads[(c, i)] = oldheads[(i + k) % 4]
    for eid in set(old):
        d.edges[eid] = [(cc, (ss - k) % 4) if cc == c else (cc, ss)
                        for (cc, ss) in d.edges[eid]]
    require(d.heads[(c, 0)] and not d.heads[(c, 2)],
            "crossing flip did not restore the PD under-strand convention")
    return d


def mirror(dg):
    d = dg
    for c in list(dg.crossings):
        d = flip_crossing(d, c)
    return d


def unknotting_search(dg):
    """Try every single crossing change; certify unknot via det==1 screen +
    greedy R1/R2 reduction to the empty diagram."""
    results = []
    certified = None
    for c in sorted(dg.crossings):
        entry = {"crossing": c, "status": "", "log": None}
        fd = flip_crossing(dg, c)
        try:
            inv = goeritz_invariants(fd)
        except Exception as ex:
            entry["status"] = "invariant error: %s" % ex
            results.append(entry)
            continue
        if inv["det"] != 1:
            entry["status"] = "screened out (det=%d != 1)" % inv["det"]
            results.append(entry)
            continue
        if inv["sigma"] != 0:
            entry["status"] = ("screened out (det=1 but sigma=%d != 0)"
                               % inv["sigma"])
            results.append(entry)
            continue
        red, log, ok = try_reduce(fd)
        if ok:
            entry["status"] = "UNKNOT CERTIFIED (reduced to 0 crossings)"
            entry["log"] = log
            certified = entry
            results.append(entry)
            break
        entry["status"] = ("probably-unknot-uncertified (det=1, sigma=0, "
                           "reduction stalled at %d crossings)" % red.n)
        entry["log"] = log
        results.append(entry)
    return certified, results

# SECTION: cli


def decide(pd):
    """Full pipeline.  Returns a result dict (see report())."""
    res = {"pd": pd, "verdict": VERDICT_UNKNOWN, "certificate": None,
           "upper_bound_certificate": None,
           "nontriviality_certificates": [],
           "obstructions": [], "notes": [], "invariants": {},
           "search": None}
    dg = parse_pd(pd)
    inv = goeritz_invariants(dg)
    res["invariants"] = {"n": dg.n, "det": inv["det"],
                         "sigma": inv["sigma"], "writhe": writhe(dg)}
    # Is the input itself already the unknot?  (then u=0, so u != 1)
    red, log, is_unknot = try_reduce(dg)
    if is_unknot:
        res["verdict"] = VERDICT_FALSE
        res["obstructions"].append(
            "input reduces to the 0-crossing diagram: K is the UNKNOT, "
            "u(K)=0 != 1 (R1/R2 certificate: %s)" % ",".join(log) if log
            else "input is the 0-crossing unknot: u(K)=0 != 1")
        return res
    # A stalled simplifier is not a certificate that the knot is nontrivial.
    # Use only invariants whose value differs from that of the unknot.
    if inv["det"] != 1:
        res["nontriviality_certificates"].append(
            "det(K)=%d != 1=det(unknot)" % inv["det"])
    if inv["sigma"] != 0:
        res["nontriviality_certificates"].append(
            "sigma(K)=%+d != 0=sigma(unknot)" % inv["sigma"])
    # obstruction battery
    if abs(inv["sigma"]) >= 4:
        res["obstructions"].append(
            "Murasugi: |sigma(K)| = %d >= 4 and |sigma| <= 2u ==> u >= 2"
            % abs(inv["sigma"]))
    lk = linking_form_analysis(inv["G"])
    res["lickorish"] = lk
    gate = lickorish_gate_ok()
    if lk["applicable"] and lk["obstructed"]:
        msg = ("Lickorish linking-form obstruction: %s ==> u >= 2"
               % lk["detail"])
        if gate:
            res["obstructions"].append(msg)
        else:
            res["notes"].append("REPORT-ONLY (gate suite failed): " + msg)
    elif lk["detail"]:
        res["notes"].append("linking form: " + lk["detail"])
    # crossing-change search
    certified, sresults = unknotting_search(dg)
    res["search"] = sresults
    if certified is not None and res["obstructions"]:
        raise ArithmeticError("CONTRADICTION: certified unknotting change "
                              "AND a u>=2 obstruction -- bug somewhere")
    if res["obstructions"]:
        res["verdict"] = VERDICT_FALSE
    elif certified is not None:
        change_certificate = (
            "change crossing #%d (1-indexed), then reduce: %s"
            % (certified["crossing"] + 1, ",".join(certified["log"])))
        if res["nontriviality_certificates"]:
            res["verdict"] = VERDICT_TRUE
            res["certificate"] = (
                "%s; input nontrivial because %s"
                % (change_certificate,
                   "; ".join(res["nontriviality_certificates"])))
        else:
            res["verdict"] = VERDICT_UNKNOWN
            res["upper_bound_certificate"] = change_certificate
            res["notes"].append(
                "the crossing change certifies u(K)<=1, but exact u(K)=1 "
                "is WITHHELD: greedy failure does not certify that the input "
                "is nontrivial; use an exact unknot oracle or another "
                "nontriviality certificate")
    else:
        res["verdict"] = VERDICT_UNKNOWN
        if any("probably-unknot" in e["status"] for e in sresults):
            res["notes"].append(
                "some crossing change passed all invariant screens but "
                "greedy R1/R2 stalled: candidate unknotting change is "
                "plausible but UNCERTIFIED")
    return res


def report(res, out=sys.stdout):
    w = out.write
    inv = res["invariants"]
    w("crossings: %d   writhe: %+d\n" % (inv["n"], inv["writhe"]))
    w("det(K)  = %d\n" % inv["det"])
    w("sigma(K)= %+d   (Goeritz + Gordon-Litherland, dual-coloring "
      "cross-checked)\n" % inv["sigma"])
    w("Alexander polynomial: NOT COMPUTED (Seifert-matrix route descoped "
      "this pass)\n")
    w("obstructions (u >= 2):\n")
    if res["obstructions"]:
        for ob in res["obstructions"]:
            w("  * %s\n" % ob)
    else:
        w("  (none fired)\n")
    if res["search"] is not None:
        w("crossing-change search (%d sites tried):\n" % len(res["search"]))
        for e in res["search"]:
            w("  #%d: %s\n" % (e["crossing"] + 1, e["status"]))
    w("VERDICT: %s\n" % res["verdict"])
    if res["certificate"]:
        w("certificate: %s\n" % res["certificate"])
    if res.get("upper_bound_certificate"):
        w("upper-bound certificate only: %s\n"
          % res["upper_bound_certificate"])
    if res.get("nontriviality_certificates"):
        w("input nontriviality certificates:\n")
        for cert in res["nontriviality_certificates"]:
            w("  * %s\n" % cert)
    for nt in res["notes"]:
        w("note: %s\n" % nt)
    w("scope: %s\n" % SCOPE_NOTES)

def selftest(verbose=True):
    """Sanity suite; returns number of failures."""
    fails = 0

    def check(name, cond, detail=""):
        nonlocal fails
        ok = bool(cond)
        fails += 0 if ok else 1
        if verbose:
            print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, detail))

    tre = decide(TREFOIL)
    ti = tre["invariants"]
    check("trefoil det=3", ti["det"] == 3, "(got %d)" % ti["det"])
    check("trefoil |sigma|=2", abs(ti["sigma"]) == 2,
          "(got %+d)" % ti["sigma"])
    check("trefoil sigma=-2 (Knot Atlas 3_1 convention)", ti["sigma"] == -2)
    check("trefoil verdict TRUE_CERTIFIED", tre["verdict"] == VERDICT_TRUE,
          "(got %s)" % tre["verdict"])
    mtre = goeritz_invariants(mirror(parse_pd(TREFOIL)))
    check("mirror trefoil sigma=+2", mtre["sigma"] == 2,
          "(got %+d)" % mtre["sigma"])
    f8 = decide(FIG8)
    fi = f8["invariants"]
    check("fig8 det=5", fi["det"] == 5, "(got %d)" % fi["det"])
    check("fig8 sigma=0", fi["sigma"] == 0, "(got %+d)" % fi["sigma"])
    check("fig8 verdict TRUE_CERTIFIED", f8["verdict"] == VERDICT_TRUE,
          "(got %s)" % f8["verdict"])
    s74 = decide(K7_4)
    si = s74["invariants"]
    check("7_4 det=15", si["det"] == 15, "(got %d)" % si["det"])
    check("7_4 |sigma|=2", abs(si["sigma"]) == 2, "(got %+d)" % si["sigma"])
    check("7_4 Lickorish gate", lickorish_gate_ok())
    check("7_4 linking form obstructed", s74["lickorish"]["obstructed"],
          "(%s)" % s74["lickorish"]["detail"])
    check("7_4 verdict FALSE_CERTIFIED", s74["verdict"] == VERDICT_FALSE,
          "(got %s)" % s74["verdict"])
    hostile = decide(HOSTILE_UNKNOT6)
    hi = hostile["invariants"]
    check("hostile unknot det=1 sigma=0",
          hi["det"] == 1 and hi["sigma"] == 0,
          "(got det=%d sigma=%+d)" % (hi["det"], hi["sigma"]))
    _, hostile_log, hostile_reduced = try_reduce(parse_pd(HOSTILE_UNKNOT6))
    check("hostile input greedy R1/R2 stalls",
          not hostile_reduced and not hostile_log,
          "(log=%s)" % hostile_log)
    check("hostile u=0 diagram is NOT TRUE_CERTIFIED",
          hostile["verdict"] == VERDICT_UNKNOWN,
          "(got %s)" % hostile["verdict"])
    check("hostile retains only u<=1 certificate",
          hostile["upper_bound_certificate"] is not None,
          "(certificate=%s)" % hostile["upper_bound_certificate"])
    return fails


def calibrate():
    """Grid-search ETA_SIGN/TYPE_T0 against Knot Atlas targets."""
    global ETA_SIGN, TYPE_T0
    for es in (1, -1):
        for t0 in (0, 1):
            ETA_SIGN, TYPE_T0 = es, t0
            row = ["ETA_SIGN=%+d TYPE_T0=%d:" % (es, t0)]
            try:
                a = goeritz_invariants(parse_pd(TREFOIL))
                b = goeritz_invariants(mirror(parse_pd(TREFOIL)))
                c = goeritz_invariants(parse_pd(FIG8))
                d = goeritz_invariants(parse_pd(K7_4))
                row.append("3_1 s=%+d d=%d | m3_1 s=%+d | 4_1 s=%+d d=%d "
                           "| 7_4 s=%+d d=%d"
                           % (a["sigma"], a["det"], b["sigma"], c["sigma"],
                              c["det"], d["sigma"], d["det"]))
                good = (a["sigma"] == -2 and a["det"] == 3
                        and b["sigma"] == 2 and c["sigma"] == 0
                        and c["det"] == 5 and abs(d["sigma"]) == 2
                        and d["det"] == 15)
                row.append("<== CALIBRATED" if good else "")
            except Exception as ex:
                row.append("EXC: %s" % ex)
            print(" ".join(row))


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: unknot1_decider.py '<PD json>' | --test | --example "
              "| --calibrate")
        return 0
    if argv[0] == "--calibrate":
        calibrate()
        return 0
    if argv[0] == "--test":
        print("sanity suite:")
        f = selftest()
        print("example (11-crossing):")
        r = decide(EXAMPLE11)
        report(r)
        print("suite failures: %d" % f)
        return 1 if f else 0
    if argv[0] == "--example":
        report(decide(EXAMPLE11))
        return 0
    report(decide(json.loads(argv[0])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
