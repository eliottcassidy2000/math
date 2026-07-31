#!/usr/bin/env python3
"""
Lane C2 (death-star-2026-07-30-coinC2): exhaustive finite-M feasibility of the
doubled-deficit ledger for the depth law d_m = floor(gamma*m) + D0, tested
against the PROVED necessary envelope at a finite set of rational biases.

Frame (THM-2966 spine normal form, laneD ledger conventions):
  row m has 0-side cells (m,k), k=0..d_m, monomial p^z q^o with
  z = m+d_m-k, o = k+1, and 1-side cells with z = k+1, o = m+d_m-k.
  Doubled deficit delta_{cell} is an integer with
      |delta| <= binom(d_m,k)   and   delta == binom(d_m,k)  (mod 2).
  N_M(p) := 2*D_M(p) = sum_{rows m<=M} sum_cells delta * p^z q^o.

PROVED necessary condition (boundary shares, THM-2966 eq (4); laneD Sec 3):
  every fair extractor with deadline T(m) = m+1+d_m satisfies, for all M and
  all p in (0,1):
      |N_M(p)| <= p^{M+1} + q^{M+1}.                          (ENV_M)

This script performs an EXACT forward dynamic program over the cells in row
order (within a row: capacity-descending), state = tuple of exact integer
partial sums N(p_j) * b_j^{A_m} at rational biases p_j = a_j/b_j
(A_m = m + d_m + 1 is the common anti-diagonal level of row m, so every
monomial of row m has exact value integer / b^{A_m}).

Prune rule (exact, integer): a partial sum S is kept iff for every bias
      |S| <= R := (mass of not-yet-assigned cells of the current row)
                  + (p^{m+1} + q^{m+1})            [tail of ALL rows > m]
(all scaled by b^{A_m}; note sum_{m'>m} (p^{m'}q + q^{m'}p) = p^{m+1}+q^{m+1}
exactly, and each row's total cell mass is p^m q + q^m p exactly).

PRUNE-SAFETY LEMMA (proved in the findings file): any assignment of rows
1..M satisfying (ENV_{M'}) for all M' <= M passes every mid-row prune, and
the DP state set at the end of row M equals EXACTLY the set of evaluation
vectors of assignments of rows 1..M satisfying (ENV_1), ..., (ENV_M).
Corollary: if the state set becomes EMPTY at row M, then NO assignment of
rows 1..M satisfies all boundary envelopes at the chosen biases; by (ENV)
necessity this is a THEOREM: no fair extractor with this depth law exists.
Minimal such M is exact (states at row M-1 were nonempty).

Symmetrization: if the bias multiset is closed under p -> 1-p, the boundary
state set is invariant under the side-swap mirror (N-vector permuted by the
conjugation involution), because boundary sets are exactly envelope-
characterized and the envelope is mirror-symmetric.  We then keep one
representative per mirror orbit at each row boundary (exact 2x saving).

Everything asserted is exact integer arithmetic; floats appear only in logs.

Usage:
  python3 amm12592_laneC2_finiteM_bnb_deathstar.py --selftest
  python3 amm12592_laneC2_finiteM_bnb_deathstar.py --gamma 1/2 --D0 0 \
      --biases 1/2,2/3,3/4 --Mmax 24 --time-budget 330 --state-cap 1200000
"""

import argparse
import sys
import time
from math import comb


def parse_frac(s):
    a, b = s.split("/")
    a, b = int(a), int(b)
    from math import gcd
    g = gcd(a, b)
    a //= g
    b //= g
    assert 0 < a < b, (a, b)
    return (a, b)


def depth(m, gnum, gden, D0):
    return (gnum * m) // gden + D0


def conj_perm(biases):
    """Permutation sigma with biases[sigma[j]] = conjugate of biases[j], or None."""
    sigma = []
    for (a, b) in biases:
        c = (b - a, b)
        if c not in biases:
            return None
        sigma.append(biases.index(c))
    return sigma


def run(gnum, gden, D0, biases, Mmax, tbudget, scap, log=print):
    r = len(biases)
    sigma = conj_perm(biases)
    log(f"# depth law d_m = floor({gnum}/{gden} m) + {D0};  biases = "
        + ",".join(f"{a}/{b}" for a, b in biases)
        + f";  Mmax={Mmax} tbudget={tbudget}s scap={scap} mirror_sym={'ON' if sigma else 'off'}")
    states = {(0,) * r}
    t0 = time.time()
    A_prev = None
    summary = {"verdict": "RUNNING", "rows": []}
    for m in range(1, Mmax + 1):
        d = depth(m, gnum, gden, D0)
        A = m + d + 1
        if A_prev is not None:
            sh = A - A_prev
            mults = tuple(b ** sh for (a, b) in biases)
            states = {tuple(N * mults[j] for j, N in enumerate(st)) for st in states}
        A_prev = A
        # cells of row m, both sides
        cells = []
        for side in (0, 1):
            for k in range(d + 1):
                cap = comb(d, k)
                if side == 0:
                    z, o = m + d - k, k + 1
                else:
                    z, o = k + 1, m + d - k
                mono = tuple(a ** z * (b - a) ** o for (a, b) in biases)
                cells.append((cap, side, k, z, o, mono))
        cells.sort(key=lambda c: (-c[0], c[1], c[2]))
        # tail mass of ALL rows > m, scaled by b^A: (a^{m+1}+c^{m+1}) * b^{d}
        tail = tuple((a ** (m + 1) + (b - a) ** (m + 1)) * b ** d for (a, b) in biases)
        # suffix masses: suf[i] = tail + sum over cells i..end of cap*mono
        suf = [tail]
        acc = list(tail)
        for c in reversed(cells):
            cap, mono = c[0], c[5]
            acc = [acc[j] + cap * mono[j] for j in range(r)]
            suf.append(tuple(acc))
        suf.reverse()
        # process cells
        for i, (cap, side, k, z, o, mono) in enumerate(cells):
            rem = suf[i + 1]
            new = set()
            add = new.add
            rng = range(r)
            for st in states:
                lo, hi = -cap, cap
                ok = True
                for j in rng:
                    mj = mono[j]
                    sj = st[j]
                    rj = rem[j]
                    h = (rj - sj) // mj
                    if h < hi:
                        hi = h
                    l = -((rj + sj) // mj)
                    if l > lo:
                        lo = l
                    if lo > hi:
                        ok = False
                        break
                if not ok:
                    continue
                if (lo - cap) & 1:
                    lo += 1
                for dv in range(lo, hi + 1, 2):
                    add(tuple(st[j] + dv * mono[j] for j in rng))
            states = new
            if not states:
                log(f"INFEASIBLE: state set EMPTY at row m={m}, after cell #{i} "
                    f"(cap={cap} side={side} k={k} z={z} o={o}).")
                log(f"THEOREM-GRADE: no assignment of rows 1..{m} satisfies "
                    f"ENV_1..ENV_{m} at the listed biases; depth law refuted.")
                summary["verdict"] = "INFEASIBLE"
                summary["min_infeasible_M"] = m
                summary["cell"] = (cap, side, k, z, o)
                return summary
            if len(states) > scap:
                log(f"STOP: state cap exceeded ({len(states)} > {scap}) at row m={m} cell #{i}; "
                    f"exhaustive verdict only through row {m-1}.")
                summary["verdict"] = f"FEASIBLE_THROUGH_{m-1}_STATECAP"
                return summary
            if time.time() - t0 > tbudget:
                log(f"STOP: time budget exceeded at row m={m} cell #{i}; "
                    f"exhaustive verdict only through row {m-1}.")
                summary["verdict"] = f"FEASIBLE_THROUGH_{m-1}_TIME"
                return summary
        # mirror-orbit canonicalization at row boundary (exact; see docstring)
        if sigma is not None:
            states = {min(st, tuple(st[sigma[j]] for j in range(r))) for st in states}
        # row-boundary stats (floats for logs only)
        env = tail  # at row boundary the prune bound IS the envelope, scaled by b^A
        stats = []
        bestrho = None
        for st in states:
            rho = max(abs(st[j]) / env[j] for j in range(r))
            if bestrho is None or rho < bestrho:
                bestrho = rho
        for j, (a, b) in enumerate(biases):
            mn = min(st[j] for st in states)
            mx = max(st[j] for st in states)
            mna = min(abs(st[j]) for st in states)
            den = float(b) ** A
            stats.append((mn / den, mx / den, mna / den))
        n = len(states)
        line = (f"ROW m={m} d={d} A={A} states={n} t={time.time()-t0:.1f}s bestrho={bestrho:.4f} | "
                + " | ".join(f"p={a}/{b}: N in [{s[0]:.3e},{s[1]:.3e}] min|N|={s[2]:.3e}"
                             for (a, b), s in zip(biases, stats)))
        log(line)
        summary["rows"].append(
            dict(m=m, d=d, A=A, states=n, bestrho=bestrho,
                 stats=[(f"{a}/{b}",) + s for (a, b), s in zip(biases, stats)]))
        # exact small-int corridor for bias 1/2 if present (checkpointable exactly)
        for j, (a, b) in enumerate(biases):
            if (a, b) == (1, 2):
                mn = min(st[j] for st in states)
                mx = max(st[j] for st in states)
                mna = min(abs(st[j]) for st in states)
                env_half = env[j]
                log(f"  exact@1/2: A={A} envN={env_half} minN={mn} maxN={mx} min|N|={mna}")
    summary["verdict"] = f"FEASIBLE_THROUGH_{Mmax}_MMAX"
    return summary


# ---------------------------------------------------------------- selftest ---

def brute_force_boundary_set(gnum, gden, D0, biases, M):
    """All evaluation vectors of assignments of rows 1..M satisfying every
    boundary envelope ENV_1..ENV_M.  Exponential; only for tiny cases."""
    r = len(biases)
    frontier = {(0,) * r}
    A_prev = None
    for m in range(1, M + 1):
        d = depth(m, gnum, gden, D0)
        A = m + d + 1
        if A_prev is not None:
            sh = A - A_prev
            mults = tuple(b ** sh for (a, b) in biases)
            frontier = {tuple(N * mults[j] for j, N in enumerate(st)) for st in frontier}
        A_prev = A
        cells = []
        for side in (0, 1):
            for k in range(d + 1):
                cap = comb(d, k)
                z, o = (m + d - k, k + 1) if side == 0 else (k + 1, m + d - k)
                mono = tuple(a ** z * (b - a) ** o for (a, b) in biases)
                cells.append((cap, mono))
        env = tuple((a ** (m + 1) + (b - a) ** (m + 1)) * b ** d for (a, b) in biases)
        # enumerate ALL row assignments (no pruning), then filter by ENV_m
        rowsums = [(0,) * r]
        for cap, mono in cells:
            nxt = []
            for st in rowsums:
                for dv in range(-cap, cap + 1, 2):
                    nxt.append(tuple(st[j] + dv * mono[j] for j in range(r)))
            rowsums = nxt
        newf = set()
        for st in frontier:
            for rs in rowsums:
                cand = tuple(st[j] + rs[j] for j in range(r))
                if all(abs(cand[j]) <= env[j] for j in range(r)):
                    newf.add(cand)
        frontier = newf
    return frontier


def selftest():
    ok = True
    # Test 1: gamma=1/2 D0=0, biases {1/2, 2/3}, M=3: DP boundary set == brute force
    biases = [(1, 2), (2, 3)]
    bf = brute_force_boundary_set(1, 2, 0, biases, 3)
    # rerun DP capturing the boundary set at M=3
    states = {(0,) * len(biases)}
    log_sink = lambda *a, **k: None
    # inline DP replicating run() but returning boundary states at M=3
    dp = _dp_boundary_states(1, 2, 0, biases, 3)
    print(f"selftest1 gamma=1/2: brute={len(bf)} dp={len(dp)} equal={bf == dp}")
    ok &= (bf == dp)
    # Test 2: gamma=1 D0=0, biases {1/2, 3/4}, M=2
    biases = [(1, 2), (3, 4)]
    bf = brute_force_boundary_set(1, 1, 0, biases, 2)
    dp = _dp_boundary_states(1, 1, 0, biases, 2)
    print(f"selftest2 gamma=1:   brute={len(bf)} dp={len(dp)} equal={bf == dp}")
    ok &= (bf == dp)
    # Test 3: mirror symmetrization consistency: symmetric bias set, orbit expansion
    biases = [(1, 2), (1, 3), (2, 3)]
    bf = brute_force_boundary_set(1, 2, 0, biases, 3)
    dp = _dp_boundary_states(1, 2, 0, biases, 3)
    sigma = conj_perm(biases)
    canon_bf = {min(st, tuple(st[sigma[j]] for j in range(len(biases)))) for st in bf}
    canon_dp = {min(st, tuple(st[sigma[j]] for j in range(len(biases)))) for st in dp}
    print(f"selftest3 mirror:    brute-orbits={len(canon_bf)} dp-orbits={len(canon_dp)} "
          f"equal={canon_bf == canon_dp} raw_equal={bf == dp}")
    ok &= (canon_bf == canon_dp and bf == dp)
    print("SELFTEST", "PASS" if ok else "FAIL")
    return 0 if ok else 1


def _dp_boundary_states(gnum, gden, D0, biases, M):
    """DP identical to run() (no symmetrization, no caps), returning the
    boundary state set after row M."""
    r = len(biases)
    states = {(0,) * r}
    A_prev = None
    for m in range(1, M + 1):
        d = depth(m, gnum, gden, D0)
        A = m + d + 1
        if A_prev is not None:
            sh = A - A_prev
            mults = tuple(b ** sh for (a, b) in biases)
            states = {tuple(N * mults[j] for j, N in enumerate(st)) for st in states}
        A_prev = A
        cells = []
        for side in (0, 1):
            for k in range(d + 1):
                cap = comb(d, k)
                z, o = (m + d - k, k + 1) if side == 0 else (k + 1, m + d - k)
                mono = tuple(a ** z * (b - a) ** o for (a, b) in biases)
                cells.append((cap, side, k, z, o, mono))
        cells.sort(key=lambda c: (-c[0], c[1], c[2]))
        tail = tuple((a ** (m + 1) + (b - a) ** (m + 1)) * b ** d for (a, b) in biases)
        suf = [tail]
        acc = list(tail)
        for c in reversed(cells):
            cap, mono = c[0], c[5]
            acc = [acc[j] + cap * mono[j] for j in range(r)]
            suf.append(tuple(acc))
        suf.reverse()
        for i, (cap, side, k, z, o, mono) in enumerate(cells):
            rem = suf[i + 1]
            new = set()
            for st in states:
                lo, hi = -cap, cap
                dead = False
                for j in range(r):
                    mj, sj, rj = mono[j], st[j], rem[j]
                    h = (rj - sj) // mj
                    if h < hi:
                        hi = h
                    l = -((rj + sj) // mj)
                    if l > lo:
                        lo = l
                    if lo > hi:
                        dead = True
                        break
                if dead:
                    continue
                if (lo - cap) & 1:
                    lo += 1
                for dv in range(lo, hi + 1, 2):
                    new.add(tuple(st[j] + dv * mono[j] for j in range(r)))
            states = new
    return states


# =========================================================================
# Lane G3 additions (death-star coinC2, resumed session).
#
# The full-frontier DP above is EXACT but its state set explodes at row 6
# (~1.6e6 states, runB/runC logs).  The additions below answer targeted
# questions with first-witness depth-first search + sound exact pruning
# (style of the main session's amm12592_laneC2_finiteM_dfs_deathstar.py,
# logic re-derived here so this file is self-contained):
#
#   g3_dfs      first-witness feasibility of ENV_1..ENV_M at rational biases,
#               optionally with per-cell CLAMPS (forced values) and/or an
#               EQ-TARGET "final scaled N(p_j0) == v exactly".
#   mode=dfs        feasibility sweep M=1.. for a depth law (par ON/OFF)
#   mode=corridor   exact achievable SET of final N_M(1/2) values (corridor)
#   mode=pin        per-cell feasible value sets (clamp tests)
#   mode=witness    extract a surviving witness, exact greedy comparison,
#                   homogenized band content + forced-parity referee assert
#
# Soundness of the prune (same lemma as the DP): a prefix with partial sum
# S dies only if for some checkpoint m' >= current row
#   |S| > (wiggle of unassigned cells of current row)
#         + sum_{rows r in (row,m']} (full row mass) + Env_{m'},
# which every ENV-satisfying completion obeys.  For the eq-target we prune
# when |v - S[j0]| > (all remaining cell mass at j0 through row M), also
# sound.  Hence "INFEASIBLE" verdicts are exhaustive-search theorems about
# the boundary-envelope relaxation; witnesses are certified by re-check.
# =========================================================================


class G3Budget(Exception):
    pass


def g3_build(g1, g2, D0, Mtry, biases):
    """Row data at common scale den^Lmax, Lmax = A(Mtry)."""
    dep = lambda m: (g1 * m) // g2 + D0
    Afun = lambda m: m + dep(m) + 1
    Lmax = Afun(Mtry)
    nP = len(biases)
    rows, ids, Arow = [], [], []
    for m in range(1, Mtry + 1):
        dm = dep(m)
        cells, cid = [], []
        for k in range(dm + 1):
            cap = comb(dm, k)
            for side, (z, o) in enumerate(((m + dm - k, k + 1),
                                           (k + 1, m + dm - k))):
                w = tuple(n ** z * (b - n) ** o * b ** (Lmax - z - o)
                          for (n, b) in biases)
                cells.append((cap, w))
                cid.append((m, side, k, z, o))
        rows.append(cells)
        ids.append(cid)
        Arow.append(Afun(m))
    R = [tuple(sum(c[0] * c[1][i] for c in rows[m]) for i in range(nP))
         for m in range(Mtry)]
    Env = [tuple((n ** (m + 2) + (b - n) ** (m + 2)) * b ** (Lmax - m - 2)
                 for (n, b) in biases) for m in range(Mtry)]
    return rows, ids, Arow, R, Env, Lmax


def g3_dfs(g1, g2, D0, Mtry, biases, parity=True, clamps=None,
           eq_target=None, node_cap=None, time_budget=None,
           value_order="abs", zero_bias=None):
    """First-witness DFS.  Returns (status, witness, nodes) with status in
    {FEASIBLE, INFEASIBLE, CAP}; witness = {(m,side,k): value} on success.
    eq_target = (bias_index j0, exact scaled integer v): additionally demand
    final N(p_j0)*b^Lmax == v.  clamps = {(m,side,k): forced value}.
    value_order: 'abs' (small |v| first) or 'desc' (saturation-first, +cap
    down to -cap) — ordering only affects speed, never the verdict.
    zero_bias: index jz — replace ENV at bias jz by the far stronger demand
    N_m(p_jz) == 0 at EVERY row end (Env[m][jz] := 0; prune stays sound
    since the modified constraint set is still row-monotone)."""
    rows, ids, Arow, R, Env, Lmax = g3_build(g1, g2, D0, Mtry, biases)
    nP = len(biases)
    clamps = clamps or {}
    if zero_bias is not None:
        Env = [tuple(0 if i == zero_bias else e[i] for i in range(nP))
               for e in Env]
    allowed = []
    for m in range(Mtry):
        best = None
        for mp in range(m, Mtry):
            tot = tuple(Env[mp][i] +
                        sum(R[r][i] for r in range(m + 1, mp + 1))
                        for i in range(nP))
            best = tot if best is None else tuple(
                min(a, b) for a, b in zip(best, tot))
        allowed.append(best)
    suffix = []
    for m in range(Mtry):
        sw = [tuple(0 for _ in range(nP))]
        for c in reversed(rows[m]):
            sw.append(tuple(sw[-1][i] + c[0] * c[1][i] for i in range(nP)))
        sw.reverse()
        suffix.append(sw)
    later = [tuple(sum(R[r][i] for r in range(m + 1, Mtry))
                   for i in range(nP)) for m in range(Mtry)]
    t0 = time.time()
    state = {"nodes": 0}
    assign = {}

    def choices(cap, key, v_hint_desc):
        par = cap % 2 if parity else -1
        if key in clamps:
            v = clamps[key]
            assert abs(v) <= cap and (par == -1 or (v - cap) % 2 == 0), \
                ("bad clamp", key, v, cap)
            return [v]
        if par == -1:
            vals = list(range(-cap, cap + 1))
        else:
            vals = [v for v in range(-cap, cap + 1) if (v - par) % 2 == 0]
        if v_hint_desc:
            return sorted(vals, reverse=True)
        return sorted(vals, key=abs)

    hint_desc = bool(eq_target and eq_target[1] > 0) or value_order == "desc"

    def dfs(m, ci, S):
        state["nodes"] += 1
        if node_cap and state["nodes"] > node_cap:
            raise G3Budget
        if time_budget and state["nodes"] % 4096 == 0 \
                and time.time() - t0 > time_budget:
            raise G3Budget
        cells = rows[m]
        if ci == len(cells):
            for i in range(nP):
                if abs(S[i]) > Env[m][i] or abs(S[i]) > allowed[m][i]:
                    return False
            if m + 1 == Mtry:
                if eq_target is not None and S[eq_target[0]] != eq_target[1]:
                    return False
                return True
            return dfs(m + 1, 0, S)
        rest = suffix[m][ci]
        for i in range(nP):
            if abs(S[i]) > rest[i] + allowed[m][i]:
                return False
        if eq_target is not None:
            j0, v = eq_target
            if abs(v - S[j0]) > rest[j0] + later[m][j0]:
                return False
        cap, w = cells[ci]
        key = ids[m][ci][:3]
        for v in choices(cap, key, hint_desc):
            assign[key] = v
            S2 = tuple(S[i] + v * w[i] for i in range(nP)) if v else S
            if dfs(m, ci + 1, S2):
                return True
        del assign[key]
        return False

    try:
        ok = dfs(0, 0, (0,) * nP)
    except G3Budget:
        return "CAP", None, state["nodes"]
    if not ok:
        return "INFEASIBLE", None, state["nodes"]
    # certify the witness by independent exact re-evaluation
    wit = dict(assign)
    S = [0] * nP
    for m in range(Mtry):
        for (cap, w), cid in zip(rows[m], ids[m]):
            v = wit.get(cid[:3], 0)
            assert abs(v) <= cap and (not parity or (v - cap) % 2 == 0)
            for i in range(nP):
                S[i] += v * w[i]
        for i in range(nP):
            assert abs(S[i]) <= Env[m][i], ("witness fails ENV", m, i)
        if eq_target is not None and m + 1 == Mtry:
            assert S[eq_target[0]] == eq_target[1]
    return "FEASIBLE", wit, state["nodes"]


def g3_beta_parity(M, A_M, d_M):
    """Coefficient parities of beta_M(x) = (1+x)^A_M + (1+x^{M+1})(1+x)^d_M."""
    beta = [0] * (A_M + 1)
    for o in range(A_M + 1):
        b = comb(A_M, o) & 1
        if o <= d_M:
            b ^= comb(d_M, o) & 1
        if 0 <= o - (M + 1) <= d_M:
            b ^= comb(d_M, o - (M + 1)) & 1
        beta[o] = b
    return beta


def g3_homogenize(g1, g2, D0, Mtry, wit):
    """Exact homogenized doubled-ledger coefficients c_O (O = q-exponent) at
    level A_M for the assignment wit; also returns (A_M, d_M)."""
    dep = lambda m: (g1 * m) // g2 + D0
    A_M = Mtry + dep(Mtry) + 1
    coeff = [0] * (A_M + 1)
    for m in range(1, Mtry + 1):
        dm = dep(m)
        t = A_M - (m + dm + 1)
        for k in range(dm + 1):
            for side, (z, o) in enumerate(((m + dm - k, k + 1),
                                           (k + 1, m + dm - k))):
                v = wit.get((m, side, k), 0)
                if v:
                    for j in range(t + 1):
                        coeff[o + j] += v * comb(t, j)
    return coeff, A_M, dep(Mtry)


def g3_greedy(g1, g2, D0, Mtry, biases, jref=0):
    """In-lane greedy: cells in the same order as g3_dfs, each cell takes the
    parity-allowed value minimizing |S'[jref]| (ties -> smaller |v|), no
    lookahead, envelopes NOT enforced.  Returns (wit, rowends, firstviol),
    firstviol = smallest 1-based row where some ENV_m fails (None if none)."""
    rows, ids, Arow, R, Env, Lmax = g3_build(g1, g2, D0, Mtry, biases)
    nP = len(biases)
    S = (0,) * nP
    wit = {}
    rowends = []
    firstviol = None
    for m in range(Mtry):
        for (cap, w), cid in zip(rows[m], ids[m]):
            par = cap % 2
            vals = [v for v in range(-cap, cap + 1) if (v - par) % 2 == 0]
            best = min(vals, key=lambda v: (abs(S[jref] + v * w[jref]), abs(v)))
            wit[cid[:3]] = best
            S = tuple(S[i] + best * w[i] for i in range(nP))
        rowends.append(S)
        if firstviol is None and any(abs(S[i]) > Env[m][i] for i in range(nP)):
            firstviol = m + 1
    return wit, rowends, firstviol, Env, Lmax


def g3_mode_extend(args, biases):
    """Incremental frontier growth: witness at M0 by full DFS, then extend
    row by row, clamping already-solved rows; on failure free the last
    `back` rows (back = 1, 2, ... up to full) before giving up.  Every
    FEASIBLE step is a certified witness at that M (re-checked in g3_dfs);
    an INFEASIBLE verdict here is only conclusive when back == M (full)."""
    gnum, gden = parse_frac(args.gamma)
    M0, M1 = args.Mlist[0], args.Mlist[-1]
    st, wit, nodes = g3_dfs(gnum, gden, args.D0, M0, biases, parity=True,
                            node_cap=args.node_cap,
                            time_budget=args.time_budget,
                            value_order=args.value_order)
    print(f"[g3 extend g={args.gamma} D0={args.D0}] base M={M0}: {st} "
          f"nodes={nodes}", flush=True)
    if st != "FEASIBLE":
        return
    for M in range(M0 + 1, M1 + 1):
        done = False
        for back in list(range(1, min(6, M) + 1)) + [M]:
            clamps = {key: v for key, v in wit.items() if key[0] <= M - back}
            t0 = time.time()
            st, w2, nodes = g3_dfs(gnum, gden, args.D0, M, biases,
                                   parity=True, clamps=clamps,
                                   node_cap=args.node_cap,
                                   time_budget=args.time_budget,
                                   value_order=args.value_order)
            print(f"[g3 extend] M={M} back={back} clamped_rows<={M-back}: "
                  f"{st} nodes={nodes} t={time.time()-t0:.1f}s", flush=True)
            if st == "FEASIBLE":
                wit = w2
                done = True
                break
            if st == "CAP":
                continue   # try freeing more rows before giving up
        if not done:
            print(f"[g3 extend] frontier stops at M={M-1} "
                  f"(no certified verdict at M={M} within budget)",
                  flush=True)
            return
        # summarize witness saturation
        nz = sat = tot = 0
        for m in range(1, M + 1):
            dm = (gnum * m) // gden + args.D0
            for k in range(dm + 1):
                for side in (0, 1):
                    v = wit.get((m, side, k), 0)
                    c = comb(dm, k)
                    tot += 1
                    if v:
                        nz += 1
                    if abs(v) == c:
                        sat += 1
        print(f"[g3 extend] M={M} witness: cells={tot} nonzero={nz} "
              f"saturated={sat}", flush=True)


def g3_mode_dfs(args, biases):
    gnum, gden = parse_frac(args.gamma)
    for Mtry in range(1, args.Mmax + 1):
        t0 = time.time()
        jz = None
        if args.zero_bias:
            jz = biases.index(parse_frac(args.zero_bias))
        st, wit, nodes = g3_dfs(gnum, gden, args.D0, Mtry, biases,
                                parity=not args.no_parity,
                                node_cap=args.node_cap,
                                time_budget=args.time_budget,
                                value_order=args.value_order,
                                zero_bias=jz)
        print(f"[g3 dfs g={args.gamma} D0={args.D0} "
              f"par={'OFF' if args.no_parity else 'ON'} "
              f"ord={args.value_order}"
              + (f" zero@{args.zero_bias}" if args.zero_bias else "")
              + f"] M={Mtry:2d}: "
              f"{st:10s} nodes={nodes} t={time.time()-t0:.1f}s", flush=True)
        if st == "INFEASIBLE":
            print(f"*** exhausted search: ENV_1..ENV_{Mtry} at these biases "
                  f"UNSATISFIABLE -> depth-law THEOREM ***", flush=True)
            break
        if st == "CAP":
            break


def g3_mode_corridor(args, biases):
    gnum, gden = parse_frac(args.gamma)
    j0 = biases.index((1, 2))
    for Mtry in args.Mlist:
        rows, ids, Arow, R, Env, Lmax = g3_build(gnum, gden, args.D0,
                                                 Mtry, biases)
        envN = Env[Mtry - 1][j0]
        # forced parity of the final scaled sum at j0 (parity ON only):
        # S = sum v*w with v == cap (mod 2), so S == sum_{w odd} cap (mod 2).
        P0 = None
        if not args.no_parity:
            P0 = sum(c[0] for m in range(Mtry) for c in rows[m]
                     if c[1][j0] % 2 == 1) % 2
        ach, inf, cap = [], [], []
        vrange = (range(envN, -1, -1) if args.scan_desc
                  else range(0, envN + 1))
        for v in vrange:
            if P0 is not None and (v - P0) % 2 != 0:
                inf.append(v)      # proved by parity, no search needed
                continue
            st, wit, nodes = g3_dfs(gnum, gden, args.D0, Mtry, biases,
                                    parity=not args.no_parity,
                                    eq_target=(j0, v),
                                    node_cap=args.node_cap,
                                    time_budget=args.time_budget)
            print(f"    M={Mtry} v={v}: {st} nodes={nodes}", flush=True)
            (ach if st == "FEASIBLE" else inf if st == "INFEASIBLE"
             else cap).append(v)
            if args.scan_desc and st == "FEASIBLE":
                break   # corridor max found; lower |v| values not scanned
        print(f"[g3 corridor g={args.gamma} D0={args.D0} "
              f"par={'OFF' if args.no_parity else 'ON'}] M={Mtry} "
              f"Lmax={Lmax} envN={envN} | achievable(v>=0)={ach} | "
              f"infeasible={inf} | undecided={cap}", flush=True)
        print(f"    (by global-negation symmetry the corridor is "
              f"-max..max; max_achievable={max(ach) if ach else None}, "
              f"width/env={ (max(ach)/envN) if ach else 0 :.3f})", flush=True)


def g3_mode_pin(args, biases):
    gnum, gden = parse_frac(args.gamma)
    Mtry = args.Mlist[0]
    rows, ids, Arow, R, Env, Lmax = g3_build(gnum, gden, args.D0,
                                             Mtry, biases)
    rowsel = set(args.pin_rows) if args.pin_rows else None
    print(f"[g3 pin g={args.gamma} D0={args.D0} M={Mtry}] per-cell feasible "
          f"|value| sets under ENV_1..ENV_{Mtry} at {len(biases)} biases "
          f"(global +- symmetry => sets are sign-symmetric; testing v>=0)"
          + (f" rows={sorted(rowsel)}" if rowsel else ""), flush=True)
    for m in range(Mtry):
        if rowsel is not None and (m + 1) not in rowsel:
            continue
        for (capc, w), cid in zip(rows[m], ids[m]):
            key = cid[:3]
            par = capc % 2
            feas, undec = [], []
            for v in range(par, capc + 1, 2):
                st, wit, nodes = g3_dfs(gnum, gden, args.D0, Mtry, biases,
                                        parity=True, clamps={key: v},
                                        node_cap=args.node_cap,
                                        time_budget=args.time_budget)
                (feas if st == "FEASIBLE" else undec
                 if st == "CAP" else []).append(v)
            tagbits = []
            if len(feas) == 1 and not undec:
                tagbits.append("PINNED|v|=%d" % feas[0])
            mm, side, k = key
            print(f"  cell m={mm} side={side} k={k} cap={capc} "
                  f"(z={cid[3]},o={cid[4]}): feasible|v|={feas}"
                  + (f" undecided={undec}" if undec else "")
                  + ("  <== " + ",".join(tagbits) if tagbits else ""),
                  flush=True)


def g3_mode_witness(args, biases):
    gnum, gden = parse_frac(args.gamma)
    Mtry = args.Mlist[0]
    st, wit, nodes = g3_dfs(gnum, gden, args.D0, Mtry, biases, parity=True,
                            node_cap=args.node_cap,
                            time_budget=args.time_budget)
    print(f"[g3 witness g={args.gamma} D0={args.D0} M={Mtry}] "
          f"status={st} nodes={nodes}", flush=True)
    if st != "FEASIBLE":
        return
    rows, ids, Arow, R, Env, Lmax = g3_build(gnum, gden, args.D0,
                                             Mtry, biases)
    nP = len(biases)
    # exact witness trajectory at row ends
    S = (0,) * nP
    print("  witness cells (nonzero):", flush=True)
    traj = []
    for m in range(Mtry):
        for (cap, w), cid in zip(rows[m], ids[m]):
            v = wit.get(cid[:3], 0)
            if v:
                print(f"    m={cid[0]} side={cid[1]} k={cid[2]} "
                      f"(z={cid[3]},o={cid[4]}) cap={cap} delta={v}",
                      flush=True)
            S = tuple(S[i] + v * w[i] for i in range(nP))
        traj.append(S)
        rho = max(abs(S[i]) / Env[m][i] for i in range(nP))
        j0 = biases.index((1, 2)) if (1, 2) in biases else 0
        print(f"  rowend m={m+1}: N(1/2)*2^Lmax={S[j0]} "
              f"env={Env[m][j0]} maxrho={rho:.4f}", flush=True)
    # greedy comparison
    gw, grow, gviol, gEnv, gL = g3_greedy(gnum, gden, args.D0, Mtry, biases)
    j0 = biases.index((1, 2))
    print(f"  greedy(min |N(1/2)| per cell, no lookahead): first ENV "
          f"violation at M={gviol}", flush=True)
    for m in range(Mtry):
        print(f"    greedy rowend m={m+1}: N(1/2)scaled={grow[m][j0]} "
              f"env={gEnv[m][j0]}"
              + ("  VIOL" if any(abs(grow[m][i]) > gEnv[m][i]
                                 for i in range(nP)) else ""), flush=True)
    diffs = [k for k in sorted(set(wit) | set(gw))
             if wit.get(k, 0) != gw.get(k, 0)]
    print(f"  witness-vs-greedy differing cells ({len(diffs)}): "
          f"{diffs[:40]}", flush=True)
    # homogenized band content + forced-parity referee (PROVED beta_M)
    for name, W in (("witness", wit), ("greedy", gw)):
        coeff, A_M, d_M = g3_homogenize(gnum, gden, args.D0, Mtry, W)
        beta = g3_beta_parity(Mtry, A_M, d_M)
        assert all((coeff[o] & 1) == beta[o] for o in range(A_M + 1)), \
            f"forced-parity referee FAILED for {name}"
        band = list(range(d_M + 2, Mtry))
        print(f"  {name}: homogenized level A_M={A_M}; forced-parity "
              f"referee PASS; band o in [{d_M+2},{Mtry-1}] coeffs "
              f"{[coeff[o] for o in band]}", flush=True)
    print("  (band parities are choice-independent =="
          " beta_M restriction — PROVED cross-check)", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gamma", default="1/2")
    ap.add_argument("--D0", type=int, default=0)
    ap.add_argument("--biases", default="1/2")
    ap.add_argument("--Mmax", type=int, default=24)
    ap.add_argument("--time-budget", type=float, default=330.0)
    ap.add_argument("--state-cap", type=int, default=1_200_000)
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--mode", default="dp",
                    choices=["dp", "dfs", "corridor", "pin", "witness",
                             "extend"])
    ap.add_argument("--value-order", default="abs", choices=["abs", "desc"])
    ap.add_argument("--scan-desc", action="store_true",
                    help="corridor mode: scan v from envN down, stop at "
                         "first feasible (= corridor max)")
    ap.add_argument("--zero-bias", default="",
                    help="dfs mode: bias (e.g. 1/2) forced to N_m == 0 at "
                         "every row end (replaces its envelope)")
    ap.add_argument("--no-parity", action="store_true")
    ap.add_argument("--node-cap", type=int, default=None)
    ap.add_argument("--Mlist", default="10",
                    help="comma list of M values for corridor/pin/witness")
    ap.add_argument("--pin-rows", default="",
                    help="comma list of 1-based rows to test in pin mode")
    args = ap.parse_args()
    args.pin_rows = [int(s) for s in args.pin_rows.split(",") if s]
    if args.selftest:
        sys.exit(selftest())
    args.Mlist = [int(s) for s in args.Mlist.split(",")]
    biases = [parse_frac(s) for s in args.biases.split(",")]
    if args.mode == "dp":
        gnum, gden = parse_frac(args.gamma)
        summary = run(gnum, gden, args.D0, biases, args.Mmax,
                      args.time_budget, args.state_cap)
        print("VERDICT:", summary["verdict"])
    elif args.mode == "dfs":
        g3_mode_dfs(args, biases)
    elif args.mode == "corridor":
        g3_mode_corridor(args, biases)
    elif args.mode == "pin":
        g3_mode_pin(args, biases)
    elif args.mode == "witness":
        g3_mode_witness(args, biases)
    elif args.mode == "extend":
        g3_mode_extend(args, biases)


if __name__ == "__main__":
    main()
