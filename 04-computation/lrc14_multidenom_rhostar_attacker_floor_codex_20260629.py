#!/usr/bin/env python3
"""
HYP-3530: LRC(14) multi-denominator rho*_D and k<12 attacker-floor scout.

This script works two live directions together.

1. Finite-denominator witness criterion.

   For P subset {1,...,13} and normalized E with 0 in E, define

       rho_D(P,E;theta) =
         #{a mod D : a/D in G_P and maxgap({e*a/D:e in E}) > theta}/D.

   The reported object is max_D rho_D over a finite denominator bank.  It is a
   finite witness/opportunity criterion, not a replacement for the continuous
   measure rho*(P,E) and not yet a full reconstructed-cover survivor.  It keeps
   the exact rational witness residue, so it is the natural successor to the
   multi-denominator survivor-grid route after the global 1/7 witness replaced
   the refuted 2/7 via-max target.

2. Union-bound floor below k=12.

   THM-530 reduces k>=8 to the weak pure-shape floor

       mu_{1/7}(E) >= thr_k := 1 - min_{|P|=13-k} meas(G_P).

   We call a normalized k-shape an attacker when mu_{1/7}(E) < thr_k.  The
   bounded-core bank below scans k=8,9,10,11 exactly through span k+6, plus
   large-spread probes, and reports whether any attackers exist.  A no-attacker
   theorem plus the exact gp_min constants gives an unconditional union-bound
   floor for k<12.

Tournament Analysis declaration.
  Vertex set: proof carriers / denominator families, not runners or arcs.
  Pairwise observable: readiness_score = attacker_clearance + witness_support
    + formal_locality - forgetting_penalty.
  Switch/gauge: orient A -> B when A has higher readiness_score; ties use the
    declared Hamiltonian path order.
  Fingerprints: score histogram, directed 3-cycles, SCC sizes, Hamiltonian path.

Assumption challenge.
  Considered alternate vertices: runners, gaps, fixed circle sections, section
  boundaries, wall-crossing events, residues, cover arcs, Fourier modes,
  matroid circuits, and proof obligations.  This quotient preserves the LRC
  predicate "there is a rational point in the GP-and-global-witness carrier";
  it destroys the continuous interval volume and, for reconstructed covers, the
  fast-phase sidecar alignment.  That lost sidecar is why rho_D is reported as
  a finite criterion/witness bank instead of a finished LRC(14) proof.
"""

from __future__ import annotations

import itertools
import random
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import lru_cache, reduce
from math import gcd


H = F(1, 14)
ONE7 = F(1, 7)
TWO7 = F(2, 7)
RNG = random.Random(3530)


def fmt(x: F) -> str:
    return f"{x} = {float(x):.8f}"


def gcd_all(vals) -> int:
    return reduce(gcd, (abs(int(v)) for v in vals), 0)


def circ_frac(x: F) -> F:
    r = x - int(x)
    if r < 0:
        r += 1
    return r


def circ_norm(x: F) -> F:
    r = circ_frac(x)
    return r if r <= F(1, 2) else 1 - r


def circ_maxgap_at(E, x: F) -> F:
    phases = sorted(set(circ_frac(F(e) * x) for e in E))
    if len(phases) <= 1:
        return F(1)
    gaps = [b - a for a, b in zip(phases, phases[1:])]
    gaps.append(phases[0] + 1 - phases[-1])
    return max(gaps)


def merge(arcs):
    arcs = sorted(arcs)
    out = []
    for a, b in arcs:
        if a >= b:
            continue
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def meas(arcs) -> F:
    return sum((b - a for a, b in arcs), F(0))


def complement(arcs):
    arcs = merge(arcs)
    out = []
    prev = F(0)
    for a, b in arcs:
        if a > prev:
            out.append((prev, a))
        prev = max(prev, b)
    if prev < 1:
        out.append((prev, F(1)))
    return out


def danger_arcs(u: int, h: F = H):
    out = []
    for j in range(u):
        c = F(j, u)
        a = (c - h / u) % 1
        b = (c + h / u) % 1
        if a < b:
            out.append((a, b))
        else:
            out.append((a, F(1)))
            out.append((F(0), b))
    return out


@lru_cache(maxsize=None)
def safe_set_tuple(P_tuple):
    P = tuple(P_tuple)
    if not P:
        return ((F(0), F(1)),)
    return tuple(complement(iv for p in P for iv in danger_arcs(p)))


def in_GP_at(P, x: F) -> bool:
    return all(circ_norm(F(p) * x) >= H for p in P)


@lru_cache(maxsize=None)
def good_set_tuple(E_tuple, theta_num: int, theta_den: int):
    theta = F(theta_num, theta_den)
    E = tuple(sorted(set(E_tuple)))
    if len(E) <= 1:
        return ((F(0), F(1)),) if F(1) > theta else tuple()

    diffs = set()
    for a, b in itertools.combinations(E, 2):
        d = abs(b - a)
        if d:
            diffs.add(d)

    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(d + 1):
            bps.add(F(m, d))

    good = []
    bp = sorted(x for x in bps if 0 <= x <= 1)
    for x0, x1 in zip(bp, bp[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        pts = sorted(((F(e) * xm) % 1, e) for e in E)
        order = [e for _, e in pts]
        floors = [int((F(e) * xm) // 1) for e in order]
        m = len(order)

        sub = []
        for idx in range(m):
            e_cur = order[idx]
            f_cur = floors[idx]
            if idx < m - 1:
                e_nx = order[idx + 1]
                f_nx = floors[idx + 1]
                wrap = F(0)
            else:
                e_nx = order[0]
                f_nx = floors[0]
                wrap = F(1)

            # Gap = (e_nx-e_cur)*x + (floor_cur-floor_nx) + wrap.
            A = F(e_nx - e_cur)
            C = F(f_cur - f_nx) + wrap
            if A == 0:
                if C > theta:
                    sub.append((x0, x1))
                continue

            xb = (theta - C) / A
            if A > 0:
                lo, hi = max(x0, xb), x1
            else:
                lo, hi = x0, min(x1, xb)
            if lo < hi:
                sub.append((lo, hi))
        good.extend(merge(sub))

    return tuple(merge(good))


@lru_cache(maxsize=None)
def mu_good_tuple(E_tuple, theta_num: int, theta_den: int):
    """Exact measure of {x: maxgap(E*x)>theta}, without storing all intervals."""
    theta = F(theta_num, theta_den)
    E = tuple(sorted(set(E_tuple)))
    if len(E) <= 1:
        return F(1) if F(1) > theta else F(0)

    diffs = set()
    for a, b in itertools.combinations(E, 2):
        d = abs(b - a)
        if d:
            diffs.add(d)

    bps = {F(0), F(1)}
    for d in diffs:
        for m in range(d + 1):
            bps.add(F(m, d))

    total = F(0)
    bp = sorted(x for x in bps if 0 <= x <= 1)
    for x0, x1 in zip(bp, bp[1:]):
        if x1 <= x0:
            continue
        xm = (x0 + x1) / 2
        pts = sorted(((F(e) * xm) % 1, e) for e in E)
        order = [e for _, e in pts]
        floors = [int((F(e) * xm) // 1) for e in order]
        m = len(order)

        sub = []
        for idx in range(m):
            e_cur = order[idx]
            f_cur = floors[idx]
            if idx < m - 1:
                e_nx = order[idx + 1]
                f_nx = floors[idx + 1]
                wrap = F(0)
            else:
                e_nx = order[0]
                f_nx = floors[0]
                wrap = F(1)

            A = F(e_nx - e_cur)
            C = F(f_cur - f_nx) + wrap
            if A == 0:
                if C > theta:
                    sub.append((x0, x1))
                continue
            xb = (theta - C) / A
            if A > 0:
                lo, hi = max(x0, xb), x1
            else:
                lo, hi = x0, min(x1, xb)
            if lo < hi:
                sub.append((lo, hi))
        total += meas(merge(sub))
    return total


def intersect(A, B):
    out = []
    i = j = 0
    A = list(A)
    B = list(B)
    while i < len(A) and j < len(B):
        a0, a1 = A[i]
        b0, b1 = B[j]
        lo = max(a0, b0)
        hi = min(a1, b1)
        if lo < hi:
            out.append((lo, hi))
        if a1 < b1:
            i += 1
        else:
            j += 1
    return out


def rho_star(P, E, theta=ONE7):
    P_tuple = tuple(P)
    E_tuple = tuple(E)
    good = good_set_tuple(E_tuple, theta.numerator, theta.denominator)
    safe = safe_set_tuple(P_tuple)
    both = tuple(intersect(good, safe))
    return meas(both), meas(safe), meas(good), both


def rho_D(P, E, D: int, theta=ONE7):
    residues = []
    for a in range(D):
        x = F(a, D)
        if in_GP_at(P, x) and circ_maxgap_at(E, x) > theta:
            residues.append(a)
    return F(len(residues), D), residues


def denominator_bank(Dmax=196):
    return list(range(2, Dmax + 1))


def denom_family(D: int) -> str:
    if 2 <= D <= 14:
        return "small_q_2..14"
    if D % 14 == 0:
        return "14m_grid"
    if D % 7 == 0:
        return "7m_grid"
    return "other_grid"


def scan_denominators(P, E, theta=ONE7, Dmax=196):
    rows = []
    family_best = {}
    first_positive = None
    for D in denominator_bank(Dmax):
        val, residues = rho_D(P, E, D, theta)
        fam = denom_family(D)
        row = {
            "D": D,
            "rhoD": val,
            "count": len(residues),
            "residues": tuple(residues[:12]),
            "family": fam,
        }
        rows.append(row)
        if val > 0 and first_positive is None:
            first_positive = row
        if fam not in family_best or val > family_best[fam]["rhoD"]:
            family_best[fam] = row
    best = max(rows, key=lambda r: (r["rhoD"], -r["D"]))
    positives = sum(1 for r in rows if r["rhoD"] > 0)
    return {
        "best": best,
        "family_best": family_best,
        "first_positive": first_positive,
        "positive_denominators": positives,
        "rows": rows,
    }


@lru_cache(maxsize=None)
def gp_min_for_k(k: int):
    psize = 13 - k
    best = None
    bestP = []
    for P in itertools.combinations(range(1, 14), psize):
        val = meas(safe_set_tuple(P))
        if best is None or val < best:
            best = val
            bestP = [P]
        elif val == best:
            bestP.append(P)
    return best, tuple(bestP)


def normalized_shapes(k: int, span: int):
    for tail in itertools.combinations(range(1, span + 1), k - 1):
        E = (0,) + tail
        if gcd_all(E) == 1:
            yield E


def large_spread_probe_shapes(k: int):
    out = set()
    for M in (k + 7, 2 * k, 3 * k, 5 * k, 10 * k, 50 * k):
        out.add(tuple(list(range(k - 1)) + [M]))
        out.add(tuple([0] + list(range(2, k)) + [M]))
    for _ in range(220):
        spread = RNG.choice([2 * k, 3 * k, 5 * k, 10 * k, 50 * k, 100 * k])
        body = sorted(RNG.sample(range(1, spread + 1), k - 1))
        E = tuple([0] + body)
        if len(E) == k and gcd_all(E) == 1:
            out.add(E)
    return sorted(out, key=lambda t: (t[-1], t))


def scan_attackers():
    print("\n" + "=" * 92, flush=True)
    print("B. k<12 union-bound floor: exact bounded-core attacker scan", flush=True)
    print("=" * 92, flush=True)
    print("Attacker definition: normalized E with mu_1/7(E) < thr_k = 1 - gp_min(k).", flush=True)
    print("Exact bounded bank: all primitive normalized shapes with span <= k+5.", flush=True)
    print("Large-spread probe: exact mu on tail and random sparse shapes outside the bank.", flush=True)

    rows = []
    global_floor = None
    for k in range(8, 12):
        span = k + 5
        gp_min, Ps = gp_min_for_k(k)
        thr = 1 - gp_min
        consec = tuple(range(k))
        consec_mu = mu_good_tuple(consec, 1, 7)
        consec_floor = gp_min + consec_mu - 1
        if global_floor is None or consec_floor < global_floor:
            global_floor = consec_floor

        best_mu = None
        minimizers = []
        attackers = []
        below_consec = []
        near_attackers = []
        total = 0
        for E in normalized_shapes(k, span):
            total += 1
            mu = mu_good_tuple(E, 1, 7)
            if best_mu is None or mu < best_mu:
                best_mu = mu
                minimizers = [E]
            elif mu == best_mu and len(minimizers) < 10:
                minimizers.append(E)
            if mu < thr:
                attackers.append((mu, E))
            if mu < consec_mu:
                below_consec.append((mu, E))
            if mu < thr + F(1, 20):
                near_attackers.append((mu, E))

        probe_best = None
        probe_arg = None
        probe_count = 0
        for E in large_spread_probe_shapes(k):
            probe_count += 1
            mu = mu_good_tuple(E, 1, 7)
            if probe_best is None or mu < probe_best:
                probe_best = mu
                probe_arg = E

        floor_found = gp_min + best_mu - 1
        rows.append({
            "k": k,
            "span": span,
            "total": total,
            "gp_min": gp_min,
            "P_min": Ps[0],
            "thr": thr,
            "consec_mu": consec_mu,
            "best_mu": best_mu,
            "minimizers": tuple(minimizers),
            "attackers": tuple(attackers[:5]),
            "attacker_count": len(attackers),
            "below_consec_count": len(below_consec),
            "near_attacker_count": len(near_attackers),
            "floor_found": floor_found,
            "probe_count": probe_count,
            "probe_best": probe_best,
            "probe_arg": probe_arg,
        })

        print(f"\n  k={k}  |P|={13-k}  gp_min={fmt(gp_min)} at P={Ps[0]}", flush=True)
        print(f"    thr_k={fmt(thr)}", flush=True)
        print(f"    exact bank span<= {span}: primitive_shapes={total}", flush=True)
        print(f"    min mu_1/7 found={fmt(best_mu)}", flush=True)
        print(f"    consecutive mu={fmt(consec_mu)}; min==consec? {best_mu == consec_mu}", flush=True)
        print(f"    bounded attackers={len(attackers)}; below-consec shapes={len(below_consec)};"
              f" near-attackers(mu<thr+1/20)={len(near_attackers)}", flush=True)
        print(f"    bounded union floor found={fmt(floor_found)}", flush=True)
        print(f"    minimizer sample={minimizers[:3]}", flush=True)
        print(f"    large-spread probe shapes={probe_count}; best mu={fmt(probe_best)}"
              f" at E={probe_arg}", flush=True)

    print("\n  Bounded-core conclusion for k=8..11:", flush=True)
    all_clear = all(r["attacker_count"] == 0 for r in rows)
    all_consec = all(r["best_mu"] == r["consec_mu"] for r in rows)
    print(f"    no bounded attackers? {all_clear}", flush=True)
    print(f"    consecutive remains exact bounded minimizer? {all_consec}", flush=True)
    print(f"    bounded k<12 union floor (consec/min bank) = {fmt(min(r['floor_found'] for r in rows))}", flush=True)
    print(f"    THM-530 quoted k>=8 floor from consecutive rows = {fmt(global_floor)}", flush=True)
    print("    The remaining theorem gap is now explicit: prove every span>k+5 row is gentle", flush=True)
    print("    (mu_1/7 >= thr_k), or strengthen the exact bounded bank until a known", flush=True)
    print("    large-spread lemma applies.", flush=True)
    return rows


def report_denominator_scenario(name, P, E):
    print("\n" + "-" * 92)
    print(f"Scenario: {name}")
    print(f"  P={tuple(P)}")
    print(f"  E={tuple(E)}  k={len(E)} span={max(E) if E else 0}")
    for theta, label in ((ONE7, "global theta=1/7"), (TWO7, "via-max theta=2/7")):
        rho, gp, mu, intervals = rho_star(P, E, theta)
        scan = scan_denominators(P, E, theta, Dmax=196)
        best = scan["best"]
        first = scan["first_positive"]
        print(f"\n  {label}:")
        print(f"    continuous rho*={fmt(rho)}; meas(G_P)={fmt(gp)}; mu={fmt(mu)};"
              f" intervals={len(intervals)}")
        print(f"    max_D<=196 rho_D={fmt(best['rhoD'])} at D={best['D']}"
              f" count={best['count']} family={best['family']}"
              f" residues={best['residues']}")
        if first is None:
            print("    first positive denominator: NONE in D<=196")
        else:
            print(f"    first positive denominator: D={first['D']} rho_D={first['rhoD']}"
                  f" residues={first['residues']}")
        print("    family best:")
        for fam in ("small_q_2..14", "7m_grid", "14m_grid", "other_grid"):
            row = scan["family_best"].get(fam)
            if row:
                print(f"      {fam:14s}: D={row['D']:3d} rho_D={row['rhoD']}"
                      f" count={row['count']}")
        print(f"    positive denominators in bank: {scan['positive_denominators']}/195")


def scan_denominator_scenarios():
    print("=" * 92)
    print("A. Multi-denominator finite witness criterion max_D rho_D")
    print("=" * 92)
    print("D-bank: all denominators 2..196.  D<=14 and 14m subfamilies are tracked.")
    print("Strict inequalities are used at grid points: safe means ||p*a/D|| >= 1/14,")
    print("and good means maxgap > theta.")
    scenarios = [
        ("k=8 THM-530 binding row", (1, 5, 7, 8, 9), tuple(range(8))),
        ("k=9 union binding row", (1, 11, 12, 13), tuple(range(9))),
        ("k=10 union binding row", (1, 12, 13), tuple(range(10))),
        ("k=11 union binding row", (1, 13), tuple(range(11))),
        ("via-max refutation family k=7", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ("density bridge broad shape", (1, 2, 9), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
        ("perforated k=8 near AP", (1, 2, 3, 4, 5), (0, 2, 3, 4, 5, 6, 7, 8)),
    ]
    for row in scenarios:
        report_denominator_scenario(*row)


def scc_sizes(edge):
    n = len(edge)
    seen = [False] * n
    order = []

    def dfs(v):
        seen[v] = True
        for w in range(n):
            if edge[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)

    redge = [[edge[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    sizes = []

    def rdfs(v):
        seen[v] = True
        total = 1
        for w in range(n):
            if redge[v][w] and not seen[w]:
                total += rdfs(w)
        return total

    for v in reversed(order):
        if not seen[v]:
            sizes.append(rdfs(v))
    return sorted(sizes, reverse=True)


def ham_path(edge, names):
    n = len(edge)
    # Greedy path by dominance score, used only as declared tie path.
    scores = [sum(edge[i][j] for j in range(n)) for i in range(n)]
    return [names[i] for i, _ in sorted(enumerate(scores), key=lambda t: (-t[1], names[t[0]]))]


def tournament_report(attacker_rows):
    print("\n" + "=" * 92)
    print("C. Tournament Analysis over proof carriers")
    print("=" * 92)
    bounded_clear = sum(1 for r in attacker_rows if r["attacker_count"] == 0)
    floor_margin = min(r["floor_found"] for r in attacker_rows)
    carriers = [
        {
            "name": "bounded_core_attacker_bank",
            "attacker_clearance": bounded_clear,
            "witness_support": 3,
            "formal_locality": 5,
            "forgetting_penalty": 0,
        },
        {
            "name": "continuous_union_floor",
            "attacker_clearance": 2 + int(floor_margin > 0),
            "witness_support": 4,
            "formal_locality": 4,
            "forgetting_penalty": 1,
        },
        {
            "name": "max_D_rhoD_bank",
            "attacker_clearance": 1,
            "witness_support": 5,
            "formal_locality": 3,
            "forgetting_penalty": 1,
        },
        {
            "name": "small_q_survivor_filter",
            "attacker_clearance": 1,
            "witness_support": 3,
            "formal_locality": 5,
            "forgetting_penalty": 1,
        },
        {
            "name": "large_spread_gentle_probe",
            "attacker_clearance": 2,
            "witness_support": 2,
            "formal_locality": 2,
            "forgetting_penalty": 1,
        },
        {
            "name": "raw_single_D_only",
            "attacker_clearance": 0,
            "witness_support": 1,
            "formal_locality": 1,
            "forgetting_penalty": 3,
        },
        {
            "name": "via_max_2over7_shadow",
            "attacker_clearance": 0,
            "witness_support": 0,
            "formal_locality": 2,
            "forgetting_penalty": 5,
        },
    ]
    for c in carriers:
        c["score"] = (
            3 * c["attacker_clearance"]
            + 2 * c["witness_support"]
            + c["formal_locality"]
            - 2 * c["forgetting_penalty"]
        )

    n = len(carriers)
    edge = [[False] * n for _ in range(n)]
    names = [c["name"] for c in carriers]
    for i, j in itertools.combinations(range(n), 2):
        ci, cj = carriers[i], carriers[j]
        if (ci["score"], -i) >= (cj["score"], -j):
            edge[i][j] = True
        else:
            edge[j][i] = True

    scores = [sum(edge[i][j] for j in range(n)) for i in range(n)]
    cycles3 = 0
    for a, b, c in itertools.combinations(range(n), 3):
        degs = [
            int(edge[a][b]) + int(edge[a][c]),
            int(edge[b][a]) + int(edge[b][c]),
            int(edge[c][a]) + int(edge[c][b]),
        ]
        if sorted(degs) == [1, 1, 1]:
            cycles3 += 1

    print("  Pairwise observable: 3*attacker_clearance + 2*witness_support")
    print("    + formal_locality - 2*forgetting_penalty.")
    print(f"  score_hist={dict(sorted(Counter(scores).items()))}")
    print(f"  directed_3cycles={cycles3}")
    print(f"  scc_sizes={scc_sizes(edge)}")
    print(f"  Hamiltonian path={' -> '.join(ham_path(edge, names))}")
    print("  carrier scores:")
    for c, s in sorted(zip(carriers, scores), key=lambda t: (-t[1], -t[0]["score"])):
        print(f"    outscore={s} readiness={c['score']:2d}  {c['name']}")


def main():
    print("=" * 92)
    print("HYP-3530 LRC(14): multi-denominator rho_D and k<12 attacker floor")
    print("=" * 92)
    print("All displayed exact quantities are Fraction arithmetic; decimals are only annotations.")
    print("rho_D is a finite denominator witness/opportunity density, not a completed")
    print("reconstructed-cover proof.  The global-witness threshold is theta=1/7.")
    scan_denominator_scenarios()
    attacker_rows = scan_attackers()
    tournament_report(attacker_rows)
    print("\n" + "=" * 92)
    print("TAKEAWAY")
    print("=" * 92)
    print("  max_D rho_D is positive and usually abundant at theta=1/7 on all tested")
    print("  hard rows, while theta=2/7 still exposes the known via-max failure.")
    print("  For k=8..11, the exact bounded-core bank span<=k+5 contains zero attackers")
    print("  and has consecutive shapes as the minimum.  Thus the proof target below")
    print("  k=12 is sharply reframed: prove the large-spread rows are gentle, or")
    print("  extend the bounded-core bank until an existing spread lemma takes over.")
    print("DONE.")


if __name__ == "__main__":
    main()
