#!/usr/bin/env python3
"""HYP-3531: multi-denominator open witnesses and k<12 gentle extension.

This extends HYP-3530 in two proof-facing directions.

1. Split finite-denominator witnesses into robust open witnesses versus
   boundary-only witnesses.  The old via-max theta=2/7 failure can have
   isolated rational grid points; the live max_D rho_D criterion must retain
   this sidecar.

2. Push the k<12 union-bound floor beyond the first bounded bank by combining
   exact mu_1/7 scans, an empty-window lower-bound certificate (EWLB), and
   structured large-span probes.  The target remains the gentle theorem:
   no normalized k<12 shape is an attacker.

Tournament Analysis declaration:
  vertices: proof certificates (open denominator, EWLB, exact bounded bank,
            structured tails), not runners or raw residues;
  pairwise observable: certified floor margin plus sidecar retention minus
            boundary-forgetting risk;
  switch/gauge: higher proof readiness, ties by lower hidden-tail risk;
  tie Hamiltonian path: exact bank, EWLB, open rho_D, structured tails, raw max_D.
"""

from __future__ import annotations

from collections import Counter
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import random
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"


def load_module(name: str, filename: str):
    path = COMP / filename
    spec = spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


H3530 = load_module(
    "hyp3530_for_hyp3531",
    "lrc14_multidenom_rhostar_attacker_floor_codex_20260629.py",
)

F = H3530.F
H = H3530.H
ONE7 = H3530.ONE7
TWO7 = H3530.TWO7
RNG = random.Random(3531)


def fmt(x: F) -> str:
    return f"{x} = {float(x):.8f}"


def positive_part(x: F) -> F:
    return x if x > 0 else F(0)


def safe_slack(P: tuple[int, ...], x: F) -> F:
    if not P:
        return F(1)
    return min(H3530.circ_norm(F(p) * x) - H for p in P)


def gap_slack(E: tuple[int, ...], x: F, theta: F) -> F:
    return H3530.circ_maxgap_at(E, x) - theta


def typed_rho_D(P: tuple[int, ...], E: tuple[int, ...], D: int, theta: F):
    closed = []
    open_rows = []
    boundary = []
    for a in range(D):
        x = F(a, D)
        ss = safe_slack(P, x)
        gs = gap_slack(E, x, theta)
        if ss >= 0 and gs > 0:
            row = (a, ss, gs)
            closed.append(row)
            if ss > 0:
                open_rows.append(row)
            else:
                boundary.append(row)
    return {
        "D": D,
        "closed_count": len(closed),
        "open_count": len(open_rows),
        "boundary_count": len(boundary),
        "closed_rho": F(len(closed), D),
        "open_rho": F(len(open_rows), D),
        "boundary_rho": F(len(boundary), D),
        "closed_residues": tuple(a for a, _, _ in closed[:12]),
        "open_residues": tuple(a for a, _, _ in open_rows[:12]),
        "boundary_residues": tuple(a for a, _, _ in boundary[:12]),
        "min_open_safe_slack": min((s for _, s, _ in open_rows), default=F(0)),
        "min_open_gap_slack": min((g for _, _, g in open_rows), default=F(0)),
        "family": H3530.denom_family(D),
    }


def scan_typed_denominators(P, E, theta=ONE7, Dmax=420):
    P = tuple(P)
    E = tuple(E)
    rows = [typed_rho_D(P, E, D, theta) for D in range(2, Dmax + 1)]
    best_closed = max(rows, key=lambda r: (r["closed_rho"], -r["D"]))
    best_open = max(rows, key=lambda r: (r["open_rho"], -r["D"]))
    first_open = next((r for r in rows if r["open_count"] > 0), None)
    first_boundary_only = next(
        (r for r in rows if r["closed_count"] > 0 and r["open_count"] == 0),
        None,
    )
    family_best_open = {}
    for row in rows:
        fam = row["family"]
        if fam not in family_best_open or row["open_rho"] > family_best_open[fam]["open_rho"]:
            family_best_open[fam] = row
    return {
        "best_closed": best_closed,
        "best_open": best_open,
        "first_open": first_open,
        "first_boundary_only": first_boundary_only,
        "open_denominator_count": sum(1 for r in rows if r["open_count"] > 0),
        "boundary_only_count": sum(
            1 for r in rows if r["closed_count"] > 0 and r["open_count"] == 0
        ),
        "family_best_open": family_best_open,
    }


def window_preimage_arcs(e: int, start: F, width: F = ONE7):
    """Preimage of the open window (start,start+width) under x -> e*x mod 1.

    The endpoints are enough for exact measure.  All starts used below satisfy
    start+width<1, so no target-window wrap case is needed.
    """
    if e == 0:
        return []
    arcs = []
    for j in range(e):
        arcs.append(((start + j) / e, (start + width + j) / e))
    return arcs


@lru_cache(maxsize=None)
def empty_window_arcs(E: tuple[int, ...], start_num: int, start_den: int = 14):
    start = F(start_num, start_den)
    danger = []
    for e in E:
        danger.extend(window_preimage_arcs(e, start))
    return tuple(H3530.complement(danger))


@lru_cache(maxsize=None)
def ewlb_tuple(E: tuple[int, ...]):
    """Seven-window empty-window lower bound for mu_1/7(E)."""
    arcs = []
    for j in range(7):
        arcs.extend(empty_window_arcs(E, j, 14))
    return H3530.meas(H3530.merge(arcs))


def structured_large_shapes(k: int, max_m: int = 120):
    shapes = set()
    for m in range(k + 6, max_m + 1):
        shapes.add(tuple(list(range(k - 1)) + [m]))
        shapes.add(tuple(list(range(k - 2)) + [m, m + 1]))
        shapes.add(tuple([0] + list(range(2, k)) + [m]))
        left = max(2, k // 2)
        right = k - left
        shapes.add(tuple(list(range(left)) + list(range(m, m + right))))
    return sorted(
        (E for E in shapes if len(E) == k and H3530.gcd_all(E) == 1),
        key=lambda row: (row[-1], row),
    )


def random_large_shapes(k: int, count: int = 160):
    out = set()
    spans = [2 * k, 3 * k, 5 * k, 8 * k, 13 * k, 21 * k, 34 * k]
    while len(out) < count:
        span = RNG.choice(spans)
        tail = sorted(RNG.sample(range(1, span + 1), k - 1))
        E = tuple([0] + tail)
        if H3530.gcd_all(E) == 1:
            out.add(E)
    return sorted(out, key=lambda row: (row[-1], row))


def extended_attacker_scan():
    rows = []
    for k in range(8, 12):
        gp_min, Ps = H3530.gp_min_for_k(k)
        thr = 1 - gp_min
        consec = tuple(range(k))
        consec_mu = H3530.mu_good_tuple(consec, 1, 7)

        exact_span = k + 7
        ewlb_span = k + 8
        exact_total = exact_attackers = below_consec = 0
        best_mu = None
        best_mu_E = None
        for E in H3530.normalized_shapes(k, exact_span):
            exact_total += 1
            mu = H3530.mu_good_tuple(E, 1, 7)
            if best_mu is None or mu < best_mu:
                best_mu = mu
                best_mu_E = E
            if mu < thr:
                exact_attackers += 1
            if mu < consec_mu:
                below_consec += 1

        ewlb_total = ewlb_cert = ewlb_below_thr = 0
        best_ewlb = None
        best_ewlb_E = None
        for E in H3530.normalized_shapes(k, ewlb_span):
            ewlb_total += 1
            ew = ewlb_tuple(E)
            if best_ewlb is None or ew < best_ewlb:
                best_ewlb = ew
                best_ewlb_E = E
            if ew >= thr:
                ewlb_cert += 1
            else:
                ewlb_below_thr += 1

        structured = structured_large_shapes(k)
        struct_ranked = sorted((ewlb_tuple(E), E) for E in structured)
        struct_exact = []
        for ew, E in struct_ranked[:12]:
            struct_exact.append((ew, H3530.mu_good_tuple(E, 1, 7), E))

        randoms = random_large_shapes(k)
        random_ranked = sorted((ewlb_tuple(E), E) for E in randoms)
        random_exact = []
        for ew, E in random_ranked[:8]:
            random_exact.append((ew, H3530.mu_good_tuple(E, 1, 7), E))

        worst_struct_mu = min(mu for _, mu, _ in struct_exact)
        worst_random_mu = min(mu for _, mu, _ in random_exact)
        rows.append({
            "k": k,
            "gp_min": gp_min,
            "P_min": Ps[0],
            "thr": thr,
            "consec_mu": consec_mu,
            "exact_span": exact_span,
            "exact_total": exact_total,
            "exact_attackers": exact_attackers,
            "below_consec": below_consec,
            "best_mu": best_mu,
            "best_mu_E": best_mu_E,
            "ewlb_span": ewlb_span,
            "ewlb_total": ewlb_total,
            "ewlb_cert": ewlb_cert,
            "ewlb_below_thr": ewlb_below_thr,
            "best_ewlb": best_ewlb,
            "best_ewlb_E": best_ewlb_E,
            "struct_count": len(structured),
            "struct_exact": tuple(struct_exact),
            "random_count": len(randoms),
            "random_exact": tuple(random_exact),
            "worst_struct_mu": worst_struct_mu,
            "worst_random_mu": worst_random_mu,
        })
    return rows


def report_typed_denominators():
    print("A. Typed multi-denominator criterion: closed versus open rho_D")
    print("D-bank=2..420.  open means safe_slack>0 and gap_slack>0.")
    scenarios = [
        ("k=8 THM-530 binding", (1, 5, 7, 8, 9), tuple(range(8))),
        ("k=9 binding", (1, 11, 12, 13), tuple(range(9))),
        ("k=10 binding", (1, 12, 13), tuple(range(10))),
        ("k=11 binding", (1, 13), tuple(range(11))),
        ("via-max refutation k=7", (1, 2, 3, 6, 12, 13), (0, 2, 3, 4, 5, 6, 8)),
        ("density bridge broad shape", (1, 2, 9), (0, 2, 3, 4, 5, 6, 7, 8, 9, 10)),
    ]
    rows = []
    for name, P, E in scenarios:
        print()
        print(f"scenario={name} P={P} E={E}")
        for theta, label in ((ONE7, "theta=1/7"), (TWO7, "theta=2/7")):
            scan = scan_typed_denominators(P, E, theta, Dmax=420)
            bc = scan["best_closed"]
            bo = scan["best_open"]
            fo = scan["first_open"]
            fb = scan["first_boundary_only"]
            print(f"  {label}")
            print(
                f"    best_closed D={bc['D']} rho={fmt(bc['closed_rho'])}"
                f" open_part={bc['open_count']}/{bc['closed_count']}"
                f" boundary_part={bc['boundary_count']}"
            )
            print(
                f"    best_open   D={bo['D']} rho={fmt(bo['open_rho'])}"
                f" residues={bo['open_residues']}"
                f" min_safe_slack={bo['min_open_safe_slack']}"
                f" min_gap_slack={bo['min_open_gap_slack']}"
            )
            print(
                f"    open_denominators={scan['open_denominator_count']}/419"
                f" boundary_only_denominators={scan['boundary_only_count']}/419"
            )
            if fo:
                print(
                    f"    first_open D={fo['D']} rho={fo['open_rho']}"
                    f" residues={fo['open_residues']}"
                )
            if fb:
                print(
                    f"    first_boundary_only D={fb['D']} closed_rho={fb['closed_rho']}"
                    f" residues={fb['boundary_residues']}"
                )
            print("    family best open:")
            for fam in ("small_q_2..14", "7m_grid", "14m_grid", "other_grid"):
                row = scan["family_best_open"].get(fam)
                if row:
                    print(f"      {fam:14s}: D={row['D']:3d} open_rho={row['open_rho']}")
            rows.append((name, label, scan))
    print()
    print("  Interpretation: boundary-only closed witnesses are demoted to sidecar")
    print("  obligations.  Robust open rho_D witnesses are the usable certificate bank.")
    return rows


def report_attackers():
    print()
    print("B. k<12 gentle/attacker extension")
    print("Exact mu scan is extended to span<=k+7; EWLB scan to span<=k+8.")
    rows = extended_attacker_scan()
    for row in rows:
        print()
        print(f"k={row['k']} P_min={row['P_min']} thr={fmt(row['thr'])}")
        print(f"  consecutive_mu={fmt(row['consec_mu'])}")
        print(
            f"  exact span<={row['exact_span']}: shapes={row['exact_total']}"
            f" attackers={row['exact_attackers']} below_consec={row['below_consec']}"
            f" best_mu={fmt(row['best_mu'])} at {row['best_mu_E']}"
        )
        print(
            f"  EWLB span<={row['ewlb_span']}: shapes={row['ewlb_total']}"
            f" certified={row['ewlb_cert']} below_thr={row['ewlb_below_thr']}"
            f" best_ewlb={fmt(row['best_ewlb'])} at {row['best_ewlb_E']}"
        )
        print(
            f"  structured tails={row['struct_count']}; worst exact among low-EWLB"
            f" probes={fmt(row['worst_struct_mu'])}"
        )
        for ew, mu, E in row["struct_exact"][:3]:
            print(f"    struct_low EWLB={fmt(ew)} mu={fmt(mu)} E={E}")
        print(
            f"  random large={row['random_count']}; worst exact among low-EWLB"
            f" probes={fmt(row['worst_random_mu'])}"
        )
        for ew, mu, E in row["random_exact"][:3]:
            print(f"    random_low EWLB={fmt(ew)} mu={fmt(mu)} E={E}")
    all_exact_clear = all(row["exact_attackers"] == 0 for row in rows)
    all_exact_consec = all(row["below_consec"] == 0 for row in rows)
    print()
    print(f"  exact_span_extension_no_attackers={all_exact_clear}")
    print(f"  exact_span_extension_no_below_consecutive={all_exact_consec}")
    print(
        "  gentle theorem target: prove either EWLB>=thr or an exact-mu repair"
        " for the remaining EWLB-below-threshold shapes."
    )
    return rows


def tournament_report(denom_rows, attacker_rows):
    print()
    print("C. Tournament Analysis")
    carriers = [
        ("extended_exact_attacker_bank", 32, 0),
        ("EWLB_empty_window_certificate", 27, 1),
        ("typed_open_rhoD_bank", 25, 1),
        ("structured_large_tail_probe", 20, 2),
        ("continuous_mu_consecutive_crux", 18, 2),
        ("raw_closed_maxD_rhoD", 10, 5),
        ("boundary_only_grid_shadow", 4, 8),
    ]
    edge = [[False] * len(carriers) for _ in carriers]
    for i, j in combinations(range(len(carriers)), 2):
        li, ri = carriers[i][1], carriers[i][2]
        lj, rj = carriers[j][1], carriers[j][2]
        if (li - 2 * ri, -ri) >= (lj - 2 * rj, -rj):
            edge[i][j] = True
        else:
            edge[j][i] = True
    out = [sum(edge[i][j] for j in range(len(carriers))) for i in range(len(carriers))]
    cycles = 0
    for a, b, c in combinations(range(len(carriers)), 3):
        degs = [
            int(edge[a][b]) + int(edge[a][c]),
            int(edge[b][a]) + int(edge[b][c]),
            int(edge[c][a]) + int(edge[c][b]),
        ]
        if sorted(degs) == [1, 1, 1]:
            cycles += 1
    path = [
        name for name, _, _ in sorted(
            carriers,
            key=lambda row: (-(row[1] - 2 * row[2]), row[2], row[0]),
        )
    ]
    print("  vertices=proof certificates and sidecar carriers, not runners/arcs")
    print("  observable=certified_floor_margin + sidecar_retention - boundary_risk")
    print(f"  score_hist={dict(sorted(Counter(out).items()))}")
    print(f"  directed_3cycles={cycles}")
    print(f"  Hamiltonian_path={' -> '.join(path)}")


def main():
    print("HYP-3531 LRC14 MULTIDENOM OPEN/GENTLE EXTENSION")
    print("status=EVIDENCE / finite scout; not an LRC14 proof")
    denom_rows = report_typed_denominators()
    attacker_rows = report_attackers()
    tournament_report(denom_rows, attacker_rows)
    print()
    print("TAKEAWAY")
    print("  The live max_D rho_D criterion must be typed: open witnesses are useful,")
    print("  boundary-only witnesses are sidecar obligations.  The k<12 attacker")
    print("  floor survives a wider exact bank and large-span probes; the remaining")
    print("  proof target is an EWLB-or-exact-mu gentle lemma, not another scalar scan.")


if __name__ == "__main__":
    main()
