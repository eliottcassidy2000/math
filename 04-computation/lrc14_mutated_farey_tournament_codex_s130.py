#!/usr/bin/env python3
"""S130: mutated Farey operators as LRC14 tournament-analysis carriers.

User prompt: compare four Farey-sequence variations:

  1. numerator multiplied by denominator,
  2. numerator and denominator added,
  3. denominator raised to the numerator,
  4. numerator raised to the denominator.

This script treats each reduced fraction p/q as a two-coordinate address and
tests the four mutations as side-channel "denominators" for LRC14.

The exact LRC14 gap identity is

    M(S)=p/q, threshold=1/14, e=14p-q, gap=M-1/14=e/(14q).

So a proof may safely replace q only with a payload W(p,q) if the resulting
proxy e/W preserves the binding-scale order.  We test this on:

  A. Farey-neighbor locality in the interval [1/14,1/13],
  B. LRC14 AP/GW/petal/single-replacement row families,
  C. tournaments whose vertices are the four mutated payloads plus q.

This is not a proof of LRC14.  It is a filter: the useful mutations are those
that retain enough binding-scale geometry to become proof side channels.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, log
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
THR = F(1, 14)


def load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    assert spec and spec.loader
    module = module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


s124 = load_module("s124_ap_gw", REPO / "04-computation" / "lrc14_ap_gw_common_conditions_codex_s124.py")
s127 = load_module("s127_tournaments", REPO / "04-computation" / "lrc14_tournament_realizability_summit_codex_s127.py")


@dataclass(frozen=True)
class Payload:
    name: str
    desc: str

    def weight_int(self, p: int, q: int) -> int:
        if self.name == "q":
            return q
        if self.name == "sum":
            return p + q
        if self.name == "product":
            return p * q
        if self.name == "denpow":
            return q**p
        if self.name == "numpow":
            return p**q
        raise ValueError(self.name)

    def log_weight(self, p: int, q: int) -> float:
        if self.name == "q":
            return log(q)
        if self.name == "sum":
            return log(p + q)
        if self.name == "product":
            return log(p) + log(q)
        if self.name == "denpow":
            return p * log(q)
        if self.name == "numpow":
            return q * log(p) if p > 0 else float("-inf")
        raise ValueError(self.name)


PAYLOADS = [
    Payload("q", "ordinary Farey denominator"),
    Payload("sum", "p+q additive scan"),
    Payload("product", "p*q multiplicative area scan"),
    Payload("denpow", "q^p denominator-power scan"),
    Payload("numpow", "p^q numerator-power scan"),
]


@dataclass(frozen=True)
class Row:
    label: str
    speeds: tuple[int, ...]
    family: str


def exact_M(row: tuple[int, ...]) -> F:
    M, _pts = s124.M_exact(row)
    return M


def farey_excess(value: F) -> int:
    return 14 * value.numerator - value.denominator


def farey_nodes(lo: F, hi: F, order: int) -> list[F]:
    nodes: list[F] = []
    for q in range(1, order + 1):
        for p in range(0, q + 1):
            if gcd(p, q) == 1:
                x = F(p, q)
                if lo <= x <= hi:
                    nodes.append(x)
    return sorted(set(nodes))


def rank_map(nodes: list[F], payload: Payload) -> dict[F, int]:
    ordered = sorted(
        nodes,
        key=lambda x: (payload.log_weight(x.numerator, x.denominator), x),
    )
    return {x: i for i, x in enumerate(ordered)}


def inversion_count_by_key(items: list, key) -> int:
    inv = 0
    for i in range(len(items)):
        ki = key(items[i])
        for j in range(i + 1, len(items)):
            if ki > key(items[j]):
                inv += 1
    return inv


def farey_locality(order: int = 120) -> dict[str, dict[str, F | int]]:
    nodes = farey_nodes(THR, F(1, 13), order)
    out: dict[str, dict[str, F | int]] = {}
    for payload in PAYLOADS:
        ranks = rank_map(nodes, payload)
        jumps = [abs(ranks[nodes[i + 1]] - ranks[nodes[i]]) for i in range(len(nodes) - 1)]
        inv = inversion_count_by_key(nodes, lambda x: ranks[x])
        out[payload.name] = {
            "nodes": len(nodes),
            "avg_jump": F(sum(jumps), len(jumps)) if jumps else F(0),
            "max_jump": max(jumps) if jumps else 0,
            "inversions": inv,
        }
    return out


def candidate_rows(max_replacement: int = 70) -> list[Row]:
    AP = tuple(range(1, 14))
    rows: list[Row] = [
        Row("AP", AP, "known tight"),
        Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "known tight"),
        Row("residue-liar 12->26", tuple(list(range(1, 12)) + [13, 26]), "q-threshold loose"),
        Row("near-miss 12->36", tuple(list(range(1, 12)) + [13, 36]), "Farey child loose"),
    ]
    for r in range(1, 14):
        for v in range(14, max_replacement + 1):
            row = tuple(sorted((set(AP) - {r}) | {v}))
            if len(row) == 13:
                rows.append(Row(f"swap {r}->{v}", row, "single AP replacement"))
    for m in range(2, 13):
        row = tuple(sorted(set(list(range(1, 12)) + [13, 12 * m])))
        if len(row) == 13:
            rows.append(Row(f"12m family m={m}", row, "12m tail"))

    seen: set[tuple[int, ...]] = set()
    uniq: list[Row] = []
    for row in rows:
        if row.speeds not in seen:
            seen.add(row.speeds)
            uniq.append(row)
    return uniq


def rank_score(true_keys: list[F], pred_keys: list[F]) -> tuple[int, int, int, F]:
    agree = flip = tie = 0
    n = len(true_keys)
    for i in range(n):
        for j in range(i + 1, n):
            td = (true_keys[i] > true_keys[j]) - (true_keys[i] < true_keys[j])
            pd = (pred_keys[i] > pred_keys[j]) - (pred_keys[i] < pred_keys[j])
            if pd == 0:
                tie += 1
            elif pd == td:
                agree += 1
            else:
                flip += 1
    total = agree + flip + tie
    score = F(2 * agree + tie, 2 * total) if total else F(1)
    return agree, flip, tie, score


def proxy_gap(value: F, payload: Payload) -> F:
    e = farey_excess(value)
    if e <= 0:
        return F(0)
    W = payload.weight_int(value.numerator, value.denominator)
    return F(e, W)


def row_bank_analysis(rows: list[Row]) -> tuple[list[F], dict[str, dict[str, object]]]:
    Ms: list[F] = []
    for row in rows:
        Ms.append(exact_M(row.speeds))
    true_keys = [M - THR for M in Ms]
    out: dict[str, dict[str, object]] = {}
    for payload in PAYLOADS:
        pred = [proxy_gap(M, payload) for M in Ms]
        agree, flip, tie, score = rank_score(true_keys, pred)
        worst: list[tuple[F, str, str]] = []
        for i in range(len(rows)):
            for j in range(i + 1, len(rows)):
                td = (true_keys[i] > true_keys[j]) - (true_keys[i] < true_keys[j])
                pd = (pred[i] > pred[j]) - (pred[i] < pred[j])
                if pd and td and pd != td:
                    sep = abs(true_keys[i] - true_keys[j])
                    worst.append((sep, rows[i].label, rows[j].label))
        worst = sorted(worst, reverse=True)[:5]
        out[payload.name] = {
            "agree": agree,
            "flip": flip,
            "tie": tie,
            "score": score,
            "worst": worst,
        }
    return Ms, out


def tournament_from_scores(names: list[str], scores: dict[str, F], larger_better: bool = True) -> int:
    mask = 0
    bit = 0
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            a = scores[names[i]]
            b = scores[names[j]]
            if (a > b) if larger_better else (a < b):
                mask |= 1 << bit
            elif a == b:
                mask |= 1 << bit
            bit += 1
    return mask


def print_child_chain() -> None:
    print("[Unit-excess Farey child chain above 1/14]")
    print("  p/q = p/(14p-1); q follows q -> q+14, the additive n+2-style lane at apex 14.")
    header = "  p  value   q   sum  product  log(q^p)  log(p^q)"
    print(header)
    for p in range(1, 9):
        q = 14 * p - 1
        x = F(p, q)
        print(
            f"  {p:1d}  {str(x):>6s} {q:3d} {p+q:5d} {p*q:8d} "
            f"{p*log(q):9.3f} {q*log(p) if p > 1 else 0.0:9.3f}"
        )
    print()


def print_locality(local: dict[str, dict[str, F | int]]) -> None:
    print("[Farey-order locality in [1/14,1/13], order 120]")
    print("  vertices: reduced Farey fractions in the interval")
    print("  observable: rank distance of adjacent Farey neighbors after mutated ordering")
    print("  switch/gauge: smaller avg_jump / fewer inversions is better")
    for name in [p.name for p in PAYLOADS]:
        data = local[name]
        print(
            f"  {name:8s} nodes={data['nodes']:3d} "
            f"avg_jump={data['avg_jump']} max_jump={data['max_jump']:3d} "
            f"inversions={data['inversions']}"
        )
    print()


def print_row_analysis(rows: list[Row], Ms: list[F], analysis: dict[str, dict[str, object]]) -> None:
    print(f"[LRC14 row-bank proxy analysis: {len(rows)} rows]")
    print("  rows: AP/GW anchors, AP single replacements up to 70, and 12m tails")
    print("  true key: M(S)-1/14 = e/(14q)")
    print("  proxy key: e/W(p,q), with W one mutated Farey payload")
    print("  agreement score = pairwise order agreement against true key")
    for name in [p.name for p in PAYLOADS]:
        d = analysis[name]
        print(
            f"  {name:8s} agree={d['agree']:6d} flip={d['flip']:6d} tie={d['tie']:4d} "
            f"score={d['score']}"
        )
        for sep, a, b in d["worst"][:2]:
            print(f"      large inversion sep={sep}: {a} vs {b}")
    print()

    interesting = [
        "AP",
        "GW 12->24",
        "residue-liar 12->26",
        "near-miss 12->36",
        "swap 10->20",
        "swap 12->24",
        "swap 12->36",
    ]
    by_label = {row.label: (row, M) for row, M in zip(rows, Ms, strict=True)}
    print("[Selected row readout]")
    print(f"  {'row':24s} {'M':>7s} {'e':>3s} {'q':>4s} {'sum':>5s} {'prod':>6s} {'q^p':>12s} {'p^q':>12s}")
    for label in interesting:
        if label in by_label:
            _row, M = by_label[label]
            p, q = M.numerator, M.denominator
            vals = {payload.name: payload.weight_int(p, q) for payload in PAYLOADS}
            print(
                f"  {label:24s} {str(M):>7s} {farey_excess(M):3d} {q:4d} "
                f"{vals['sum']:5d} {vals['product']:6d} "
                f"{vals['denpow']:12d} {vals['numpow']:12d}"
            )
    print()


def print_variant_tournaments(local: dict[str, dict[str, F | int]], analysis: dict[str, dict[str, object]]) -> None:
    names = [p.name for p in PAYLOADS]
    risk_scores = {name: analysis[name]["score"] for name in names}
    local_scores = {name: F(1, 1 + int(local[name]["inversions"])) for name in names}
    risk_mask = tournament_from_scores(names, risk_scores, larger_better=True)
    loc_mask = tournament_from_scores(names, local_scores, larger_better=True)
    majority_scores = {
        name: risk_scores[name] + local_scores[name]
        for name in names
    }
    maj_mask = tournament_from_scores(names, majority_scores, larger_better=True)
    print("[Tournament Analysis on payload variants]")
    print("  vertices: q, sum, product, denpow, numpow")
    print("  risk observable: row-bank pairwise order agreement")
    print("  locality observable: inverse Farey-order inversion count")
    print("  tie Hamiltonian path: q -> sum -> product -> denpow -> numpow")
    for label, mask in [("risk", risk_mask), ("locality", loc_mask), ("majority", maj_mask)]:
        fp = s127.tournament_fingerprint(mask, len(names))
        print(
            f"  {label:8s} fingerprint: score_hist={fp['score_hist']} "
            f"c3={fp['c3']} scc={fp['scc']} hp={fp['hp']}"
        )
    order = sorted(names, key=lambda n: (majority_scores[n], -names.index(n)), reverse=True)
    print(f"  majority Hamiltonian order: {' > '.join(order)}")
    print()


def print_proof_readout() -> None:
    print("[Proof readout]")
    print("  The ordinary denominator q is not just a scale; it is the binding-scale label.")
    print("  The sum payload p+q is the least destructive mutation: on the Farey child chain it")
    print("  is exactly linear, and on the row bank it is close to q but still not exact.")
    print("  The product payload p*q is the additive/multiplicative bridge: it remembers both")
    print("  coordinates and echoes the repo's pair-sum/2n-1 shell arithmetic, but it creates")
    print("  more global inversions than sum.")
    print("  The power payloads are useful as magnitude amplifiers for unbounded false positives,")
    print("  not as local proof denominators: they destroy too much Farey adjacency.")
    print("  LRC14 target sharpened by this scout:")
    print("    use q/binding-scale for the theorem; use sum/product as side-channel ledgers;")
    print("    reserve power variants for detecting unbounded magnitude leakage in fixed-rule")
    print("    tournament discriminators such as floor-odd or CF-parity.")


def main() -> None:
    print("S130 LRC14 MUTATED FAREY TOURNAMENT SCOUT")
    print("=" * 78)
    print("[Assumption challenge]")
    print("  considered vertices: Farey fractions, row families, payload variants, runners,")
    print("    optimum witness scales, and proof obligations.")
    print("  chosen binary relations: mutated-payload order on Farey nodes; proxy-risk order")
    print("    on LRC rows; meta-tournaments on payload variants.")
    print("  preserved predicate: binding-scale closeness to 1/14 via e/(14q).")
    print("  destroyed information: exact safe neighborhoods and wall-state packets.")
    print()
    print_child_chain()
    local = farey_locality()
    print_locality(local)
    rows = candidate_rows()
    Ms, analysis = row_bank_analysis(rows)
    print_row_analysis(rows, Ms, analysis)
    print_variant_tournaments(local, analysis)
    print_proof_readout()


if __name__ == "__main__":
    main()
