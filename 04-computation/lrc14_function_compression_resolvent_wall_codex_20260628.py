#!/usr/bin/env python3
"""HYP-3150 scout: factor-through compression and the k=8 resolvent wall.

The prompt asks to merge Worpitzky pair functions, the K3 edge-flip graph, the
two K4 modeling schemes, and the k=8 biquadratic resolvent into one creative
LRC14 proof push.  This script turns that into a small exact audit:

  compression q: X -> Y
  observable  f: X -> Z
  legal iff f is constant on q-fibers,
  otherwise the missing payload must be named as a sidecar.

Tournament Analysis declaration:
  vertices: proof carriers / functions / fibers / sidecars, not runners;
  pairwise observable: majority vote over retained payload axes;
  switch: orient A -> B when A wins more axes, ties follow a declared path;
  tie path: factor-through audit, resolvent fold, K4 OR compression, K3 kernel,
            PGF curve, sidecars, raw scalar, quintic-wall alarm.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations, product
from math import comb
from typing import Callable, Dict, Hashable, Iterable, List, Sequence, Tuple


Value = Hashable


def fmt_word(word: frozenset[str]) -> str:
    return "E" if not word else "".join(sorted(word))


def factors_through(
    domain: Iterable[Value],
    quotient: Callable[[Value], Value],
    observable: Callable[[Value], Value],
) -> Tuple[bool, Dict[Value, List[Tuple[Value, Value]]]]:
    """Return whether observable is constant on quotient fibers."""
    fibers: Dict[Value, List[Tuple[Value, Value]]] = defaultdict(list)
    for item in domain:
        fibers[quotient(item)].append((item, observable(item)))
    collisions = {
        key: rows
        for key, rows in fibers.items()
        if len({value for _, value in rows}) > 1
    }
    return not collisions, collisions


def print_factor_audit(
    title: str,
    domain: Iterable[Value],
    quotient: Callable[[Value], Value],
    observables: Sequence[Tuple[str, Callable[[Value], Value], str]],
    item_fmt: Callable[[Value], str] = str,
) -> None:
    print(title)
    for name, observable, why in observables:
        ok, collisions = factors_through(domain, quotient, observable)
        verdict = "FACTORS" if ok else "SIDE_CAR_REQUIRED"
        print(f"  {name}: {verdict}; {why}")
        if collisions:
            for key, rows in sorted(collisions.items(), key=lambda kv: str(kv[0]))[:4]:
                rendered = ", ".join(f"{item_fmt(item)}->{value}" for item, value in rows)
                print(f"    fiber {key}: {rendered}")
    print()


# ---------------------------------------------------------------------------
# Pair functions: unordered quotient versus ordered channels.


def pair_function_audit() -> None:
    pairs = [(a, b) for a in range(1, 6) for b in range(1, 6) if a != b]
    quotient = lambda pair: tuple(sorted(pair))
    observables = [
        ("a+b", lambda pair: pair[0] + pair[1], "commutative; unordered edge shadow"),
        ("a*b", lambda pair: pair[0] * pair[1], "commutative; unordered edge shadow"),
        ("a^b", lambda pair: pair[0] ** pair[1], "base/exponent order matters"),
        ("b^a", lambda pair: pair[1] ** pair[0], "base/exponent order matters"),
        (
            "{a^b,b^a}",
            lambda pair: tuple(sorted((pair[0] ** pair[1], pair[1] ** pair[0]))),
            "unordered exponent multiset is safe but weaker than an oriented channel",
        ),
    ]
    print_factor_audit("1. PAIR FUNCTIONS UNDER UNORDERED-PAIR QUOTIENT", pairs, quotient, observables)

    print("pair_function_examples:")
    for a, b in [(2, 3), (2, 4), (3, 5), (4, 5)]:
        print(
            f"  ({a},{b}) vs ({b},{a}): "
            f"sum={a+b}/{b+a} product={a*b}/{b*a} "
            f"a^b={a**b}/{b**a} b^a={b**a}/{a**b}"
        )
    print("readout=sum/product can be quotient currency; exponentials need an ordered endpoint sidecar.")
    print()


# ---------------------------------------------------------------------------
# K3 edge-flip kernel and state-level PGF curves.


K3_EDGES = ((0, 1), (0, 2), (1, 2))


def k3_adj(bits: int) -> List[List[bool]]:
    adj = [[False] * 3 for _ in range(3)]
    for idx, (i, j) in enumerate(K3_EDGES):
        if (bits >> idx) & 1:
            adj[j][i] = True
        else:
            adj[i][j] = True
    return adj


def k3_scores(bits: int) -> Tuple[int, int, int]:
    adj = k3_adj(bits)
    return tuple(sum(1 for j in range(3) if adj[i][j]) for i in range(3))


def k3_class(bits: int) -> str:
    seq = tuple(sorted(k3_scores(bits)))
    if seq == (0, 1, 2):
        return "T"
    if seq == (1, 1, 1):
        return "C"
    raise AssertionError(seq)


def flip_bit(bits: int, idx: int) -> int:
    return bits ^ (1 << idx)


def k3_exit_edges(bits: int) -> Tuple[str, ...]:
    base = k3_class(bits)
    exits = []
    for idx, edge in enumerate(K3_EDGES):
        if k3_class(flip_bit(bits, idx)) != base:
            exits.append(f"{edge[0]}{edge[1]}")
    return tuple(exits)


def k3_transition_signature(bits: int) -> Tuple[Tuple[str, int], ...]:
    counts = Counter(k3_class(flip_bit(bits, idx)) for idx in range(3))
    return tuple(sorted(counts.items()))


def k3_forward_edge_poly(bits: int) -> Tuple[int, int, int]:
    adj = k3_adj(bits)
    poly = [0, 0, 0]
    for path in permutations(range(3)):
        if all(adj[path[i]][path[i + 1]] for i in range(2)):
            ascents = sum(1 for i in range(2) if path[i] < path[i + 1])
            poly[ascents] += 1
    return tuple(poly)


def k3_worpitzky_descents(bits: int) -> Value:
    if k3_class(bits) != "T":
        return "cyclic"
    scores = k3_scores(bits)
    order = tuple(vertex for vertex, _ in sorted(enumerate(scores), key=lambda item: -item[1]))
    return sum(1 for a, b in zip(order, order[1:]) if a > b)


def k3_audit() -> None:
    states = list(range(8))
    print_factor_audit(
        "2. K3 SCORE-CLASS COMPRESSION",
        states,
        k3_class,
        [
            ("class_transition_signature", k3_transition_signature, "the two-state Markov kernel is legal"),
            ("exit_edge_set", k3_exit_edges, "minority/exit edge identity is not legal after class quotient"),
            ("state_forward_edge_PGF", k3_forward_edge_poly, "aggregate PGF can hide state-level curves"),
            ("Worpitzky_descent_word", k3_worpitzky_descents, "descents refine the transitive fiber"),
        ],
        item_fmt=lambda bits: format(bits, "03b"),
    )

    class_counts = Counter(k3_class(bits) for bits in states)
    raw_counts: Dict[str, Counter[str]] = {cls: Counter() for cls in ("C", "T")}
    for bits in states:
        for idx in range(3):
            raw_counts[k3_class(bits)][k3_class(flip_bit(bits, idx))] += 1
    averaged = {
        cls: {dst: Fraction(raw_counts[cls][dst], class_counts[cls]) for dst in ("C", "T")}
        for cls in ("C", "T")
    }

    print("K3_flip_kernel:")
    print(f"  class_counts={dict(sorted(class_counts.items()))}")
    print(f"  raw_counts={{{', '.join(f'{cls}:{dict(raw_counts[cls])}' for cls in ('C','T'))}}}")
    print(f"  averaged_rows_C_T={averaged}")
    print("  normalized_matrix_rows_C_T=[[0,1],[1/3,2/3]], eigenvalues=(1,-1/3)")

    pgf_by_class: Dict[str, Counter[Tuple[int, int, int]]] = defaultdict(Counter)
    aggregate: Dict[str, List[int]] = {cls: [0, 0, 0] for cls in ("C", "T")}
    for bits in states:
        poly = k3_forward_edge_poly(bits)
        cls = k3_class(bits)
        pgf_by_class[cls][poly] += 1
        aggregate[cls] = [aggregate[cls][i] + poly[i] for i in range(3)]
    print("state_level_forward_edge_PGF_curves:")
    for cls in ("C", "T"):
        print(f"  {cls}: curves={dict(sorted(pgf_by_class[cls].items()))} aggregate={tuple(aggregate[cls])}")
    print("readout=K3 class compression preserves the kernel but not the edge role or full PGF curve.")
    print()

    print("Worpitzky_n3_check:")
    for x in range(0, 8):
        rhs = comb(x + 2, 3) + 4 * comb(x + 1, 3) + comb(x, 3)
        print(f"  x={x}: x^3={x**3} rhs={rhs} ok={x**3 == rhs}")
    print()


# ---------------------------------------------------------------------------
# K4 fixed-path cube, OR compression, and canary/fiber sidecars.


K4_CLASS_BY_WORD = {
    frozenset(): "T",
    frozenset({"a"}): "+",
    frozenset({"b"}): "-",
    frozenset({"c"}): "S",
    frozenset({"a", "b"}): "S",
    frozenset({"a", "c"}): "S",
    frozenset({"b", "c"}): "S",
    frozenset({"a", "b", "c"}): "S",
}


def all_k4_words() -> List[frozenset[str]]:
    out = []
    for bits in product((0, 1), repeat=3):
        out.append(frozenset(name for name, bit in zip(("a", "b", "c"), bits) if bit))
    return out


def k4_class(word: frozenset[str]) -> str:
    return K4_CLASS_BY_WORD[word]


def k4_or_compression(word: frozenset[str]) -> Tuple[int, int]:
    return (int("a" in word or "c" in word), int("b" in word or "c" in word))


def k4_xy_class(xy: Tuple[int, int]) -> str:
    return {(0, 0): "T", (1, 0): "+", (0, 1): "-", (1, 1): "S"}[xy]


def k4_delete_one_stable_in_class(word: frozenset[str]) -> bool:
    if not word:
        return False
    base = k4_class(word)
    return all(k4_class(frozenset(set(word) - {letter})) == base for letter in word)


def k4_audit() -> None:
    words = all_k4_words()
    print_factor_audit(
        "3. K4 FIXED-PATH CUBE THROUGH OR COMPRESSION",
        words,
        k4_or_compression,
        [
            ("score_class", k4_class, "class factors through x=a OR c, y=b OR c"),
            ("flip_weight", lambda w: len(w), "fiber PGF/order is destroyed by OR compression"),
            ("flip_word", lambda w: tuple(sorted(w)), "the concrete fixed-path word is destroyed"),
            ("c_canary_status", lambda w: "c" in w, "canary/filler coordinate is not recoverable"),
            (
                "delete_one_stable_status",
                k4_delete_one_stable_in_class,
                "deletion robustness is a sidecar, not a class quotient",
            ),
        ],
        item_fmt=fmt_word,
    )

    fibers: Dict[Tuple[int, int], List[frozenset[str]]] = defaultdict(list)
    for word in words:
        fibers[k4_or_compression(word)].append(word)
    print("K4_OR_fibers:")
    for xy in sorted(fibers):
        members = sorted(fibers[xy], key=lambda w: (len(w), fmt_word(w)))
        pgf = Counter(len(w) for w in members)
        print(
            f"  xy={xy} class={k4_xy_class(xy)} "
            f"members={[fmt_word(w) for w in members]} weight_pgf={dict(sorted(pgf.items()))}"
        )
    class_pgf: Dict[str, Counter[int]] = defaultdict(Counter)
    for word in words:
        class_pgf[k4_class(word)][len(word)] += 1
    print(f"fixed_path_class_weight_PGFs={{{', '.join(f'{cls}:{dict(sorted(class_pgf[cls].items()))}' for cls in ('T','+','-','S'))}}}")
    print("two_bit_scaffold=E->T, x->+, y->-, xy->S is exact; it trades fiber memory for quotient legality.")
    print("readout=the OR map is a legal class compression, but only if fiber PGF and canary/deletion debt are named elsewhere.")
    print()


# ---------------------------------------------------------------------------
# k=8 resolvent degree drop.


def poly_from_roots(roots: Sequence[int]) -> List[int]:
    coeffs = [1]
    for root in roots:
        nxt = [0] * (len(coeffs) + 1)
        for i, coeff in enumerate(coeffs):
            nxt[i] += coeff
            nxt[i + 1] -= coeff * root
        coeffs = nxt
    return coeffs


def eval_poly(coeffs: Sequence[int], x: int) -> int:
    value = 0
    for coeff in coeffs:
        value = value * x + coeff
    return value


def k8_resolvent_audit() -> None:
    roots_t = (1, 2, 4, 5)
    coeff_t = poly_from_roots(roots_t)
    roots_u = tuple(root - 3 for root in roots_t)
    coeff_u = poly_from_roots(roots_u)
    roots_v = tuple(sorted({u * u for u in roots_u}))
    coeff_v = poly_from_roots(roots_v)
    disc_v = coeff_v[1] * coeff_v[1] - 4 * coeff_v[0] * coeff_v[2]

    print("4. K=8 RESOLVENT DEGREE DROP")
    print(f"  g(t) roots={roots_t} coeffs={coeff_t}")
    print(f"  u=t-3 roots={roots_u} coeffs={coeff_u}")
    print("  g(u+3)=u^4-5*u^2+4 has zero odd coefficients")
    print(f"  v=u^2 roots={roots_v} coeffs={coeff_v} discriminant={disc_v}")
    print("  effective_degree: raw_quartic=4, even_variable_v_quadratic=2")
    print("  root_pairs: u=+-1 and u=+-2; the sidecar debt is the lost sign/odd coordinate.")

    domain = list(range(-5, 6))
    print_factor_audit(
        "4b. EVEN RESOLVENT FACTOR-THROUGH TEST",
        domain,
        lambda u: u * u,
        [
            ("g_even(u)", lambda u: u**4 - 5 * u * u + 4, "factors through v=u^2"),
            ("u", lambda u: u, "sign/odd coordinate is destroyed"),
            ("u^3", lambda u: u**3, "odd leakage is destroyed"),
            ("abs(u)", lambda u: abs(u), "a weaker even sidecar is legal"),
        ],
    )

    print("sample_even_values:")
    for u in range(-3, 4):
        print(f"  u={u}: v={u*u} g={u**4 - 5*u*u + 4}")
    print("readout=the current k=8 wall is quartic but not generic quartic; it is a biquadratic with a degree-two visible variable.")
    print()


# ---------------------------------------------------------------------------
# Tournament Analysis over proof carriers.


AXES = (
    "factor_through_exactness",
    "sidecar_hygiene",
    "lrc_transfer",
    "gf_root_curve_retention",
    "degree_control",
    "tournament_kernel_signal",
    "canary_deletion_signal",
    "raw_scalar_warning",
)


@dataclass(frozen=True)
class Carrier:
    name: str
    scores: Dict[str, int]


def carrier(name: str, values: Sequence[int]) -> Carrier:
    return Carrier(name, dict(zip(AXES, values)))


CARRIERS = [
    carrier("factor_through_compression_audit", (10, 10, 10, 7, 8, 7, 7, 10)),
    carrier("k8_even_resolvent_v_u2_fold", (9, 8, 10, 9, 10, 5, 4, 9)),
    carrier("k4_OR_class_compression", (9, 8, 9, 5, 7, 8, 9, 8)),
    carrier("k3_edge_flip_kernel", (8, 7, 8, 5, 6, 10, 4, 7)),
    carrier("fiber_PGF_curve_sidecar", (7, 9, 9, 10, 5, 5, 6, 10)),
    carrier("ordered_pair_exponent_sidecar", (6, 10, 8, 4, 4, 6, 3, 9)),
    carrier("worpitzky_descent_function_sidecar", (7, 9, 8, 7, 5, 7, 3, 8)),
    carrier("c_canary_deletion_sidecar", (7, 10, 9, 4, 5, 7, 10, 8)),
    carrier("symmetric_sum_product_shadow", (8, 4, 5, 2, 4, 3, 2, 5)),
    carrier("raw_scalar_or_class_value", (2, 1, 2, 0, 1, 2, 1, 0)),
    carrier("generic_quintic_wall_alarm", (3, 8, 7, 5, 9, 2, 2, 9)),
]

TIE_PATH = [c.name for c in CARRIERS]


def compare(a: Carrier, b: Carrier, axes: Sequence[str] = AXES) -> str:
    aw = bw = 0
    for axis in axes:
        if a.scores[axis] > b.scores[axis]:
            aw += 1
        elif b.scores[axis] > a.scores[axis]:
            bw += 1
    if aw > bw:
        return a.name
    if bw > aw:
        return b.name
    return a.name if TIE_PATH.index(a.name) < TIE_PATH.index(b.name) else b.name


def tournament_edges(axes: Sequence[str] = AXES) -> Dict[str, set[str]]:
    adj = {c.name: set() for c in CARRIERS}
    for a, b in combinations(CARRIERS, 2):
        winner = compare(a, b, axes)
        loser = b.name if winner == a.name else a.name
        adj[winner].add(loser)
    return adj


def tarjan_scc(adj: Dict[str, set[str]]) -> List[List[str]]:
    index = 0
    stack: List[str] = []
    on_stack: set[str] = set()
    indices: Dict[str, int] = {}
    low: Dict[str, int] = {}
    out: List[List[str]] = []

    def strongconnect(v: str) -> None:
        nonlocal index
        indices[v] = index
        low[v] = index
        index += 1
        stack.append(v)
        on_stack.add(v)
        for w in adj[v]:
            if w not in indices:
                strongconnect(w)
                low[v] = min(low[v], low[w])
            elif w in on_stack:
                low[v] = min(low[v], indices[w])
        if low[v] == indices[v]:
            comp = []
            while True:
                w = stack.pop()
                on_stack.remove(w)
                comp.append(w)
                if w == v:
                    break
            out.append(sorted(comp))

    for v in adj:
        if v not in indices:
            strongconnect(v)
    return sorted(out, key=lambda comp: (-len(comp), comp[0]))


def hamiltonian_path_count(adj: Dict[str, set[str]]) -> int:
    names = list(adj)
    count = 0
    for order in permutations(names):
        if all(order[i + 1] in adj[order[i]] for i in range(len(order) - 1)):
            count += 1
    return count


def directed_3cycle_count(adj: Dict[str, set[str]]) -> int:
    count = 0
    for a, b, c in combinations(adj, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            count += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            count += 1
    return count


def tournament_analysis() -> None:
    adj = tournament_edges()
    score_hist = Counter(len(adj[name]) for name in adj)
    selected_path = sorted(adj, key=lambda name: (-len(adj[name]), TIE_PATH.index(name)))

    raw_axes = ("raw_scalar_warning",)
    raw_adj = tournament_edges(raw_axes)
    edge_flips = 0
    for a, b in combinations(adj, 2):
        proof_winner = a if b in adj[a] else b
        raw_winner = a if b in raw_adj[a] else b
        if proof_winner != raw_winner:
            edge_flips += 1

    print("5. TOURNAMENT ANALYSIS OVER COMPRESSION CARRIERS")
    print("vertices=proof carriers/functions/fibers/sidecars/resolvent variables, not runners")
    print("pairwise_observable=majority over factor-through exactness, sidecar hygiene, LRC transfer, GF/root retention, degree control, tournament kernel signal, canary/deletion signal, raw-scalar warning")
    print("switch=A beats B when A wins more axes; ties follow the declared carrier path")
    print(f"score_hist={dict(sorted(score_hist.items()))}")
    print(f"directed_3cycles={directed_3cycle_count(adj)}")
    print(f"SCCs={tarjan_scc(adj)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(adj)}")
    print(f"edge_flips_against_raw_scalar_warning_gauge={edge_flips}")
    print("selected_path=")
    for name in selected_path:
        print(f"  -> {name} score={len(adj[name])}")
    print()


def hypotheses() -> None:
    print("6. HYPOTHESES GENERATED BY THE AUDIT")
    rows = [
        (
            "H1 factor-through lemma",
            "Every LRC14 quotient step should declare q, f, and whether f is constant on q-fibers before scalarization.",
        ),
        (
            "H2 sidecar taxonomy",
            "The recurring missing data types are ordered endpoint channel, minority edge, PGF curve, canary/deletion status, and odd resolvent coordinate.",
        ),
        (
            "H3 OR compression as filler normal form",
            "The K4 map x=a OR c, y=b OR c is the smallest finite example where a non-linear compression is legal for class but illegal for fiber arithmetic.",
        ),
        (
            "H4 Worpitzky as fiber derivative",
            "Eulerian/Worpitzky data should be treated as a derivative coordinate of the transitive fiber, not as a replacement for the class kernel.",
        ),
        (
            "H5 bounded-degree proof spine",
            "The current LRC14 endgame appears to use degree 2 kernels, two-bit OR compressions, and one quartic that folds to a quadratic; no generic degree-5 object is visible.",
        ),
        (
            "H6 quintic alarm",
            "If a future compression needs five genuinely unordered branches with no even/reflection/filler sidecar, it should be treated as a proof-route alarm rather than a solvable core.",
        ),
        (
            "H7 PGF whole-curve rule",
            "Any equality at a single value, including an aggregate PGF value, must be rechecked at the state/fiber-curve level before it enters HYP-3140.",
        ),
    ]
    for name, text in rows:
        print(f"  {name}: {text}")
    print()


def main() -> None:
    print("HYP-3150 LRC14 function-compression resolvent-degree wall")
    print("source=codex-2026-06-28")
    print()
    pair_function_audit()
    k3_audit()
    k4_audit()
    k8_resolvent_audit()
    tournament_analysis()
    hypotheses()
    print("bottom_line=legal compression is factor-through plus named sidecars; the current hard core stays at effective degree <=4 in this finite audit.")


if __name__ == "__main__":
    main()
