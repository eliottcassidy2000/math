#!/usr/bin/env python3
"""
HYP-3144 / codex-2026-06-27-S274

Three-edge K3 quotient scout for the LRC14 Worpitzky/pair-function prompt.

The script is intentionally small and exact:
- enumerate all labelled tournaments on 3 vertices;
- quotient by score class T=(0,1,2) versus C=(1,1,1);
- flip each edge and compute the quotient kernel;
- compare with the three-coin quotient all-same versus two-to-one mix;
- test which pair functions survive unordered-pair quotienting;
- check the n=3 Worpitzky identity and inspect K3 forward-edge PGFs;
- run a proof-carrier Tournament Analysis.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import comb


EDGES = ((0, 1), (0, 2), (1, 2))
TOURNAMENT_CLASS_ORDER = ("T", "C")
COIN_CLASS_ORDER = ("mix", "same")


def adj_from_bits(bits):
    adj = [[False for _ in range(3)] for _ in range(3)]
    for i, (u, v) in enumerate(EDGES):
        if (bits >> i) & 1:
            adj[u][v] = True
        else:
            adj[v][u] = True
    return adj


def edge_word(bits):
    parts = []
    for i, (u, v) in enumerate(EDGES):
        if (bits >> i) & 1:
            parts.append(f"{u}{v}")
        else:
            parts.append(f"{v}{u}")
    return ",".join(parts)


def out_scores(adj):
    return tuple(sum(1 for w in range(3) if adj[v][w]) for v in range(3))


def score_class(bits):
    scores = tuple(sorted(out_scores(adj_from_bits(bits))))
    if scores == (0, 1, 2):
        return "T"
    if scores == (1, 1, 1):
        return "C"
    raise ValueError(f"unexpected score sequence {scores}")


def flip_bit(bits, edge_index):
    return bits ^ (1 << edge_index)


def topological_order_for_transitive(bits):
    adj = adj_from_bits(bits)
    scores = out_scores(adj)
    return tuple(sorted(range(3), key=lambda v: (-scores[v], v)))


def flip_role(bits, edge_index):
    c0 = score_class(bits)
    c1 = score_class(flip_bit(bits, edge_index))
    if c0 == "C":
        return "cycle_edge_breaks_to_T"
    order = topological_order_for_transitive(bits)
    edge_set = set(EDGES[edge_index])
    if edge_set == {order[0], order[2]}:
        return "long_source_sink_edge_exits_to_C"
    return "adjacent_order_edge_returns_to_T" if c1 == "T" else "unexpected"


def quotient_kernel(states, class_func, class_order):
    class_sizes = Counter(class_func(s) for s in states)
    raw = {a: Counter() for a in class_order}
    for state in states:
        a = class_func(state)
        for i in range(3):
            b = class_func(flip_bit(state, i))
            raw[a][b] += 1

    averaged = {}
    for a in class_order:
        averaged[a] = {}
        for b in class_order:
            averaged[a][b] = Fraction(raw[a][b], class_sizes[a])
    return class_sizes, raw, averaged


def coin_class(bits):
    weight = bits.bit_count()
    return "same" if weight in (0, 3) else "mix"


def normalize_row(row):
    total = sum(row)
    return [Fraction(x, total) for x in row]


def forward_edge_poly(bits):
    adj = adj_from_bits(bits)
    poly = [0, 0, 0]
    paths = []
    for path in permutations(range(3)):
        if all(adj[path[i]][path[i + 1]] for i in range(2)):
            ascents = sum(1 for i in range(2) if path[i] < path[i + 1])
            poly[ascents] += 1
            paths.append((path, ascents))
    return tuple(poly), paths


def eval_poly(coeffs, x):
    return sum(coeffs[i] * (x ** i) for i in range(len(coeffs)))


def solve_fraction_system(mat, rhs):
    n = len(rhs)
    aug = [[Fraction(mat[i][j]) for j in range(n)] + [Fraction(rhs[i])] for i in range(n)]
    for col in range(n):
        pivot = None
        for row in range(col, n):
            if aug[row][col] != 0:
                pivot = row
                break
        if pivot is None:
            raise ValueError("singular system")
        aug[col], aug[pivot] = aug[pivot], aug[col]
        scale = aug[col][col]
        aug[col] = [x / scale for x in aug[col]]
        for row in range(n):
            if row == col:
                continue
            factor = aug[row][col]
            if factor == 0:
                continue
            aug[row] = [aug[row][j] - factor * aug[col][j] for j in range(n + 1)]
    return [aug[i][n] for i in range(n)]


def worpitzky_basis_coeffs_degree2(poly):
    # p(x) = sum_{k=0}^2 w_k * C(x+k, 2)
    points = [0, 1, 2]
    mat = [[comb(x + k, 2) for k in range(3)] for x in points]
    rhs = [eval_poly(poly, x) for x in points]
    return solve_fraction_system(mat, rhs)


def eulerian_number(n, k):
    return sum(((-1) ** j) * comb(n + 1, j) * ((k + 1 - j) ** n) for j in range(k + 2))


def format_fraction(x):
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def format_matrix(mat, rows, cols):
    lines = []
    header = "          " + "  ".join(f"to {c:>5}" for c in cols)
    lines.append(header)
    for r in rows:
        vals = "  ".join(f"{format_fraction(mat[r][c]):>8}" for c in cols)
        lines.append(f"from {r:<4} {vals}")
    return "\n".join(lines)


def pair_function_rows():
    samples = [(2, 3), (2, 4), (4, 2), (3, 5), (5, 3)]
    functions = [
        ("a+b", lambda a, b: a + b, "symmetric_sum", "survives_unordered_pair"),
        ("a*b", lambda a, b: a * b, "symmetric_product", "survives_unordered_pair"),
        ("a^b", lambda a, b: a ** b, "ordered_exponent", "requires_ordered_sidecar"),
        ("b^a", lambda a, b: b ** a, "ordered_exponent_swapped", "requires_ordered_sidecar"),
    ]
    rows = []
    for name, fn, payload, verdict in functions:
        invariant_count = 0
        examples = []
        for a, b in samples:
            v1 = fn(a, b)
            v2 = fn(b, a)
            if v1 == v2:
                invariant_count += 1
            examples.append(f"({a},{b}):{v1}/{v2}")
        rows.append((name, payload, verdict, invariant_count, examples))
    return rows


def build_carrier_tournament():
    weights = {
        "exact_finite": 4,
        "keeps_order": 5,
        "keeps_curve": 4,
        "lrc_transfer": 5,
        "quotient_guard": 4,
        "connects_hyp3143": 3,
        "minimal_toy": 2,
        "scalar_only": -4,
    }
    carriers = [
        ("fiber_pgf_conditional_curve", dict(exact_finite=3, keeps_order=3, keeps_curve=5, lrc_transfer=5, quotient_guard=4, connects_hyp3143=3, minimal_toy=0, scalar_only=0)),
        ("ordered_pair_exponent_sidecar", dict(exact_finite=4, keeps_order=5, keeps_curve=2, lrc_transfer=4, quotient_guard=5, connects_hyp3143=4, minimal_toy=4, scalar_only=0)),
        ("worpitzky_ascent_payload", dict(exact_finite=4, keeps_order=4, keeps_curve=5, lrc_transfer=4, quotient_guard=4, connects_hyp3143=3, minimal_toy=3, scalar_only=0)),
        ("three_edge_flip_kernel", dict(exact_finite=5, keeps_order=3, keeps_curve=2, lrc_transfer=4, quotient_guard=5, connects_hyp3143=5, minimal_toy=5, scalar_only=0)),
        ("packet_subbasis_order_audit", dict(exact_finite=4, keeps_order=4, keeps_curve=1, lrc_transfer=5, quotient_guard=5, connects_hyp3143=5, minimal_toy=2, scalar_only=0)),
        ("tip_tail_commutator_shadow", dict(exact_finite=3, keeps_order=5, keeps_curve=1, lrc_transfer=4, quotient_guard=4, connects_hyp3143=3, minimal_toy=1, scalar_only=0)),
        ("symmetric_sum_product_quotient", dict(exact_finite=5, keeps_order=1, keeps_curve=1, lrc_transfer=3, quotient_guard=2, connects_hyp3143=2, minimal_toy=4, scalar_only=0)),
        ("score_class_scalar", dict(exact_finite=5, keeps_order=0, keeps_curve=0, lrc_transfer=2, quotient_guard=1, connects_hyp3143=2, minimal_toy=5, scalar_only=2)),
        ("raw_single_value", dict(exact_finite=3, keeps_order=0, keeps_curve=0, lrc_transfer=1, quotient_guard=0, connects_hyp3143=0, minimal_toy=3, scalar_only=5)),
    ]
    scores = {}
    for name, features in carriers:
        scores[name] = sum(weights[k] * features.get(k, 0) for k in weights)

    names = [name for name, _ in carriers]
    index = {name: i for i, name in enumerate(names)}
    adj = {name: set() for name in names}
    for a, b in combinations(names, 2):
        if scores[a] > scores[b] or (scores[a] == scores[b] and index[a] < index[b]):
            adj[a].add(b)
        else:
            adj[b].add(a)
    return names, scores, adj


def directed_triangle_count(names, adj):
    count = 0
    for a, b, c in combinations(names, 3):
        if b in adj[a] and c in adj[b] and a in adj[c]:
            count += 1
        if c in adj[a] and b in adj[c] and a in adj[b]:
            count += 1
    return count


def scc_sizes(names, adj):
    index = 0
    stack = []
    on_stack = set()
    indices = {}
    low = {}
    sizes = []

    def strongconnect(v):
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
            size = 0
            while True:
                w = stack.pop()
                on_stack.remove(w)
                size += 1
                if w == v:
                    break
            sizes.append(size)

    for name in names:
        if name not in indices:
            strongconnect(name)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(names, adj):
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    dp = {(1 << i, i): 1 for i in range(n)}
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if val == 0:
                continue
            for nxt_name in adj[names[last]]:
                nxt = idx[nxt_name]
                if mask & (1 << nxt):
                    continue
                dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    return sum(dp.get((full, last), 0) for last in range(n))


def main():
    print("HYP-3144 / codex-2026-06-27-S274")
    print("LRC14 Worpitzky pair-function three-edge quotient scout")
    print()

    states = list(range(8))
    print("1. LABELLED K3 TOURNAMENT STATES")
    class_rows = defaultdict(list)
    f_aggregate = defaultdict(lambda: [0, 0, 0])
    f_distinct = defaultdict(Counter)
    for bits in states:
        cls = score_class(bits)
        adj = adj_from_bits(bits)
        scores = out_scores(adj)
        fpoly, paths = forward_edge_poly(bits)
        class_rows[cls].append(bits)
        f_distinct[cls][fpoly] += 1
        for i, v in enumerate(fpoly):
            f_aggregate[cls][i] += v
        print(f"state={bits:03b} edges={edge_word(bits):8s} scores={scores} class={cls} F={fpoly} H={sum(fpoly)}")
    print()

    print("2. THREE-EDGE SCORE-CLASS FLIP KERNEL")
    sizes, raw, avg = quotient_kernel(states, score_class, TOURNAMENT_CLASS_ORDER)
    print(f"class_sizes={dict(sizes)}")
    print("raw labelled flip counts:")
    print(format_matrix(raw, TOURNAMENT_CLASS_ORDER, TOURNAMENT_CLASS_ORDER))
    print("per-source quotient multiplicities:")
    print(format_matrix(avg, TOURNAMENT_CLASS_ORDER, TOURNAMENT_CLASS_ORDER))
    normalized = {
        "T": {"T": Fraction(2, 3), "C": Fraction(1, 3)},
        "C": {"T": Fraction(1, 1), "C": Fraction(0, 1)},
    }
    print("normalized Markov kernel:")
    print(format_matrix(normalized, TOURNAMENT_CLASS_ORDER, TOURNAMENT_CLASS_ORDER))
    print("unnormalized_eigenvalues=(3,-1)")
    print("normalized_eigenvalues=(1,-1/3)")
    print("stationary_distribution={T:3/4,C:1/4}, matching labelled class split 6:2")
    print()

    print("edge flip roles:")
    role_counts = Counter()
    for bits in states:
        for edge_index in range(3):
            role_counts[flip_role(bits, edge_index)] += 1
    for role, count in sorted(role_counts.items()):
        print(f"  {role}: {count}")
    print("readout=the transitive class has a unique order-sensitive long-edge exit; score class forgets which edge it is")
    print()

    print("3. THREE-COIN QUOTIENT ISOMORPHISM")
    coin_sizes, coin_raw, coin_avg = quotient_kernel(states, coin_class, COIN_CLASS_ORDER)
    print(f"coin_class_sizes={dict(coin_sizes)}")
    print(format_matrix(coin_avg, COIN_CLASS_ORDER, COIN_CLASS_ORDER))
    print("identification: T <-> mix, C <-> same")
    print("readout=three independent flips collapse to the same two-node/three-edge quotient kernel")
    print()

    print("4. FORWARD-EDGE PGF / WORPITZKY PAYLOAD")
    for cls in TOURNAMENT_CLASS_ORDER:
        print(f"class={cls} aggregate_F={tuple(f_aggregate[cls])} distinct_F={dict(f_distinct[cls])}")
        for fpoly in sorted(f_distinct[cls]):
            w = worpitzky_basis_coeffs_degree2(fpoly)
            print(f"  F={fpoly} Worpitzky_basis_coeffs={[format_fraction(x) for x in w]}")
    print("readout=score class and aggregate F can agree while state-level F-curves differ; use the whole PGF curve when order matters")
    print()

    print("n=3 Worpitzky identity check:")
    euler = [eulerian_number(3, k) for k in range(3)]
    print(f"Eulerian_A(3,k)={euler}")
    ok = True
    checks = []
    for x in range(0, 8):
        left = x ** 3
        right = sum(euler[k] * comb(x + k, 3) for k in range(3))
        checks.append(f"x={x}:{left}={right}")
        ok = ok and left == right
    print(" ".join(checks))
    print(f"identity_status={'PASS' if ok else 'FAIL'}")
    print("warning=do not equate the Eulerian/Worpitzky ascent distribution with the score-class count; use it as the ordered-function coordinate system")
    print()

    print("5. PAIR FUNCTION ORDER TEST")
    for name, payload, verdict, invariant_count, examples in pair_function_rows():
        print(f"{name:3s} payload={payload:28s} verdict={verdict:28s} swap_equal_samples={invariant_count}/5 examples={' ; '.join(examples)}")
    print("readout=sum/product are legal unordered-pair scalars; exponentials are oriented-edge data except for accidental equalities")
    print()

    print("6. LRC14 TRANSFER SIGNALS")
    signals = [
        ("three_edge_flip_kernel", "exact K3 quotient matrix [[2,1],[3,0]]"),
        ("antisymmetric_flip_eigenmode", "normalized eigenvalue -1/3"),
        ("ordered_pair_exponent_sidecar", "base/exponent role must survive quotienting"),
        ("score_class_F_curve_collision", "aggregate F can hide state-level F variation"),
        ("fiber_pgf_order_loss_alarm", "HYP-3140 F_R and F_R,Q should not be replaced by a single value before sidecar audit"),
        ("hyp3143_packet_order_preflight", "K3 unique-exit edge is the baby version of K4 lower-order leakage"),
        ("tip_tail_commutator_shadow", "HYP-3141 edge direction is a function argument, not just an arc"),
    ]
    for key, desc in signals:
        print(f"{key}: {desc}")
    print()

    print("7. TOURNAMENT ANALYSIS OVER PROOF CARRIERS")
    names, scores, adj = build_carrier_tournament()
    out_scores_carrier = {name: len(adj[name]) for name in names}
    hist = Counter(out_scores_carrier.values())
    selected_path = sorted(names, key=lambda name: (-scores[name], name))
    print("pairwise_observable=weighted retained proof payload: exactness, order, curve, LRC transfer, quotient guard, HYP-3143 compatibility, minimality")
    print("binary_gauge=A->B iff weighted payload score(A)>score(B); ties use lexical Hamiltonian path")
    print(f"carrier_scores={scores}")
    print(f"score_hist={dict(sorted(hist.items()))}")
    print(f"directed_3cycles={directed_triangle_count(names, adj)}")
    print(f"scc_sizes={scc_sizes(names, adj)}")
    print(f"hamiltonian_path_count={hamiltonian_path_count(names, adj)}")
    print("selected_path=" + " -> ".join(selected_path))
    print()

    print("8. BOLD HYPOTHESES TO TEST NEXT")
    print("H1: the -1/3 antisymmetric K3 mode is a reusable order-loss contraction signature for one-coordinate quotienting in HYP-3140 fiber PGFs.")
    print("H2: HYP-3143's n=4 lower-order leakage is the K3 unique long-edge exit unfolded by one more free bit.")
    print("H3: LRC14 scalar extremality should be phrased as a legal function-evaluation theorem: symmetric functions may quotient; ordered exponent-like functions require a sidecar or named debt.")
    print("H4: Worpitzky ascent payloads are the right coordinate for detecting when a PGF equality at x=1 hides an order-sensitive curve difference.")


if __name__ == "__main__":
    main()
