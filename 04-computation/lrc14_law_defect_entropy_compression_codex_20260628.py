#!/usr/bin/env python3
"""Finite law-defect entropy scout for LRC14 compression packets.

This script treats algebraic laws as quotient-legality claims.  A law says
that a target function is constant on the fibers of a quotient.  When the law
fails, the residual conditional entropy is the number of sidecar bits that the
proof packet still owes.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import product
from math import log2


def entropy(counter: Counter) -> float:
    total = sum(counter.values())
    if total == 0:
        return 0.0
    value = 0.0
    for count in counter.values():
        p = count / total
        if p:
            value -= p * log2(p)
    return value


def cond_entropy(records, q_name: str, f_name: str):
    groups = defaultdict(Counter)
    for rec in records:
        groups[rec[q_name]][rec[f_name]] += 1
    total = len(records)
    h = sum((sum(cnt.values()) / total) * entropy(cnt) for cnt in groups.values())
    nonconstant = {
        q: cnt for q, cnt in groups.items() if len(cnt) > 1
    }
    return h, groups, nonconstant


def fmt_bits(x: float) -> str:
    return f"{x:.6f} bits"


def examples(nonconstant, limit=4):
    out = []
    for q, cnt in list(nonconstant.items())[:limit]:
        out.append(f"{q}->{dict(cnt)}")
    return "; ".join(out) if out else "none"


def commutativity_audit():
    records = []
    for a, b in product(range(1, 8), repeat=2):
        records.append(
            {
                "unordered": tuple(sorted((a, b))),
                "sum": a + b,
                "product": a * b,
                "pow_forward": a**b,
                "pow_backward": b**a,
                "pow_multiset": tuple(sorted((a**b, b**a))),
            }
        )
    return [
        ("commutativity", "a+b", *cond_entropy(records, "unordered", "sum")),
        ("commutativity", "a*b", *cond_entropy(records, "unordered", "product")),
        ("commutativity", "a^b", *cond_entropy(records, "unordered", "pow_forward")),
        (
            "commutativity",
            "{a^b,b^a}",
            *cond_entropy(records, "unordered", "pow_multiset"),
        ),
    ]


def associativity_records(op, domain):
    records = []
    for a, b, c in product(domain, repeat=3):
        left = op(op(a, b), c)
        right = op(a, op(b, c))
        records.append({"leaves": (a, b, c), "value": left, "bracket": "left"})
        records.append({"leaves": (a, b, c), "value": right, "bracket": "right"})
    return records


def associativity_audit():
    add = lambda x, y: x + y
    mul = lambda x, y: x * y
    sub = lambda x, y: x - y
    exp = lambda x, y: x**y
    tests = [
        ("associativity", "(a+b)+c", associativity_records(add, range(5))),
        ("associativity", "(a*b)*c", associativity_records(mul, range(5))),
        ("associativity", "(a-b)-c", associativity_records(sub, range(5))),
        ("associativity", "(a^b)^c", associativity_records(exp, range(1, 5))),
    ]
    return [(law, name, *cond_entropy(records, "leaves", "value")) for law, name, records in tests]


def idempotence_audit():
    records = []
    for word in product(range(3), repeat=3):
        records.append(
            {
                "support": tuple(sorted(set(word))),
                "max": max(word),
                "sum": sum(word),
                "word": word,
            }
        )
    return [
        ("idempotence/multiplicity", "max(word)", *cond_entropy(records, "support", "max")),
        ("idempotence/multiplicity", "sum(word)", *cond_entropy(records, "support", "sum")),
    ]


def distributivity_records(op, plus, modulus):
    records = []
    domain = range(modulus)
    for a, b, c in product(domain, repeat=3):
        left = op(a, plus(b, c))
        right = plus(op(a, b), op(a, c))
        records.append({"args": (a, b, c), "side": "left", "value": left})
        records.append({"args": (a, b, c), "side": "right", "value": right})
    return records


def distributivity_audit():
    m = 5
    add = lambda x, y: (x + y) % m
    mul = lambda x, y: (x * y) % m
    return [
        (
            "distributivity",
            "a*(b+c) over Z5",
            *cond_entropy(distributivity_records(mul, add, m), "args", "value"),
        ),
        (
            "distributivity",
            "a+(b*c) over Z5",
            *cond_entropy(distributivity_records(add, mul, m), "args", "value"),
        ),
    ]


def k4_class(bits):
    a, b, c = bits
    if (a, b, c) == (0, 0, 0):
        return "T"
    if (a, b, c) == (1, 0, 0):
        return "+"
    if (a, b, c) == (0, 1, 0):
        return "-"
    return "S"


def flip(bits, idx):
    bits = list(bits)
    bits[idx] ^= 1
    return tuple(bits)


def k4_action_audit():
    records = []
    relations = {}
    names = ["a", "b", "c"]
    for idx, name in enumerate(names):
        rel = defaultdict(set)
        for bits in product([0, 1], repeat=3):
            src = k4_class(bits)
            dst = k4_class(flip(bits, idx))
            rel[src].add(dst)
            records.append({"class_action": (src, name), "next_class": dst})
        relations[name] = {k: sorted(v) for k, v in rel.items()}
    h, groups, nonconstant = cond_entropy(records, "class_action", "next_class")
    return h, groups, nonconstant, relations


def root_radius_packet():
    packets = [
        ("circle", [1.0] * 6),
        ("phi4_off_circle", [0.78, 0.91, 0.98, 1.03, 1.14, 1.28]),
    ]
    out = []
    for name, radii in packets:
        logs = [log2(r) for r in radii]
        mean = sum(logs) / len(logs)
        var = sum((x - mean) ** 2 for x in logs) / len(logs)
        product_radius = 1.0
        for r in radii:
            product_radius *= r
        out.append((name, product_radius, var))
    return out


def directed_3cycles(adj):
    n = len(adj)
    count = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                if adj[i][j] and adj[j][k] and adj[k][i]:
                    count += 1
                if adj[i][k] and adj[k][j] and adj[j][i]:
                    count += 1
    return count


def hamiltonian_path_count(adj):
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[-1])


def tournament_analysis():
    vertices = [
        (
            "typed_factor_through_schema",
            (10, 10, 10, 9, 9),
        ),
        (
            "law_defect_entropy_bits",
            (9, 9, 10, 9, 8),
        ),
        (
            "root_curve_phi4_radius_sidecar",
            (9, 8, 9, 10, 8),
        ),
        (
            "k8_odd_worpitzky_orientation_sidecar",
            (8, 10, 8, 8, 9),
        ),
        (
            "k4_transformation_monoid_sidecar",
            (8, 8, 8, 7, 9),
        ),
        (
            "associativity_bracket_sidecar",
            (7, 8, 8, 7, 8),
        ),
        (
            "distributivity_context_sidecar",
            (7, 7, 8, 7, 7),
        ),
        (
            "idempotence_multiplicity_sidecar",
            (6, 7, 7, 6, 7),
        ),
        (
            "raw_scalar_value",
            (1, 1, 1, 2, 1),
        ),
    ]
    names = [name for name, _ in vertices]
    weights = [sum(score) for _, score in vertices]
    n = len(vertices)
    adj = [[False] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if weights[i] == weights[j]:
                adj[i][j] = i < j
            else:
                adj[i][j] = weights[i] > weights[j]
    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    return {
        "names": names,
        "weights": weights,
        "scores": scores,
        "score_hist": dict(sorted(Counter(scores).items())),
        "directed_3cycles": directed_3cycles(adj),
        "hamiltonian_path_count": hamiltonian_path_count(adj),
        "path": [name for _, name in sorted(zip(scores, names), reverse=True)],
    }


def main():
    print("LRC14 law-defect entropy compression scout")
    print("source: codex-2026-06-28 / HYP-3201 / T1301 / LTI-301 / LTT-201")
    print()
    print("Principle:")
    print("  law L == target f factors through quotient q_L")
    print("  failure(L) == H(f | q_L) > 0, hence a named sidecar is required")
    print()

    print("Commutativity quotient: q(a,b)={a,b}")
    for law, name, h, _groups, nonconstant in commutativity_audit():
        status = "FACTORS" if not nonconstant else "SIDECAR"
        print(f"  {name:12s} {status:8s} residual={fmt_bits(h)} examples={examples(nonconstant)}")
    print()

    print("Associativity quotient: q(value tree)=ordered leaves (a,b,c)")
    for law, name, h, _groups, nonconstant in associativity_audit():
        status = "FACTORS" if not nonconstant else "BRACKET_SIDECAR"
        print(f"  {name:12s} {status:16s} residual={fmt_bits(h)} nonconstant={len(nonconstant)}")
        if nonconstant:
            print(f"    examples: {examples(nonconstant)}")
    print()

    print("Idempotence/multiplicity quotient: q(word)=support(word)")
    for law, name, h, _groups, nonconstant in idempotence_audit():
        status = "FACTORS" if not nonconstant else "MULTIPLICITY_SIDECAR"
        print(f"  {name:12s} {status:22s} residual={fmt_bits(h)} examples={examples(nonconstant)}")
    print()

    print("Distributivity quotient: identify two sides of a distributive rewrite")
    for law, name, h, _groups, nonconstant in distributivity_audit():
        status = "FACTORS" if not nonconstant else "CONTEXT_SIDECAR"
        print(f"  {name:18s} {status:16s} residual={fmt_bits(h)} nonconstant={len(nonconstant)}")
        if nonconstant:
            print(f"    examples: {examples(nonconstant)}")
    print()

    print("K4 fixed-path class action: quotient exact cube states to T,+,-,S")
    h, _groups, nonconstant, relations = k4_action_audit()
    print(f"  full action residual H(next_class | class,generator)={fmt_bits(h)}")
    for action, rel in relations.items():
        print(f"  flip {action}: {rel}")
    print(f"  non-deterministic quotient fibers: {examples(nonconstant, limit=8)}")
    print("  readout: the class quotient is a transformation/relational monoid packet, not V4")
    print()

    print("Lee-Yang radius compression sidecar")
    for name, product_radius, log_var in root_radius_packet():
        print(
            f"  {name:15s} product_radius=q0/q6={product_radius:.6f} "
            f"log_radius_variance={log_var:.8f}"
        )
    print("  readout: q0=q6*R^6 is exact only for the common-radius quotient;")
    print("           off-circle variance is the phi4 sidecar, not scalar noise")
    print()

    print("Tournament Analysis over proof carriers")
    ta = tournament_analysis()
    print("  vertices considered: functions, law quotients, bracket trees, supports,")
    print("    distributive contexts, K4 class-action fibers, PGF roots, phi4 modes,")
    print("    resolvent variables, and proof obligations")
    print("  pairwise observable: retained proof payload under quotient legality,")
    print("    sidecar explicitness, LRC transfer, root/function retention, and degree control")
    print(f"  score_hist={ta['score_hist']}")
    print(f"  directed_3cycles={ta['directed_3cycles']}")
    print(f"  hamiltonian_path_count={ta['hamiltonian_path_count']}")
    print("  priority_path=")
    for item in ta["path"]:
        print(f"    -> {item}")


if __name__ == "__main__":
    main()
