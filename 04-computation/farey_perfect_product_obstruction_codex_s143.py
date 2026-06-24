#!/usr/bin/env python3
"""Farey product/perfect-number/Kuratowski guardrail scout.

This is a synthesis script, not an LRC14 proof.  It tests the prompt's
observation that a reduced Farey term a/b can carry the product payload
ab = |E(K_{a,b})|, then compares the ordinary n=2 unit-excess chain that
contains even perfect numbers with the LRC14 n=14 unit-excess chain.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from fractions import Fraction
from itertools import combinations
from math import gcd, isqrt


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def factorization(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d = 3 if d == 2 else d + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def fmt_factor(n: int) -> str:
    parts = []
    for p, e in sorted(factorization(n).items()):
        parts.append(str(p) if e == 1 else f"{p}^{e}")
    return "*".join(parts) if parts else "1"


def sigma(n: int) -> int:
    ans = 1
    for p, e in factorization(n).items():
        ans *= (p ** (e + 1) - 1) // (p - 1)
    return ans


def fmt_fraction(fr: Fraction) -> str:
    if fr.denominator == 1:
        return str(fr.numerator)
    return f"{fr.numerator}/{fr.denominator}"


def farey_terms(order: int) -> list[Fraction]:
    terms = {
        Fraction(a, b)
        for b in range(1, order + 1)
        for a in range(0, b + 1)
        if gcd(a, b) == 1
    }
    return sorted(terms)


def kpb_rank(a: int, b: int) -> str:
    m = min(a, b)
    if m == 0:
        return "zero-edge"
    if m == 1:
        return "star/planar"
    if m == 2:
        return "two-block/planar"
    return "K33-wall"


def first_reduced_factor_pair(n: int) -> tuple[int, int] | None:
    best: tuple[int, int] | None = None
    for a in range(1, isqrt(n) + 1):
        if n % a:
            continue
        b = n // a
        if gcd(a, b) != 1:
            continue
        cand = (a, b)
        if best is None or (cand[1], -cand[0]) < (best[1], -best[0]):
            best = cand
    return best


def perfect_product_rows(max_den: int) -> list[tuple[int, Fraction, str, int]]:
    rows = []
    for fr in farey_terms(max_den):
        if fr.numerator == 0:
            continue
        prod = fr.numerator * fr.denominator
        if sigma(prod) == 2 * prod:
            rows.append((prod, fr, kpb_rank(fr.numerator, fr.denominator), fr.denominator))
    return sorted(rows, key=lambda row: (row[0], row[3], row[1].numerator))


def euclid_euler_rows(exponents: list[int]) -> list[dict[str, object]]:
    rows = []
    for p in exponents:
        q = 2**p - 1
        a = 2 ** (p - 1)
        n = a * q
        rows.append(
            {
                "p": p,
                "a": a,
                "q": q,
                "mersenne_prime": is_prime(q),
                "N": n,
                "first_pair": first_reduced_factor_pair(n),
                "rank": kpb_rank(a, q),
                "farey_det_vs_half": q - 2 * a,
                "abundancy": Fraction(sigma(n), n),
            }
        )
    return rows


def unit_excess_prime_q_shadow(apex: int, a: int) -> tuple[int, Fraction, Fraction, bool]:
    """Return q=apex*a-1, prime-q abundancy model, defect from 2, primality.

    The formula is the actual abundancy of a*q when a is a power of two and
    q is prime.  Otherwise it is only the two-factor prime-q shadow.
    """

    q = apex * a - 1
    model = Fraction(apex * (2 * a - 1), q)
    defect = Fraction(2, 1) - model
    return q, model, defect, is_prime(q)


def tournament_from_criteria(
    vertices: list[str],
    criteria_orders: dict[str, list[str]],
    tie_path: list[str],
) -> dict[str, object]:
    pos = {name: {v: i for i, v in enumerate(order)} for name, order in criteria_orders.items()}
    tie_pos = {v: i for i, v in enumerate(tie_path)}
    edges: dict[tuple[str, str], str] = {}
    out = {v: set() for v in vertices}
    for a, b in combinations(vertices, 2):
        av = 0
        bv = 0
        votes = []
        for cname in criteria_orders:
            if pos[cname][a] < pos[cname][b]:
                av += 1
                votes.append(cname + ":A")
            else:
                bv += 1
                votes.append(cname + ":B")
        if av > bv:
            winner, loser = a, b
        elif bv > av:
            winner, loser = b, a
        elif tie_pos[a] < tie_pos[b]:
            winner, loser = a, b
        else:
            winner, loser = b, a
        edges[(winner, loser)] = ",".join(votes)
        out[winner].add(loser)
    scores = {v: len(out[v]) for v in vertices}
    c3 = 0
    cycles = []
    for a, b, c in combinations(vertices, 3):
        if b in out[a] and c in out[b] and a in out[c]:
            c3 += 1
            cycles.append((a, b, c))
        elif c in out[a] and b in out[c] and a in out[b]:
            c3 += 1
            cycles.append((a, c, b))

    sccs = strongly_connected_components(vertices, out)
    hp = hamiltonian_path_count(vertices, out)
    return {
        "out": out,
        "scores": scores,
        "c3": c3,
        "cycles": cycles,
        "sccs": sccs,
        "hamiltonian_paths": hp,
    }


def strongly_connected_components(vertices: list[str], out: dict[str, set[str]]) -> list[list[str]]:
    rev = {v: set() for v in vertices}
    for v, ws in out.items():
        for w in ws:
            rev[w].add(v)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in out[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return sorted(comps, key=lambda c: (-len(c), c))


def hamiltonian_path_count(vertices: list[str], out: dict[str, set[str]]) -> int:
    idx = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp: dict[tuple[int, int], int] = {}
    for v in vertices:
        dp[(1 << idx[v], idx[v])] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp.get((mask, last), 0)
            if not cur:
                continue
            v = vertices[last]
            for w in out[v]:
                j = idx[w]
                if mask & (1 << j):
                    continue
                dp[(mask | (1 << j), j)] = dp.get((mask | (1 << j), j), 0) + cur
    full = (1 << n) - 1
    return sum(dp.get((full, j), 0) for j in range(n))


def main() -> None:
    print("Farey perfect-product obstruction scout (codex S143)")
    print("=" * 64)
    print("Quotient declaration:")
    print("  fraction a/b -> product ab -> complete bipartite graph K_{a,b}")
    print("  preserves: reduced factor-pair incidence, Farey level b, K33 threshold")
    print("  destroys: graph isomorphism type, minor core identity, divisor lattice")
    print()

    print("[1] Ordinary F3/F4 product seam")
    for order in [3, 4]:
        print(f"F_{order} nonzero terms:")
        for fr in farey_terms(order):
            if fr.numerator == 0:
                continue
            a, b = fr.numerator, fr.denominator
            marker = ""
            prod = a * b
            if sigma(prod) == 2 * prod:
                marker += " PERFECT_PRODUCT"
            if kpb_rank(a, b) == "K33-wall":
                marker += " FIRST_K33_WALL" if fr == Fraction(3, 4) else " K33_WALL"
            print(
                f"  {fmt_fraction(fr):>3} -> product={prod:<2} "
                f"K_{{{a},{b}}} {kpb_rank(a,b)}{marker}"
            )
        print()

    print("[2] Perfect numbers as Farey-product rows")
    rows = euclid_euler_rows([2, 3, 5, 7, 13])
    for row in rows:
        a, b = row["a"], row["q"]
        pair = row["first_pair"]
        fr = Fraction(a, b)
        print(
            f"  p={row['p']:>2}  {fmt_fraction(fr):>8}  product={row['N']:<10} "
            f"first_F={pair[1] if pair else 'none':>5}  "
            f"det(q-2a)={row['farey_det_vs_half']:+d}  "
            f"rank={row['rank']:<16}  A={fmt_fraction(row['abundancy'])}"
        )
    print()
    print("  Reduced Farey terms with perfect-number product through denominator 140:")
    grouped: dict[int, list[tuple[Fraction, str, int]]] = defaultdict(list)
    for prod, fr, rank, level in perfect_product_rows(140):
        grouped[prod].append((fr, rank, level))
    for prod in sorted(grouped):
        first_level = min(level for _, _, level in grouped[prod])
        parts = [
            f"{fmt_fraction(fr)} in F_{level} ({rank})"
            for fr, rank, level in grouped[prod]
        ]
        print(f"    product={prod:<5} first_F={first_level:<3} terms: " + "; ".join(parts))
    print()

    print("[3] Unit-excess prime-q shadow: n=2 is exact perfection; n=14 is deficient")
    powers = [1, 2, 4, 8, 16, 64]
    for apex in [2, 14]:
        print(f"  apex n={apex}: fractions a/(n*a-1)")
        for a in powers:
            q, model, defect, qprime = unit_excess_prime_q_shadow(apex, a)
            actual = ""
            if qprime:
                actual_N = a * q
                actual = f" actual_A={fmt_fraction(Fraction(sigma(actual_N), actual_N))}"
            print(
                f"    a={a:<3} q={q:<4} prime={str(qprime):<5} "
                f"prime-q-A={fmt_fraction(model):>8}  2-A={fmt_fraction(defect):>8}{actual}"
            )
        print()
    print("  Identity: for a=2^k and prime q=n*a-1,")
    print("    sigma(a*q)/(a*q) = n(2a-1)/(na-1) = 2 + (2-n)/(na-1).")
    print("  Therefore exact perfection is forced by n=2; LRC14 has defect 12/(14a-1).")
    print()

    print("[4] Kuratowski/Wagner guardrails for the product ledger")
    print("  K5 edges = 10, but Farey term 2/5 gives K_{2,5} with 10 edges and it is planar.")
    print("  K3,3 edges = 9, but 3/3 is not a reduced Farey term; 1/9 is a star.")
    print("  The first reduced Farey complete-bipartite term containing K3,3 is 3/4 -> K_{3,4}.")
    print("  K5 disjoint-union K3,3 has (edges,vertices)=(19,11), the unreduced density")
    print("  mixture of (10,5) and (9,6), but it is nonplanar only because it already")
    print("  contains the two irreducible obstruction cores.  Minors/subdivisions compose;")
    print("  mediants and edge-density averages do not manufacture new forbidden minors.")
    print()

    print("[5] Tournament Analysis over carrier roles")
    vertices = [
        "minor_closure",
        "kuratowski_core",
        "K33_incidence",
        "perfect_abundancy",
        "unit_excess_chain",
        "farey_level",
        "product_edges",
        "raw_edge_density",
    ]
    criteria_orders = {
        "planarity_obstruction": [
            "minor_closure",
            "kuratowski_core",
            "K33_incidence",
            "raw_edge_density",
            "product_edges",
            "farey_level",
            "unit_excess_chain",
            "perfect_abundancy",
        ],
        "perfect_fixed_point": [
            "perfect_abundancy",
            "unit_excess_chain",
            "product_edges",
            "farey_level",
            "K33_incidence",
            "raw_edge_density",
            "kuratowski_core",
            "minor_closure",
        ],
        "lrc_farey_transfer": [
            "farey_level",
            "unit_excess_chain",
            "K33_incidence",
            "product_edges",
            "perfect_abundancy",
            "minor_closure",
            "kuratowski_core",
            "raw_edge_density",
        ],
        "anti_scalar_guard": [
            "minor_closure",
            "kuratowski_core",
            "K33_incidence",
            "perfect_abundancy",
            "unit_excess_chain",
            "farey_level",
            "product_edges",
            "raw_edge_density",
        ],
        "density_bookkeeping": [
            "raw_edge_density",
            "product_edges",
            "farey_level",
            "unit_excess_chain",
            "K33_incidence",
            "perfect_abundancy",
            "kuratowski_core",
            "minor_closure",
        ],
    }
    tie_path = [
        "minor_closure",
        "kuratowski_core",
        "K33_incidence",
        "perfect_abundancy",
        "unit_excess_chain",
        "farey_level",
        "product_edges",
        "raw_edge_density",
    ]
    tour = tournament_from_criteria(vertices, criteria_orders, tie_path)
    print("  Alternate vertex sets considered:")
    print("    Farey fractions, products, graph isomorphism types, minor cores,")
    print("    factor pairs, divisor atoms, density ratios, wall crossings, and proof obligations.")
    print("  Chosen vertices are carrier roles.  The quotient preserves proof-role pressure")
    print("  and destroys exact runner geometry and exact graph embeddings.")
    print("  Pairwise observable: earlier role in each criterion beats later role.")
    print("  Gauge: majority of five criteria; tie Hamiltonian path:")
    print("    " + " > ".join(tie_path))
    print("  Score histogram:", dict(sorted(Counter(tour["scores"].values()).items())))
    print("  Scores:")
    for v, sc in sorted(tour["scores"].items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"    {v:<20} {sc}")
    print(f"  directed_3_cycles={tour['c3']}")
    for cyc in tour["cycles"][:8]:
        print("    cycle " + " -> ".join(cyc) + f" -> {cyc[0]}")
    print("  SCCs:")
    for comp in tour["sccs"]:
        print("    {" + ", ".join(comp) + "}")
    print(f"  Hamiltonian-path count={tour['hamiltonian_paths']}")
    print()

    print("[6] Working hypotheses produced by this scout")
    leads = [
        "Even perfect numbers are the n=2 unit-excess product chain a/(2a-1) with a power of two and 2a-1 prime.",
        "LRC14's a/(14a-1) chain is the deficient prime-q shadow of that chain, with defect 12/(14a-1).",
        "F3/F4 is the first useful seam: 2/3 gives product 6 and K_{2,3}; 3/4 gives the first reduced K33 wall.",
        "The graph-product ledger is valid only with minor labels attached; edge counts and mediants are too lossy.",
        "Tournament vertices should be proof obligations/carrier roles here, not runners, because the question is which quotient preserves the obstruction predicate.",
    ]
    for i, lead in enumerate(leads, 1):
        print(f"  H{i}. {lead}")


if __name__ == "__main__":
    main()
