#!/usr/bin/env python3
"""S636: the quadratic carrier T^2 - S*T + P around e and pi.

This continues S635.  The new unconditional gain is to add the
antisymmetric shadow

    D = e - pi.

Any two of S=e+pi, P=e*pi, and D=e-pi being algebraic would force both
e and pi algebraic.  Therefore at least two of the three shadows are
transcendental.  This is still far from proving which of e+pi or e*pi is
irrational, but it is a sharper carrier statement than the usual two-shadow
lemma.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

import mpmath as mp


mp.mp.dps = 100
E = mp.e
PI = mp.pi
S = E + PI
P = E * PI
D = E - PI


def scalar_relation(x: mp.mpf, degree: int, maxcoeff: int = 10**8):
    vec = [mp.mpf(1)]
    for _ in range(degree):
        vec.append(vec[-1] * x)
    return mp.pslq(vec, tol=mp.mpf("1e-80"), maxcoeff=maxcoeff)


def monomials_2(x: mp.mpf, y: mp.mpf, total_degree: int):
    mons = []
    vals = []
    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            mons.append((i, j))
            vals.append((x**i) * (y**j))
    return mons, vals


def monomials_3(x: mp.mpf, y: mp.mpf, z: mp.mpf, total_degree: int):
    mons = []
    vals = []
    for i in range(total_degree + 1):
        for j in range(total_degree + 1 - i):
            for k in range(total_degree + 1 - i - j):
                mons.append((i, j, k))
                vals.append((x**i) * (y**j) * (z**k))
    return mons, vals


def pair_relation(x: mp.mpf, y: mp.mpf, degree: int, maxcoeff: int = 10**7):
    mons, vals = monomials_2(x, y, degree)
    rel = mp.pslq(vals, tol=mp.mpf("1e-75"), maxcoeff=maxcoeff)
    return mons, rel


def triple_relation(x: mp.mpf, y: mp.mpf, z: mp.mpf, degree: int, maxcoeff: int = 10**7):
    mons, vals = monomials_3(x, y, z, degree)
    rel = mp.pslq(vals, tol=mp.mpf("1e-75"), maxcoeff=maxcoeff)
    return mons, rel


def format_poly_2(mons, coeffs, names=("X", "Y")) -> str:
    if coeffs is None:
        return "none"
    terms = []
    for coeff, (i, j) in zip(coeffs, mons):
        if coeff == 0:
            continue
        bits = []
        if i:
            bits.append(names[0] if i == 1 else f"{names[0]}^{i}")
        if j:
            bits.append(names[1] if j == 1 else f"{names[1]}^{j}")
        mono = "*".join(bits) if bits else "1"
        terms.append(f"{coeff}*{mono}")
    return " + ".join(terms).replace("+ -", "- ")


def format_poly_3(mons, coeffs, names=("S", "P", "D")) -> str:
    if coeffs is None:
        return "none"
    terms = []
    for coeff, (i, j, k) in zip(coeffs, mons):
        if coeff == 0:
            continue
        bits = []
        if i:
            bits.append(names[0] if i == 1 else f"{names[0]}^{i}")
        if j:
            bits.append(names[1] if j == 1 else f"{names[1]}^{j}")
        if k:
            bits.append(names[2] if k == 1 else f"{names[2]}^{k}")
        mono = "*".join(bits) if bits else "1"
        terms.append(f"{coeff}*{mono}")
    return " + ".join(terms).replace("+ -", "- ")


@dataclass(frozen=True)
class Route:
    name: str
    unconditional: int
    strength: int
    connects_repo: int
    computable: int
    risk: int


ROUTES = [
    Route("three_shadow_pair_obstruction", 5, 4, 4, 5, 1),
    Route("quadratic_trace_norm_carrier", 5, 4, 4, 5, 1),
    Route("schanuel_ai_completion", 1, 5, 4, 3, 4),
    Route("lrc_clock_descent_analogy", 3, 3, 5, 4, 2),
    Route("exponential_lift_relations", 2, 3, 4, 3, 3),
    Route("log_commensurability_ln_pi", 1, 3, 5, 3, 4),
    Route("pslq_height_sieve", 2, 1, 3, 5, 2),
    Route("raw_decimal_near_misses", 1, 1, 2, 5, 5),
]

CRITERIA = ["unconditional", "strength", "connects_repo", "computable"]


def route_tournament():
    idx = {r.name: i for i, r in enumerate(ROUTES)}
    edges = {}
    wins = {r.name: 0 for r in ROUTES}
    for a, b in combinations(ROUTES, 2):
        av = 0
        bv = 0
        for c in CRITERIA:
            if getattr(a, c) > getattr(b, c):
                av += 1
            elif getattr(b, c) > getattr(a, c):
                bv += 1
        if a.risk < b.risk:
            av += 1
        elif b.risk < a.risk:
            bv += 1
        winner = a if av >= bv else b
        i, j = idx[a.name], idx[b.name]
        edges[(i, j)] = 1 if winner is a else 0
        wins[winner.name] += 1
    return edges, wins


def beats(edges, a, b):
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_count(edges, n):
    @lru_cache(maxsize=None)
    def dp(mask, last):
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if (prev_mask >> prev) & 1 and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def directed_3cycles(edges, n):
    total = 0
    for a, b, c in combinations(range(n), 3):
        out = Counter()
        for i, j in [(a, b), (a, c), (b, c)]:
            winner = i if beats(edges, i, j) else j
            out[winner] += 1
        if sorted(out.values()) == [1, 1, 1]:
            total += 1
    return total


def main() -> str:
    lines = []
    lines.append("S636 Quadratic pi/e Carrier")
    lines.append("===========================")
    lines.append("")
    lines.append("Definitions")
    lines.append("-----------")
    lines.append(f"S = e + pi = {mp.nstr(S, 30)}")
    lines.append(f"P = e*pi   = {mp.nstr(P, 30)}")
    lines.append(f"D = e - pi = {mp.nstr(D, 30)}")
    lines.append("")
    lines.append("Unconditional sharpening")
    lines.append("------------------------")
    lines.append("Any two of S, P, D reconstruct e and pi algebraically:")
    rows = [
        ("S,D", "e=(S+D)/2, pi=(S-D)/2"),
        ("S,P", "e and pi are roots of T^2 - S*T + P"),
        ("D,P", "e and -pi are roots of T^2 - D*T - P"),
    ]
    for pair, formula in rows:
        lines.append(f"- {pair}: {formula}.")
    lines.append(
        "Therefore no pair among {S,P,D} can both be algebraic, since that "
        "would force e and pi algebraic.  Hence at least two of S, P, D are "
        "transcendental."
    )
    lines.append("")
    lines.append("Discriminant carrier")
    lines.append("--------------------")
    disc = S * S - 4 * P
    lines.append(f"S^2 - 4P = {mp.nstr(disc, 30)}")
    lines.append(f"D^2      = {mp.nstr(D * D, 30)}")
    lines.append(
        "So D is the hidden square-root side channel of the trace/norm pair. "
        "The quadratic carrier is not only S and P; it is S, P, and the "
        "discriminant sheet."
    )
    lines.append("")
    lines.append("Conditional completion")
    lines.append("----------------------")
    lines.append(
        "If Schanuel's conjecture holds, then e and pi are algebraically "
        "independent.  Then every nonconstant polynomial f(e,pi) with "
        "algebraic coefficients is transcendental; otherwise f(e,pi)=alpha "
        "would give the nonzero relation f(X,Y)-alpha."
    )
    lines.append(
        "Under that conditional, S, P, D are all transcendental, and the pairs "
        "(S,P), (S,D), and (P,D) are algebraically independent."
    )
    lines.append("")
    lines.append("Exponential lift obstruction attempts")
    lines.append("-------------------------------------")
    lines.append(
        "If S were algebraic, then exp(S)=exp(e)*exp(pi).  Since exp(pi) is "
        "known transcendental and exp(S) is transcendental for nonzero "
        "algebraic S, this gives no contradiction without a theorem about "
        "the algebraic independence of exp(e), exp(pi), or exp(S)."
    )
    lines.append(
        "If P were algebraic, then pi=P/e and exp(pi)=exp(P/e).  Current "
        "standard theorems do not rule this out because the exponent P/e "
        "is not algebraic."
    )
    lines.append(
        "So the direct attack still needs algebraic independence, not just "
        "individual transcendence."
    )
    lines.append("")
    lines.append("PSLQ height sieve")
    lines.append("-----------------")
    lines.append(
        "This is non-proof evidence only.  It asks whether small integer "
        "polynomial relations are visible at 100 digits."
    )
    for name, value in [("S", S), ("P", P), ("D", D)]:
        hit = None
        for deg in range(1, 9):
            rel = scalar_relation(value, deg)
            if rel is not None:
                hit = (deg, rel)
                break
        if hit is None:
            lines.append(f"- {name}: no scalar relation degree <= 8, maxcoeff <= 1e8.")
        else:
            lines.append(f"- {name}: relation at degree {hit[0]}: {hit[1]}")
    for label, x, y in [("(S,P)", S, P), ("(S,D)", S, D), ("(P,D)", P, D)]:
        mons, rel = pair_relation(x, y, 4)
        if rel is None:
            lines.append(f"- {label}: no pair relation total degree <= 4, maxcoeff <= 1e7.")
        else:
            lines.append(f"- {label}: {format_poly_2(mons, rel)}")
    mons3, rel3 = triple_relation(S, P, D, 2)
    lines.append(f"- (S,P,D), degree <= 2: {format_poly_3(mons3, rel3)}")
    lines.append("")
    lines.append("Repo interpretation")
    lines.append("-------------------")
    lines.append(
        "The repo's old pi/e comma thread asks whether ln(pi) is rational: "
        "multiplicative commensurability.  S636 is the additive/quadratic "
        "sibling: whether trace S or norm P descends to Qbar."
    )
    lines.append(
        "In LRC terms, S is like a reset/trace clock, P is like a norm or "
        "multiplicative shell, and D is the discriminant sheet that records "
        "which runner/constant lives on which branch.  Forgetting D is exactly "
        "the kind of observer-blind collapse S634 warned against."
    )
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append(
        "Vertices are proof routes.  Pairwise observable: which route gives a "
        "stronger and less risky obstruction while connecting back to repo "
        "carriers.  Switch/gauge: majority vote over unconditionality, strength, "
        "repo connectivity, computability, and lower risk.  Tie Hamiltonian path "
        "is the route list order."
    )
    edges, wins = route_tournament()
    ranked = sorted((r.name for r in ROUTES), key=lambda n: (-wins[n], n))
    for name in ranked:
        lines.append(f"  {wins[name]:2d}  {name}")
    lines.append(f"score histogram: {dict(sorted(Counter(wins.values()).items()))}")
    lines.append(f"directed 3-cycles: {directed_3cycles(edges, len(ROUTES))}")
    lines.append(f"Hamiltonian paths: {hamiltonian_count(edges, len(ROUTES))}")
    lines.append("")
    lines.append("Assumption challenge")
    lines.append("--------------------")
    lines.append(
        "Alternate vertices considered: constants e/pi, shadows S/P/D, "
        "polynomial relations, exponential lifts, log commensurability, "
        "height-sieve searches, LRC clocks, relation lattices, and proof routes. "
        "This lab chooses proof routes as vertices."
    )
    lines.append(
        "Preserved predicate: whether a route can rule out algebraic descent. "
        "Destroyed data: exact branch identity, individual constants, and "
        "which shadow carries transcendence."
    )
    lines.append("")
    lines.append("Next tests")
    lines.append("----------")
    lines.append(
        "1. Generalize the three-shadow lemma: for a generically finite map, "
        "which coordinate pairs force algebraic descent of the hidden source?"
    )
    lines.append(
        "2. Add discriminant-sheet labels to LRC reset-period ledgers: trace/norm "
        "quotients need an observer branch."
    )
    lines.append(
        "3. Search repo constants for other triples where any two shadows reconstruct "
        "a forbidden algebraic source."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
