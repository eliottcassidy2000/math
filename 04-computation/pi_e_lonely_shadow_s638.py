#!/usr/bin/env python3
"""S638: lonely-shadow extensions of the pi/e quadratic carrier.

S636 used the Vieta carrier

    T^2 - S*T + P,  with S=e+pi and P=e*pi,

to prove that S and P cannot both be algebraic, and that adding
D=e-pi makes any two of {S,P,D} reconstruct e and pi.  This addendum records
two further ways the carrier bites:

1. If log(pi) were algebraic, Lindemann-Weierstrass would force S, P, and D all
   transcendental.  Therefore any algebraic S/P/D would also prove log(pi)
   transcendental.
2. If one elementary symmetric shadow descended to Qbar, transverse symmetric
   shadows become forbidden.  For example, if S were algebraic then every power
   sum e^k+pi^k for k>=2 would be transcendental; if P were algebraic then
   every e^k+pi^k for k>=1 would be transcendental.

The script prints exact power-sum polynomials, a small PSLQ sanity check, and
a proof-route Tournament Analysis whose vertices are theorem routes rather than
numbers.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations

import mpmath as mp


mp.mp.dps = 110
E = mp.e
PI = mp.pi
S = E + PI
P = E * PI
D = E - PI
L = mp.log(PI)


Poly = dict[tuple[int, int], int]


def clean(poly: Poly) -> Poly:
    return {mono: coeff for mono, coeff in poly.items() if coeff}


def add(a: Poly, b: Poly, scale_b: int = 1) -> Poly:
    out = dict(a)
    for mono, coeff in b.items():
        out[mono] = out.get(mono, 0) + scale_b * coeff
    return clean(out)


def shift(poly: Poly, ds: int = 0, dp: int = 0, scale: int = 1) -> Poly:
    return clean({(i + ds, j + dp): scale * coeff for (i, j), coeff in poly.items()})


def power_sum_polys(n: int) -> list[Poly]:
    """Return p_k(S,P)=e^k+pi^k for k=0..n.

    Newton recurrence for two roots:

        p_0 = 2, p_1 = S, p_k = S*p_{k-1} - P*p_{k-2}.
    """

    polys: list[Poly] = [{(0, 0): 2}, {(1, 0): 1}]
    for _k in range(2, n + 1):
        polys.append(add(shift(polys[-1], ds=1), shift(polys[-2], dp=1), scale_b=-1))
    return polys[: n + 1]


def poly_degree(poly: Poly, var: str) -> int:
    idx = 0 if var == "S" else 1
    return max((mono[idx] for mono in poly), default=0)


def format_poly(poly: Poly) -> str:
    terms = []
    for (s_pow, p_pow), coeff in sorted(
        poly.items(), key=lambda item: (sum(item[0]), item[0][0], item[0][1]), reverse=True
    ):
        if coeff == 0:
            continue
        sign = "-" if coeff < 0 else "+"
        mag = abs(coeff)
        bits = []
        if s_pow:
            bits.append("S" if s_pow == 1 else f"S^{s_pow}")
        if p_pow:
            bits.append("P" if p_pow == 1 else f"P^{p_pow}")
        mono = "*".join(bits)
        if not mono:
            chunk = str(mag)
        elif mag == 1:
            chunk = mono
        else:
            chunk = f"{mag}*{mono}"
        terms.append((sign, chunk))
    if not terms:
        return "0"
    first_sign, first_chunk = terms[0]
    out = f"-{first_chunk}" if first_sign == "-" else first_chunk
    for sign, chunk in terms[1:]:
        out += f" {sign} {chunk}"
    return out


def scalar_relation(x: mp.mpf, degree: int = 3, maxcoeff: int = 10**7):
    vals = [mp.mpf(1)]
    for _ in range(degree):
        vals.append(vals[-1] * x)
    return mp.pslq(vals, tol=mp.mpf("1e-85"), maxcoeff=maxcoeff)


def format_relation(coeffs, names: list[str]) -> str:
    if coeffs is None:
        return "none"
    chunks = []
    for coeff, name in zip(coeffs, names):
        if coeff:
            chunks.append(f"{coeff}*{name}")
    return " + ".join(chunks).replace("+ -", "- ") + " = 0"


@dataclass(frozen=True)
class Route:
    name: str
    unconditional: int
    strength: int
    repo_fit: int
    computable: int
    risk: int


ROUTES = [
    Route("log_pi_lindemann_gate", 5, 5, 4, 4, 1),
    Route("lonely_power_sum_fallout", 5, 4, 5, 5, 1),
    Route("transverse_symmetric_polynomial_template", 5, 4, 5, 4, 1),
    Route("quadratic_trace_norm_vieta", 5, 3, 5, 5, 1),
    Route("h21_side_condition_analogy", 4, 3, 5, 4, 1),
    Route("known_e_to_pi_gelfond_constant", 3, 2, 3, 3, 3),
    Route("pslq_height_sieve", 2, 1, 3, 5, 2),
    Route("schanuel_full_completion", 1, 5, 4, 3, 4),
    Route("raw_rational_case_split", 1, 1, 2, 5, 4),
]

CRITERIA = ["unconditional", "strength", "repo_fit", "computable"]


def route_tournament():
    idx = {route.name: i for i, route in enumerate(ROUTES)}
    edges: dict[tuple[int, int], int] = {}
    wins = Counter({route.name: 0 for route in ROUTES})
    for a, b in combinations(ROUTES, 2):
        av = 0
        bv = 0
        for criterion in CRITERIA:
            if getattr(a, criterion) > getattr(b, criterion):
                av += 1
            elif getattr(b, criterion) > getattr(a, criterion):
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


def beats(edges: dict[tuple[int, int], int], a: int, b: int) -> bool:
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
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
    lines: list[str] = []
    lines.append("S638 Pi/e Lonely Shadow Addendum")
    lines.append("================================")
    lines.append("")
    lines.append("Carrier constants")
    lines.append("-----------------")
    lines.append(f"S=e+pi        = {mp.nstr(S, 35)}")
    lines.append(f"P=e*pi        = {mp.nstr(P, 35)}")
    lines.append(f"D=e-pi        = {mp.nstr(D, 35)}")
    lines.append(f"L=log(pi)     = {mp.nstr(L, 35)}")
    lines.append("")
    lines.append("New unconditional implications")
    lines.append("------------------------------")
    lines.append(
        "1. If L=log(pi) were algebraic, Lindemann-Weierstrass would make "
        "1, e=e^1, and pi=e^L linearly independent over Qbar.  The exponents "
        "0, 1, L are distinct because pi is neither 1 nor e.  Therefore S=e+pi "
        "and D=e-pi could not be algebraic.  Also P=e*pi=e^(1+L) would be "
        "transcendental by Hermite-Lindemann."
    )
    lines.append(
        "   Contrapositive: algebraicity of any one of S, P, or D would prove "
        "log(pi) transcendental."
    )
    lines.append(
        "2. If one elementary symmetric shadow descends, transverse symmetric "
        "shadows are forbidden.  If S is algebraic, then every p_k=e^k+pi^k "
        "for k>=2 is transcendental.  If P is algebraic, then every p_k for "
        "k>=1 is transcendental."
    )
    lines.append("")
    lines.append("Power-sum carrier p_k=e^k+pi^k in S,P")
    lines.append("--------------------------------------")
    polys = power_sum_polys(10)
    for k, poly in enumerate(polys[1:], start=1):
        lines.append(
            f"p_{k:<2} = {format_poly(poly):<48} "
            f"deg_S={poly_degree(poly, 'S')} deg_P={poly_degree(poly, 'P')}"
        )
    lines.append("")
    lines.append("Why the power-sum fallout is a theorem")
    lines.append("--------------------------------------")
    lines.append(
        "- With S algebraic, p_k=f_k(S,P) algebraic for any k>=2 would make P "
        "algebraic over Qbar because f_k(S,X) is nonconstant in X.  Then "
        "T^2-S*T+P makes e and pi algebraic, impossible."
    )
    lines.append(
        "- With P algebraic, p_k=f_k(S,P) algebraic for any k>=1 would make S "
        "algebraic over Qbar because f_k(X,P) is nonconstant in X.  Again "
        "Vieta reconstructs e and pi."
    )
    lines.append(
        "- More generally: if S is algebraic then every symmetric polynomial "
        "F(S,P) nonconstant in P is transcendental; if P is algebraic then "
        "every symmetric polynomial nonconstant in S is transcendental."
    )
    lines.append("")
    lines.append("Rationality fallout table")
    lines.append("-------------------------")
    rows = [
        (
            "Assume S=e+pi algebraic",
            "P, D, log(pi), and all p_k for k>=2 are transcendental",
        ),
        (
            "Assume P=e*pi algebraic",
            "S, D, log(pi), and all p_k for k>=1 are transcendental",
        ),
        (
            "Assume D=e-pi algebraic",
            "S, P, log(pi), and alternating power sums e^k+(-pi)^k for k>=2 are transcendental",
        ),
        (
            "Assume log(pi) algebraic",
            "S, P, and D are all transcendental without Schanuel",
        ),
        (
            "Assume Schanuel for 1 and i*pi",
            "e and pi are algebraically independent, so every nonconstant polynomial shadow is transcendental",
        ),
    ]
    for left, right in rows:
        lines.append(f"- {left}: {right}.")
    lines.append("")
    lines.append("PSLQ height sanity checks")
    lines.append("------------------------")
    scalars = [("S", S), ("P", P), ("D", D), ("L", L)]
    for name, val in scalars:
        lines.append(
            f"degree<=3 scalar relation for {name}, coeff<=1e7: "
            f"{format_relation(scalar_relation(val), ['1', name, name+'^2', name+'^3'])}"
        )
    p_values = [E**k + PI**k for k in range(1, 6)]
    rel = mp.pslq([mp.mpf(1), *p_values], tol=mp.mpf("1e-80"), maxcoeff=10**7)
    lines.append(
        "linear relation among 1,p_1,...,p_5, coeff<=1e7: "
        f"{format_relation(rel, ['1', 'p1', 'p2', 'p3', 'p4', 'p5'])}"
    )
    sheet_rel = mp.pslq([S * S, D * D, P], tol=mp.mpf("1e-90"), maxcoeff=100)
    lines.append(
        "structural sheet PSLQ on [S^2,D^2,P]: "
        f"{format_relation(sheet_rel, ['S^2', 'D^2', 'P'])}"
    )
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append(
        "Challenged vertex assumption: the useful tournament vertices are not "
        "S, P, D, or numeric guesses.  I considered shadows, assumptions, "
        "implication arrows, theorem machines, branch sheets, and repo carrier "
        "analogies; the chosen vertices are proof routes.  The quotient "
        "preserves field-descent obstruction data and destroys decimal "
        "near-miss information."
    )
    edges, wins = route_tournament()
    n = len(ROUTES)
    order = sorted(range(n), key=lambda i: wins[ROUTES[i].name], reverse=True)
    score_hist = Counter(wins.values())
    lines.append(f"vertices={n}")
    lines.append(f"score_hist={dict(sorted(score_hist.items()))}")
    lines.append(f"directed_3cycles={directed_3cycles(edges, n)}")
    lines.append(f"Hamiltonian_paths={hamiltonian_count(edges, n)}")
    lines.append("ranking:")
    for rank, idx in enumerate(order, start=1):
        route = ROUTES[idx]
        lines.append(f"  {rank}. {route.name} (score={wins[route.name]})")
    lines.append("")
    lines.append("Repo transfer")
    lines.append("-------------")
    lines.append(
        "The new pattern mirrors the incoming H=21 window closure: scalar H=21 "
        "is decided only after retaining strongness and c3 side conditions.  "
        "For pi/e, scalar S or P is also too thin; the retained side conditions "
        "are log-commensurability, the discriminant sheet, and transverse "
        "symmetric-polynomial fallout."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
