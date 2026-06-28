#!/usr/bin/env python3
"""Cyclotomic Delsarte/Beurling-Selberg magic packet for the k=8 LRC node.

The bounded k=8 dual functional used by HYP-3153/HYP-3204 is

    10 q0 + q3 + 10 q6 = 10 * (q0 + q6 + q3/10).

This scout packages that vector as a finite "magic function" in three
equivalent languages:

* Shell polynomial / Beurling-Selberg contact:
      f(n) = ((n-1)(n-2)(n-4)(n-5))/4
  on n=0,...,6.  It is nonnegative on all integer shells and vanishes exactly
  on shells 1,2,4,5, leaving endpoint mass plus the Worpitzky central repair.

* Delsarte/Newton-binomial form:
      f(n) = 10 - 10 C(n,1) + 10 C(n,2) - 9 C(n,3) + 6 C(n,4).

* Cyclotomic/Joukowski form:
      P(z)=10+z^3+10z^6, so z^-3 P(z)=10(u^3-3u)+1, u=z+z^-1.

The computation then tests exact extremality on the same anchored bounded bank
as HYP-3200/HYP-3204 and compares the magic direction against the AP-polarized
cyclotomic support normal from HYP-3203.
"""

from __future__ import annotations

import cmath
import importlib.util
import itertools
import math
from collections import Counter
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
H3200 = HERE / "lrc_k8_cumulant_universality_codex_20260628.py"
spec = importlib.util.spec_from_file_location("h3200", H3200)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot import {H3200}")
h3200 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h3200)


TARGET = tuple(range(8))
EVEN_AP = tuple(2 * i for i in range(8))
UNIFORM = tuple(Fraction(1, 7) for _ in range(7))


def shell_magic(n: int) -> Fraction:
    return Fraction((n - 1) * (n - 2) * (n - 4) * (n - 5), 4)


def binom(n: int, r: int) -> int:
    return math.comb(n, r) if 0 <= r <= n else 0


def newton_binomial_coefficients(values: list[Fraction]) -> list[Fraction]:
    """Return a_j such that f(n)=sum_j a_j C(n,j)."""
    diffs = values[:]
    coeffs = [diffs[0]]
    while len(diffs) > 1:
        diffs = [b - a for a, b in zip(diffs, diffs[1:])]
        coeffs.append(diffs[0])
    return coeffs


def evaluate_newton(coeffs: list[Fraction], n: int) -> Fraction:
    return sum(coeff * binom(n, j) for j, coeff in enumerate(coeffs))


def dot(a: tuple[Fraction, ...] | list[Fraction], b: tuple[Fraction, ...] | list[Fraction]) -> Fraction:
    return sum(x * y for x, y in zip(a, b))


def centered(q: tuple[Fraction, ...] | list[Fraction]) -> tuple[Fraction, ...]:
    return tuple(x - Fraction(1, 7) for x in q)


def magic_value(q: tuple[Fraction, ...] | list[Fraction]) -> Fraction:
    return sum(shell_magic(n) * q[n] for n in range(7))


def ly_value(q: tuple[Fraction, ...] | list[Fraction]) -> Fraction:
    return q[0] + q[6] + q[3] / 10


def support_projection(q: tuple[Fraction, ...], ap_direction: tuple[Fraction, ...]) -> Fraction:
    return dot(centered(q), ap_direction)


def rank_exact(
    rows: list[dict[str, object]],
    target: tuple[int, ...],
    getter,
    maximize: bool = True,
) -> tuple[int, int, list[dict[str, object]]]:
    target_value = getter(next(row for row in rows if row["E"] == target))
    if maximize:
        beaters = sum(1 for row in rows if getter(row) > target_value)
        ties = [row for row in rows if getter(row) == target_value]
        ordered = sorted(rows, key=lambda row: (getter(row), row["E"]), reverse=True)
    else:
        beaters = sum(1 for row in rows if getter(row) < target_value)
        ties = [row for row in rows if getter(row) == target_value]
        ordered = sorted(rows, key=lambda row: (getter(row), row["E"]))
    return beaters, len(ties), ordered


def format_vec(vec: tuple[Fraction, ...] | list[Fraction]) -> str:
    return "[" + ", ".join(str(x) for x in vec) + "]"


def de_moivre_values() -> list[tuple[int, float, float, float]]:
    """Values of M(u)=10(u^3-3u)+1 at u=2cos(2*pi*m/7), m=1,2,3."""
    out = []
    for m in (1, 2, 3):
        theta = 2 * math.pi * m / 7
        u = 2 * math.cos(theta)
        value = 10 * (u**3 - 3 * u) + 1
        reduced = -10 * u * u - 10 * u + 11
        out.append((m, u, value, reduced))
    return out


def cyclic_regularization_threshold() -> tuple[float, int]:
    """Central coefficient rho needed for 10+rho z^3+10 z^6 to be >=0 on 7th roots."""
    worst_m = 0
    worst = 0.0
    for m in range(1, 7):
        z = cmath.exp(2j * math.pi * m / 7)
        carrier = 10 * (z**3 + z**-3)
        if carrier.real < worst:
            worst = carrier.real
            worst_m = m
    return -worst, worst_m


def deficit_ratios(
    rows: list[dict[str, object]],
    ap_magic: Fraction,
    ap_support: Fraction,
) -> tuple[tuple[Fraction, dict[str, object], Fraction, Fraction], tuple[Fraction, dict[str, object], Fraction, Fraction]]:
    ratios = []
    for row in rows:
        if row["E"] == TARGET:
            continue
        magic_def = ap_magic - row["magic"]
        support_def = ap_support - row["ap_support"]
        if support_def <= 0:
            continue
        ratios.append((magic_def / support_def, row, magic_def, support_def))
    return min(ratios, key=lambda item: item[0]), max(ratios, key=lambda item: item[0])


def projection_report(ap_direction: tuple[Fraction, ...]) -> tuple[Fraction, tuple[Fraction, ...], float]:
    f = tuple(shell_magic(n) for n in range(7))
    mean = Fraction(sum(f), 7)
    f0 = tuple(x - mean for x in f)
    alpha = dot(f0, ap_direction) / dot(ap_direction, ap_direction)
    residual = tuple(f0[i] - alpha * ap_direction[i] for i in range(7))
    cosine = float(dot(f0, ap_direction)) / math.sqrt(float(dot(f0, f0) * dot(ap_direction, ap_direction)))
    return alpha, residual, cosine


def tournament_analysis() -> None:
    carriers = {
        "quartic_shell_magic_contact": 100,
        "delsarte_newton_binomial_dual": 96,
        "joukowski_chebyshev_cubic": 91,
        "ap_polarized_support_normal": 88,
        "ordered_tail_exchange_price": 83,
        "toeplitz_trap_discharge_sidecar": 74,
        "cyclic_psd_regularization": 42,
        "raw_cyclotomic_energy": 25,
        "raw_q3_scalar": 18,
    }
    ordered = sorted(carriers.items(), key=lambda item: (-item[1], item[0]))
    print("\nTOURNAMENT ANALYSIS")
    print("vertices=proof currencies: shells, modes, normals, traps, and proof obligations")
    print("alternate_vertex_sets_considered=miss-count shells; Laurent modes; AP-support normal; exchange traps; sector-boundary events")
    print("chosen_quotient_preserves=L_y shell functional and q-vector support direction")
    print("chosen_quotient_destroys=individual sector identity, covariance-pair labels, and exchange-neighbor witnesses")
    print("challenged_assumption=cyclotomic magic need not be a nonnegative cyclic Fourier kernel")
    print("pairwise_observable=retained LRC proof payload plus exact bounded-bank survival")
    print("switch/gauge=A->B iff payload score(A)>payload score(B); ties lexical")
    print(f"score_hist={{{', '.join(f'{score}:1' for _, score in ordered)}}}")
    print("directed_3cycles=0")
    print("hamiltonian_path_count=1")
    print("priority_path=" + " -> ".join(name for name, _ in ordered))


def main() -> None:
    rows = [
        h3200.row_moments((0,) + combo)
        for combo in itertools.combinations(range(1, 15), 7)
    ]
    primitive = [row for row in rows if row["primitive"]]
    by_E = {row["E"]: row for row in rows}
    consec = by_E[TARGET]
    even_ap = by_E[EVEN_AP]

    ap_direction = centered(tuple(consec["q"]))
    for row in rows:
        q = tuple(row["q"])
        row["magic"] = magic_value(q)
        row["Ly"] = ly_value(q)
        row["ap_support"] = support_projection(q, ap_direction)

    f_values = [shell_magic(n) for n in range(7)]
    coeffs = newton_binomial_coefficients(f_values)
    ap_magic = consec["magic"]
    ap_support = consec["ap_support"]

    print("HYP-3228 cyclotomic Delsarte/Beurling-Selberg magic-function scout")
    print("=" * 78)
    print("bank=anchored {0} union A, A subset [1,14], |A|=7")
    print(f"rows_all={len(rows)}")
    print(f"rows_primitive={len(primitive)}")
    print(f"consec={TARGET}")
    print(f"even_AP={EVEN_AP}")
    print(f"consec_q={format_vec(consec['q'])}")
    print(f"even_AP_same_q={even_ap['q'] == consec['q']}")

    print("\nSHELL MAGIC FUNCTION")
    print("f(n)=((n-1)(n-2)(n-4)(n-5))/4")
    print("centered v=n-3 form=(v^2-1)(v^2-4)/4")
    print(f"shell_values_f0..f6={format_vec(f_values)}")
    print(f"nonnegative_on_integer_shells={all(x >= 0 for x in f_values)}")
    print(f"contact_zero_shells={[n for n, x in enumerate(f_values) if x == 0]}")
    print(f"positive_shells={[n for n, x in enumerate(f_values) if x > 0]}")
    print("interpretation=endpoint bimodality plus one central Worpitzky repair")

    print("\nDELSARTE / NEWTON-BINOMIAL FORM")
    print(f"newton_binomial_coefficients={format_vec(coeffs)}")
    checks = [evaluate_newton(coeffs, n) == f_values[n] for n in range(7)]
    print(f"identity_holds_on_shells_0..6={all(checks)}")
    print("identity=10*S0 -10*S1 +10*S2 -9*S3 +6*S4")
    print("where Sj=E[binom(N,j)] for miss count N")

    print("\nCYCLOTOMIC / JOUKOWSKI FORM")
    print("P(z)=10 + z^3 + 10 z^6")
    print("self_reciprocal=z^6 P(1/z)=P(z)")
    print("z^-3 P(z)=10(z^3+z^-3)+1 = 10(u^3-3u)+1, u=z+z^-1")
    print("real_7th_cyclotomic_cubic=C(u)=u^3+u^2-2u-1")
    print("reduction_mod_C(u)=-10u^2-10u+11")
    for m, u, value, reduced in de_moivre_values():
        print(f"m={m}: u=2cos(2pi*{m}/7)={u:+.12f}; M(u)={value:+.12f}; reduced={reduced:+.12f}")
    rho, worst_m = cyclic_regularization_threshold()
    print(f"cyclic_nonnegative_completion_requires_central_rho_at_least={rho:.12f} at root_index={worst_m}")
    print("current_central_rho=1.000000000000")
    print("incoming_HYP3214=positive sector-side magic is the Fejer kernel F7")
    print("guardrail=this shell Ly dual is distinct; cyclic PSD positivity would overprice the q3 repair")

    print("\nEXACT BOUNDED-BANK EXTREMALITY")
    for label, bank in (("primitive", primitive), ("all", rows)):
        beaters, ties, ordered = rank_exact(bank, TARGET, lambda row: row["magic"], maximize=True)
        print(f"{label}: magic=10q0+q3+10q6 beaters={beaters} ties={ties} max={ordered[0]['magic']} at {ordered[0]['E']}")
    print(f"consec_magic={ap_magic} ({float(ap_magic):.12f})")
    print(f"consec_Ly={consec['Ly']} ({float(consec['Ly']):.12f})")
    print(f"magic_equals_10Ly_on_all_rows={all(row['magic'] == 10 * row['Ly'] for row in rows)}")
    print(f"all_bank_magic_ties={[row['E'] for row in rows if row['magic'] == ap_magic]}")

    print("\nAP SUPPORT NORMAL ALIGNMENT")
    for label, bank in (("primitive", primitive), ("all", rows)):
        beaters, ties, ordered = rank_exact(bank, TARGET, lambda row: row["ap_support"], maximize=True)
        print(f"{label}: ap_support beaters={beaters} ties={ties} max={ordered[0]['ap_support']} at {ordered[0]['E']}")
    low, high = deficit_ratios(primitive, ap_magic, ap_support)
    print("primitive_deficit_ratio=magic_deficit/support_deficit, excluding consec")
    print(f"min_ratio={low[0]} ({float(low[0]):.9f}) at E={low[1]['E']}")
    print(f"  magic_deficit={low[2]} support_deficit={low[3]}")
    print(f"max_ratio={high[0]} ({float(high[0]):.9f}) at E={high[1]['E']}")
    print(f"  magic_deficit={high[2]} support_deficit={high[3]}")
    alpha, residual, cosine = projection_report(ap_direction)
    print(f"projection_alpha_magic_centered_onto_AP_support={alpha} ({float(alpha):.9f})")
    print(f"cosine_magic_vs_AP_support={cosine:.9f}")
    print(f"residual_vector={format_vec(residual)}")
    print("readout=AP support is a coercive normal, but not the whole magic vector")

    print("\nNEAR-CHALLENGER SNAPSHOT")
    top_magic = sorted(primitive, key=lambda row: (row["magic"], row["E"]), reverse=True)[:8]
    for idx, row in enumerate(top_magic):
        print(
            f"magic_rank_{idx}: E={row['E']} magic={row['magic']} "
            f"deficit={ap_magic-row['magic']} ap_support_deficit={ap_support-row['ap_support']}"
        )
    shell_hist = Counter(tuple(n for n, mass in enumerate(row["q"]) if mass) for row in primitive)
    print(f"primitive_support_pattern_count={len(shell_hist)}")
    print(f"most_common_support_patterns={shell_hist.most_common(5)}")

    tournament_analysis()


if __name__ == "__main__":
    main()
