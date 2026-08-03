#!/usr/bin/env python3
"""Exact companion for THM-3222's universal reset upper-filter collar.

The script reuses the maintained THM-3110 signed product-Gamma banks and
THM-3120 reduced row fractions.  For every one of the 115 supports and both
banks it proves that the distinguished Ferrers alphabet Q is a submultiset
of the reduced pole alphabet P.  At the first live degree five, common-prefix
lambda subtraction then gives, for every nonempty tau <= P-Q,

    G_5^(Q+tau)((5)) = -Phi(h_5) * sum_(r in tau) r^5 < 0.

Only exact integer arithmetic is used.  The all-tau quantifier is algebraic;
the finite run checks its bank hypotheses and a broad set of direct virtual-
alphabet controls in every case.
"""

import ast
from collections import Counter
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
ROW_SCRIPT = ROOT / "04-computation/gmc_product_gamma_row_pole_flag_thm3120.py"
SCHUR_SCRIPT = ROOT / (
    "04-computation/gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
)
THM3110 = ROOT / (
    "01-canon/theorems/"
    "THM-3110-arbitrary-anchored-product-gamma-dominant-tail-and-low-"
    "histogram-reduction.md"
)
THM3120 = ROOT / (
    "01-canon/theorems/"
    "THM-3120-row-pole-prefix-newton-flag-positivity.md"
)

PINS = {
    ROW_SCRIPT: "93828710e859ef2215ee164708b835d339594d241dcd70a090050ad3187901dd",
    SCHUR_SCRIPT: "15b94691d53afbcdc6aefda89fc7cd5497534ca70fb780a686892dabb5646d6f",
    THM3110: "97113139be9ea16b1434ca9bc076920700920b20d0c08b6fb082ca06c6af48c4",
    THM3120: "399f031ce40c3e0e6caa2da3c81e94068c8286ba9013a5fe7631984af9f76d8b",
}


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for pinned_path, expected_hash in PINS.items():
    require(lf_hash(pinned_path) == expected_hash,
            "dependency pin drift: " + pinned_path.name)

syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)),
        "optimization-sensitive assert found")
require(not any(isinstance(node, ast.Constant) and isinstance(node.value, float)
                for node in ast.walk(syntax)), "floating literal found")

spec = spec_from_file_location("thm3120_companion", ROW_SCRIPT)
row_flag = module_from_spec(spec)
spec.loader.exec_module(row_flag)


def dominant_row(invariant):
    """THM-3110's unique Ferrers-dominant row Q_1 or Q_2."""

    if invariant == 0:
        return ((0, 3), (2, 0), (2, 0), (1, 0))
    require(invariant == 1, "unknown bank")
    return ((0, 3), (1, 1), (2, 0), (1, 0), (1, 0))


def complete_virtual(positive, negative, degree):
    """Return h_degree[positive-negative] by truncated products."""

    series = [1] + [0] * degree
    for root in positive:
        old = series
        series = [sum(old[n - k] * root**k for k in range(n + 1))
                  for n in range(degree + 1)]
    for root in negative:
        old = series
        series = [old[n] - (root * old[n - 1] if n else 0)
                  for n in range(degree + 1)]
    return series[degree]


def direct_prefix_response(invariant, a, b, q_roots, tau, c5):
    """Evaluate the prefixed degree-five singleton response from the bank."""

    sigma = tuple(q_roots) + tuple(tau)
    phi_h = 0
    phi_p = 0
    for coefficient, row in row_flag.g.BANKS[invariant]:
        roots = row_flag.g.residual_roots(invariant, row, a, b)
        phi_h += coefficient * complete_virtual(roots, sigma, 5)
        phi_p += coefficient * (
            sum(root**5 for root in roots) - sum(root**5 for root in sigma)
        )
    q_h = complete_virtual(q_roots, sigma, 5)
    q_p = sum(root**5 for root in q_roots) - sum(root**5 for root in sigma)
    require(phi_h == c5, "complete-row common-prefix invariance failed")
    require(phi_p == 0, "power-row common-prefix invariance failed")
    response = phi_h * q_p - phi_p * q_h
    expected = -c5 * sum(root**5 for root in tau)
    require(response == expected and response < 0,
            "direct upper-filter response formula failed")
    return response


records = []
direct_controls = 0
full_nonzero_h_controls = 0
case_payload = []
smallest_tangent = None
largest_full_magnitude = 0

for a, b in row_flag.support_universe():
    for invariant in (0, 1):
        numerator, denominator, remaining, _, _, _ = (
            row_flag.reduced_row_fraction(invariant, a, b)
        )
        require(numerator[:5] == [0] * 5 and numerator[5] > 0,
                "first-live numerator coefficient failed")
        c5 = numerator[5]
        base_series = row_flag.series_coefficients(numerator, denominator, 5)
        require(base_series[:5] == (0,) * 5 and base_series[5] == c5,
                "first-live complete rows failed")

        # Signed row marginal cancellation gives Phi(1)=Phi(p_j)=0.
        require(sum(coefficient for coefficient, _
                    in row_flag.g.BANKS[invariant]) == 0,
                "signed mass failed")
        for degree in range(1, 6):
            power = sum(
                coefficient * sum(root**degree for root in
                                  row_flag.g.residual_roots(
                                      invariant, row, a, b))
                for coefficient, row in row_flag.g.BANKS[invariant]
            )
            require(power == 0, "signed power marginal failed")

        p_counter = Counter(remaining)
        q_roots = tuple(sorted(
            root for root in row_flag.g.residual_roots(
                invariant, dominant_row(invariant), a, b) if root
        ))
        q_counter = Counter(q_roots)
        require(q_counter <= p_counter,
                "dominant reset alphabet is not physically available")
        r_counter = p_counter - q_counter
        remainder = tuple(sorted(r_counter.elements()))
        require(len(remainder) >= 5 and all(root > 0 for root in remainder),
                "physical reset collar lost its positive remainder")

        state_count = prod(multiplicity + 1
                           for multiplicity in r_counter.values())
        records.append((a, b, invariant + 1, sum(p_counter.values()),
                        len(q_roots), len(remainder), c5, state_count))

        # Direct controls include every distinct tangent, a five-letter case,
        # an alternating completion, and the full completion.
        tau_controls = {(root,) for root in r_counter}
        tau_controls.add(tuple(remainder[:5]))
        tau_controls.add(tuple(remainder[::2]))
        tau_controls.add(remainder)
        for tau in sorted(tau_controls):
            direct_prefix_response(invariant, a, b, q_roots, tau, c5)
            direct_controls += 1

        full_h = complete_virtual((), remainder, 5)
        require(full_h < 0, "full negative alphabet h5 control vanished")
        full_nonzero_h_controls += 1

        tangent = min(c5 * root**5 for root in r_counter)
        smallest_tangent = tangent if smallest_tangent is None else min(
            smallest_tangent, tangent)
        largest_full_magnitude = max(
            largest_full_magnitude,
            c5 * sum(root**5 for root in remainder),
        )
        case_payload.append(
            f"{a},{b},I{invariant + 1}|P={tuple(sorted(p_counter.elements()))}"
            f"|Q={q_roots}|R={remainder}|c={c5}|states={state_count}"
        )

require(len(records) == 230, "case census drift")
require(min(record[3] for record in records) == 11
        and max(record[3] for record in records) == 133, "pole depth range")
require(min(record[4] for record in records) == 6
        and max(record[4] for record in records) == 62, "reset depth range")
require(min(record[5] for record in records) == 5
        and max(record[5] for record in records) == 71, "collar width range")
require(min(record[6] for record in records) == 32
        and max(record[6] for record in records) == 410879700000,
        "first-live coefficient range")
require(len({record[6] for record in records}) == 225,
        "first-live distinct-value census")
require(min(record[7] for record in records) == 16
        and max(record[7] for record in records) == 187406683791040512,
        "filter state range")
require(sum(record[7] for record in records) == 476976523236717752,
        "aggregate filter-state census")
require(smallest_tangent == 32, "sharp tangent response")
require(full_nonzero_h_controls == 230, "full h5 hostile census")

case_digest = sha256("\n".join(case_payload).encode("ascii")).hexdigest()
print("dependency_pins=4:PASS")
print("support_universe=115:banks=2:cases=230")
print("physical_reset_alphabet=Q<=P:230/230")
print("depths=P:11..133,Q:6..62,R:5..71")
print("first_live_c5=positive:230/230:min32:max410879700000:distinct225")
print("filter_states=min16:max187406683791040512:sum476976523236717752")
print("common_prefix_identity=G5(Q+tau,(5))=-c5*sum(tau^5)")
print("direct_virtual_alphabet_controls=%d:PASS" % direct_controls)
print("full_nonzero_h5_negative_alphabet_controls=230/230")
print("sharp_nonempty_response_magnitude=32")
print("largest_full_response_magnitude=%d" % largest_full_magnitude)
print("case_digest_sha256=%s" % case_digest)
print("scope=principal_upper_filters_only:not_global_all_state_selector_cones")
print("status=PASS")
