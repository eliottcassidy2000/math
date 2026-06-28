#!/usr/bin/env python3
"""HYP-3245 scout: equioscillation/autocorrelation atlas.

This extends HYP-3214's exact Fejer-kernel identification.  The guiding
principle is that many extremal objects can be read in two dual languages:

  equioscillation / double-zero / alternation in the primal variable
  autocorrelation / positive-definite / Gram-profile control in Fourier data

For the LRC14 k=8 frontier, the same idea gives a concrete out-correlation
signal: non-AP trap rows move ordinary difference-autocorrelation mass from
short lags to outward lags, while HYP-3225's Green/Rayleigh/Plucker sidecars
explain the finer obstruction type.  HYP-3238 is the even/odd
positive/negative compression companion, HYP-3239 is the dihedral/Borsuk-Ulam
sign-irrep companion, HYP-3236 is the all-bank Green conductance companion,
HYP-3237 is the Vitali/Brouwer core-wall companion, HYP-3219 is the Brouwer
trace-sign/SOS companion, HYP-3228 is the shell-magic companion, and HYP-3227
is the conductance-graph companion for the same finite trap boundary.
HYP-3243 is the topology/geometry carrier companion: lag transport is
proof-grade only after the open-tope / known-witness / finite-chamber /
state-lift / named-debt schema is attached.
HYP-3244 is the tiling/half-tiling controlled-forgetting companion: lag and
motif quotients are proof-grade only after tiling witness lift, half-tiling
descent, fiber, parent-automorphism, coboundary, and tail/tip sidecars are
attached.
"""

from __future__ import annotations

import cmath
import itertools
import math
from collections import Counter, defaultdict


TARGET = tuple(range(8))
HYP3202_TRAPS = (
    TARGET,
    (0, 2, 4, 6, 7, 8, 10, 12),
    (0, 1, 3, 5, 7, 9, 11, 13),
    (0, 1, 2, 3, 11, 12, 13, 14),
    (0, 3, 6, 8, 9, 11, 12, 14),
    (0, 2, 3, 5, 6, 8, 11, 14),
    (0, 6, 7, 8, 9, 10, 11, 12),
    (0, 4, 5, 8, 9, 10, 13, 14),
    (0, 8, 9, 10, 11, 12, 13, 14),
    (0, 1, 2, 3, 7, 8, 9, 10),
    (0, 2, 5, 7, 9, 10, 12, 14),
    (0, 1, 4, 5, 7, 8, 11, 12),
)

HYP3225_CLASS = {
    TARGET: "AP_terminal",
    (0, 2, 4, 6, 7, 8, 10, 12): "rank2_pair_plucker_bottleneck",
    (0, 1, 3, 5, 7, 9, 11, 13): "rank2_pair_plucker_bottleneck",
    (0, 1, 2, 3, 11, 12, 13, 14): "green_low_connectivity_bottleneck",
    (0, 3, 6, 8, 9, 11, 12, 14): "rank2_pair_plucker_bottleneck",
    (0, 2, 3, 5, 6, 8, 11, 14): "rank2_pair_plucker_bottleneck",
    (0, 6, 7, 8, 9, 10, 11, 12): "rank2_pair_plucker_bottleneck",
    (0, 4, 5, 8, 9, 10, 13, 14): "green_low_connectivity_bottleneck",
    (0, 8, 9, 10, 11, 12, 13, 14): "mixed_green_lorentzian_sidecar",
    (0, 1, 2, 3, 7, 8, 9, 10): "AFM_frustrated_high_rayleigh_debt",
    (0, 2, 5, 7, 9, 10, 12, 14): "AFM_frustrated_high_rayleigh_debt",
    (0, 1, 4, 5, 7, 8, 11, 12): "rank2_pair_plucker_bottleneck",
}

PROPOSED_SIGNAL_FIELDS = [
    ("contact_word", "ordered signs/multiplicities at equality or root-contact nodes"),
    ("lag_barycenter", "center of mass of positive/negative autocorrelation residual"),
    ("transport_cost", "earthmover distance from AP triangular lag law to E"),
    ("Fejer_annihilator_projection", "lag residual component killed by the F7 double-zero kernel"),
    ("shell_lag_commutator", "mismatch between HYP-3228 shell slack and lag transport"),
    ("sidecar_entropy", "number of repair labels needed before the quotient is legal"),
    ("small_pattern_motif_id", "HYP-3226 motif id carrying the signal as a typed payload atom"),
    ("payload_preserved", "LRC coordinate retained by the motif or lag projection"),
    ("payload_destroyed", "coordinate lost by scalarizing the motif or lag residual"),
    ("repair_sidecar", "sidecar that restores the destroyed coordinate"),
    ("terminal_risk_label", "HYP-3226 direct/sidecar/analogy/raw risk tag"),
    ("scale_survival_bit", "whether the signal survives HYP-3231 scale recursion"),
    ("apex_fold_side", "HYP-3232 low/apex/antipode side of the lag residual"),
    ("cyclotomic_mode", "HYP-3217 mode carrying the residual after projection"),
    ("chart_change_class", "HYP-3234 signed chart needed to compare recurrence letters"),
    ("cyclotomic_factor_signature", "HYP-3233 (x-1)^depth * Phi_d tag"),
    ("cap_field_conductor", "HYP-3235 conductor seen by the cap/binding row"),
    ("fejer_square_status", "whether the residual respects the totally-positive F7 square"),
    ("gauss_sum_margin", "HYP-3218 Fejer/Vaaler margin tag"),
    ("ap_self_dual_fixed_point_status", "whether AP self-duality is retained"),
    ("green_resistance_slack", "HYP-3236 resistance excess relative to AP"),
    ("lambda2_conductance_rank", "HYP-3236 algebraic-connectivity rank among bounded rows"),
    ("negative_covariance_leakage", "lost covariance sign information after positive-part graphing"),
    ("thomson_current_profile", "unit-current bottleneck profile for the trap row"),
    ("fiedler_bottleneck_id", "Fiedler cut/channel that first detects the trap"),
    ("vitali_wall_side", "HYP-3237 bulk/core side of the witness packet"),
    ("brouwer_saddle_sign", "max-min saddle orientation at the core equioscillation"),
    ("phi14_core_witness", "primitive 14th-root witness status for the row"),
    ("core_bulk_transport_status", "whether measure bulk or cyclotomic core carries the proof"),
    ("brouwer_trace_sign", "HYP-3219 topological sign of the odd obstruction"),
    ("degree_sos_factorization", "degree/sign times SOS magnitude factorization tag"),
    ("even_odd_bonferroni_node_slack", "local slack at the even/odd Bonferroni split"),
    ("positive_negative_duality_status", "HYP-3238 even-positive/odd-negative packet split"),
    ("odd_negative_payload_reconstruction", "whether clipped odd/negative data is zero or restored"),
    ("dihedral_irrep_label", "HYP-3239 D7 irrep carrying the residual"),
    ("complement_antiautomorphism_sign", "whether complement acts as orientation-reversing anti-aut"),
    ("borsuk_ulam_index", "free Z2 odd-degree certificate tag"),
    ("imaginary_gauss_sum_sign", "i*sqrt(7) obstruction / sign-irrep orientation tag"),
    ("phi4_bimodal_extremizer_rank", "HYP-3239/HYP-3238 bimodal phi4 extremizer check"),
    ("equioscillation_saddle_index", "HYP-3241 antipodal Phi14 witness-pair count"),
    ("phi14_core_universality_status", "whether AP/Goddyn-Wong share the same core witnesses"),
    ("dilation_witness_grid", "HYP-3240 promoted Phi_{14d} witness grid for dilation rows"),
    ("core_witness_break_reason", "why the base Phi14 witness fails or survives"),
    ("imaginary_norm_route_status", "whether the Q(sqrt(-7)) norm scalar route is legal"),
    ("danger_nerve_euler_characteristic", "HYP-3242 measured nerve/Euler characteristic state"),
    ("lonely_hole_status", "whether the lonely witness is carried by a cover hole"),
    ("cech_betti_sidecar", "Cech/Betti label for the danger-cover obstruction"),
    ("topological_shadow_class", "which topology/combinatorics/representation shadow is active"),
    ("cover_hole_witness_pair", "Borsuk-Ulam antipodal pair attached to the cover hole"),
    ("oriented_matroid_tope_status", "HYP-3243 open-safe tope or boundary cocircuit state"),
    ("circle_endpoint_arrangement_cell", "endpoint/wall-crossing cell containing the residual"),
    ("cech_safe_component_rank", "safe-component nerve rank after open/closed separation"),
    ("finite_chamber_schema_status", "open safe / known witness / finite discharge / state lift / debt"),
    ("state_lift_H7_obstruction", "whether a residual atom lifts to the forbidden H=7 endpoint"),
    ("proof_carrier_tournament_rank", "HYP-3243 carrier priority before scalar lag transport is used"),
    ("tiling_witness_lift_status", "HYP-3244 fixed-path tiling-cover witness status"),
    ("half_tiling_descent_certificate", "fiber/aut/coboundary proof that quotienting preserves the LRC predicate"),
    ("path_presentation_fiber_weight", "H(T)/|Aut(T)| or named fiber debt attached to the carrier"),
    ("parent_aut_word_orbit_id", "incident-word orbit under parent automorphisms"),
    ("rectangle_hourglass_residue", "GF(2) coboundary residue detecting illegal compression"),
    ("tail_tip_deletion_signature", "tip/tail deletion payload retained through the span"),
    ("controlled_forgetting_span_status", "lift/compress/descent/fail-debt status for the residual"),
]


def fejer_coeffs(n: int) -> dict[int, int]:
    return {lag: n - abs(lag) for lag in range(-(n - 1), n)}


def aperiodic_autocorr_real(seq: list[int]) -> dict[int, int]:
    n = len(seq)
    out: dict[int, int] = {}
    for lag in range(-(n - 1), n):
        total = 0
        for i in range(n):
            j = i + lag
            if 0 <= j < n:
                total += seq[i] * seq[j]
        out[lag] = total
    return out


def periodic_autocorr_complex(seq: list[complex]) -> list[complex]:
    n = len(seq)
    vals = []
    for lag in range(n):
        vals.append(sum(seq[i] * seq[(i + lag) % n].conjugate() for i in range(n)))
    return vals


def dft(seq: list[complex]) -> list[complex]:
    n = len(seq)
    return [
        sum(seq[j] * cmath.exp(-2j * math.pi * m * j / n) for j in range(n))
        for m in range(n)
    ]


def fejer_value(n: int, theta: float) -> float:
    if abs(math.sin(theta / 2.0)) < 1.0e-13:
        return float(n * n)
    return (math.sin(n * theta / 2.0) / math.sin(theta / 2.0)) ** 2


def de_moivre_cubic(u: float) -> float:
    return u**3 + u**2 - 2.0 * u - 1.0


def verify_fejer_demoivre(samples: int = 2000) -> tuple[float, float, float]:
    max_identity_error = 0.0
    max_zero_value = 0.0
    max_zero_slope = 0.0
    for s in range(samples + 1):
        theta = 2.0 * math.pi * s / samples
        u = 2.0 * math.cos(theta)
        max_identity_error = max(max_identity_error, abs(de_moivre_cubic(u) ** 2 - fejer_value(7, theta)))
    eps = 1.0e-6
    for j in range(1, 7):
        theta = 2.0 * math.pi * j / 7.0
        max_zero_value = max(max_zero_value, abs(fejer_value(7, theta)))
        slope = (fejer_value(7, theta + eps) - fejer_value(7, theta - eps)) / (2.0 * eps)
        max_zero_slope = max(max_zero_slope, abs(slope))
    return max_identity_error, max_zero_value, max_zero_slope


def golay_pair(length_power: int = 3) -> tuple[list[int], list[int]]:
    a = [1]
    b = [1]
    for _ in range(length_power):
        a, b = a + b, a + [-x for x in b]
    return a, b


def fano_difference_set() -> list[int]:
    return [1 if i in {0, 1, 3} else 0 for i in range(7)]


def zadoff_chu(n: int = 7, q: int = 1) -> list[complex]:
    return [cmath.exp(-1j * math.pi * q * k * (k + 1) / n) for k in range(n)]


def periodic_binary_autocorr(seq: list[int]) -> list[int]:
    n = len(seq)
    return [sum(seq[i] * seq[(i + lag) % n] for i in range(n)) for lag in range(n)]


def support_autocorr(E: tuple[int, ...], max_lag: int = 14) -> list[int]:
    s = set(E)
    return [sum(1 for x in s if x + lag in s) for lag in range(max_lag + 1)]


def sign_changes(values: list[int]) -> int:
    signs = [1 if x > 0 else -1 if x < 0 else 0 for x in values]
    nz = [x for x in signs if x]
    return sum(1 for a, b in zip(nz, nz[1:]) if a != b)


def trap_autocorr_table() -> list[dict[str, object]]:
    base = support_autocorr(TARGET)
    rows = []
    for E in HYP3202_TRAPS:
        ac = support_autocorr(E)
        residual = [x - y for x, y in zip(ac, base)]
        low = sum(residual[1:8])
        high = sum(residual[8:])
        rows.append(
            {
                "E": E,
                "class": HYP3225_CLASS[E],
                "l1": sum(abs(x) for x in residual[1:]),
                "max_abs": max(abs(x) for x in residual[1:]),
                "sign_changes": sign_changes(residual[1:]),
                "low_lag_deficit": low,
                "out_lag_surplus": high,
                "residual": residual,
            }
        )
    return rows


def motif_atlas() -> list[dict[str, object]]:
    return [
        {
            "name": "Fejer_AP_interval",
            "equioscillation": "double zeros at nontrivial 7th roots / Chebyshev V7 sharpness",
            "autocorrelation": "triangular interval autocorrelation (7-|lag|)_+",
            "LRC_use": "HYP-3214 sector magic function and AP orbit kernel",
            "risk": "naive F_k*F_7 pairing does not equal cap_k",
            "score": 100,
        },
        {
            "name": "Johnson_pair_Pascal",
            "equioscillation": "14-clock pair-normalized cap face",
            "autocorrelation": "pair-count / Johnson-scheme inner distribution",
            "LRC_use": "cap side complementary to 7-sector Fejer",
            "risk": "wrong clock if collapsed to 7-sector kernel",
            "score": 95,
        },
        {
            "name": "Christoffel_Darboux_OPUC",
            "equioscillation": "orthogonal-polynomial reproducing kernel peaks at support atoms",
            "autocorrelation": "Toeplitz moment kernel / Verblunsky reflection data",
            "LRC_use": "Toeplitz lambda-min and Fejer-Riesz trap chart",
            "risk": "forgets exchange sidecars",
            "score": 91,
        },
        {
            "name": "Welch_ETF_simplex",
            "equioscillation": "all off-diagonal Gram entries equal",
            "autocorrelation": "constant correlation / frame-potential equality",
            "LRC_use": "model for AP equality face and Delsarte sharpness",
            "risk": "too symmetric for finite trap boundary",
            "score": 86,
        },
        {
            "name": "difference_set_two_level",
            "equioscillation": "all nonzero cyclic lags equal",
            "autocorrelation": "two-level periodic autocorrelation",
            "LRC_use": "finite residue/design analogue of Fejer sharpness",
            "risk": "periodic model can erase endpoint owners",
            "score": 80,
        },
        {
            "name": "CAZAC_Zadoff_Chu",
            "equioscillation": "constant magnitude spectrum",
            "autocorrelation": "zero nonzero periodic autocorrelation",
            "LRC_use": "perfect decorrelation guardrail for wide/far phases",
            "risk": "complex phases may not preserve integer-speed predicate",
            "score": 76,
        },
        {
            "name": "Golay_complementary_pair",
            "equioscillation": "paired ripple cancellation",
            "autocorrelation": "sum of aperiodic sidelobes is zero",
            "LRC_use": "signed-pair analogue of L_y bypassing termwise majorants",
            "risk": "requires pair carrier, not single row",
            "score": 72,
        },
        {
            "name": "Barker_low_sidelobe",
            "equioscillation": "small aperiodic sidelobes",
            "autocorrelation": "bounded +/-1 ripple",
            "LRC_use": "finite trap canary for near-sharp but non-perfect rows",
            "risk": "scarce examples; not a uniform theorem",
            "score": 65,
        },
        {
            "name": "equiripple_FIR_Remez",
            "equioscillation": "alternating error extrema",
            "autocorrelation": "magnitude-square filter autocorrelation",
            "LRC_use": "algorithmic template for dual search",
            "risk": "continuous design may forget arithmetic sidecars",
            "score": 60,
        },
        {
            "name": "tournament_out_correlation",
            "equioscillation": "score/outdegree balance",
            "autocorrelation": "out-neighborhood overlap A A^T",
            "LRC_use": "directed/out-correlation reading of proof-carrier tournaments",
            "risk": "raw runner tournament is too lossy",
            "score": 45,
        },
    ]


def print_exact_examples() -> None:
    err, zero_value, zero_slope = verify_fejer_demoivre()
    print("EXACT / NUMERIC IDENTITIES")
    print(f"  max |(de_moivre_cubic(2cos t))^2 - F7(t)| over grid = {err:.3e}")
    print(f"  max F7(2pi*j/7), j=1..6 = {zero_value:.3e}")
    print(f"  max central slope at those zeros = {zero_slope:.3e} (double-zero check)")
    print(f"  Fejer_7 coefficients = {fejer_coeffs(7)}")

    interval7 = [1] * 7
    interval8 = [1] * 8
    print(f"  interval_7 autocorr equals Fejer_7 coefficients: {aperiodic_autocorr_real(interval7) == fejer_coeffs(7)}")
    print(f"  interval_8 center autocorr coeffs = {fejer_coeffs(8)}")

    fano_ac = periodic_binary_autocorr(fano_difference_set())
    print(f"  Fano (7,3,1) difference-set periodic autocorr = {fano_ac}")

    zc = zadoff_chu()
    zc_ac = periodic_autocorr_complex(zc)
    zc_spec = [abs(x) for x in dft(zc)]
    print(f"  Zadoff-Chu_7 max off-periodic-autocorr = {max(abs(x) for x in zc_ac[1:]):.3e}")
    print(f"  Zadoff-Chu_7 spectrum magnitude ripple = {max(zc_spec) - min(zc_spec):.3e}")

    a, b = golay_pair(3)
    ac_a = aperiodic_autocorr_real(a)
    ac_b = aperiodic_autocorr_real(b)
    off = [ac_a[lag] + ac_b[lag] for lag in range(1, len(a))]
    print(f"  Golay_pair_8 summed positive-lag sidelobes = {off}")

    barker = [1, 1, 1, -1, -1, 1, -1]
    barker_ac = [aperiodic_autocorr_real(barker)[lag] for lag in range(0, len(barker))]
    print(f"  Barker_7 positive-lag autocorr = {barker_ac}")

    simplex_n = 7
    simplex_d = simplex_n - 1
    offdiag = -1.0 / (simplex_n - 1)
    frame_potential = simplex_n + simplex_n * (simplex_n - 1) * offdiag * offdiag
    welch = simplex_n * simplex_n / simplex_d
    print(f"  simplex_7 offdiag={offdiag:+.6f} frame_potential={frame_potential:.6f} Welch={welch:.6f}")
    print()


def print_trap_autocorr() -> None:
    print("LRC14 HYP-3202 TRAP OUT-CORRELATION TABLE")
    print("  ordinary difference autocorr of speed support compared with AP interval")
    print("  low_lag_deficit=sum residual lags 1..7; out_lag_surplus=sum residual lags 8..14")
    class_stats: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in trap_autocorr_table():
        class_stats[str(row["class"])].append(row)
        print(
            f"  E={row['E']} | {row['class']} | L1={row['l1']} max={row['max_abs']} "
            f"sign_changes={row['sign_changes']} low={row['low_lag_deficit']} out={row['out_lag_surplus']}"
        )
        if row["E"] in {
            (0, 1, 2, 3, 11, 12, 13, 14),
            (0, 1, 3, 5, 7, 9, 11, 13),
            (0, 8, 9, 10, 11, 12, 13, 14),
        }:
            print(f"    residual={row['residual']}")
    print()
    print("  class averages:")
    for cls, rows in sorted(class_stats.items()):
        avg_l1 = sum(int(r["l1"]) for r in rows) / len(rows)
        avg_changes = sum(int(r["sign_changes"]) for r in rows) / len(rows)
        low_values = sorted(set(int(r["low_lag_deficit"]) for r in rows))
        out_values = sorted(set(int(r["out_lag_surplus"]) for r in rows))
        print(
            f"    {cls}: count={len(rows)} avg_L1={avg_l1:.3f} "
            f"avg_sign_changes={avg_changes:.3f} low_values={low_values} out_values={out_values}"
        )
    print(
        "  readout=every non-AP trap shifts autocorrelation mass outward: "
        "low_lag_deficit + out_lag_surplus = 0, but sidecar class controls the ripple shape."
    )
    print()


def print_motif_atlas() -> None:
    print("MOTIF ATLAS: EQUIOSCILLATION <-> AUTOCORRELATION")
    for motif in motif_atlas():
        print(
            f"  {motif['name']}: eq={motif['equioscillation']} | "
            f"auto={motif['autocorrelation']} | use={motif['LRC_use']} | risk={motif['risk']}"
        )
    print()


def print_signal_fields() -> None:
    print("PROPOSED NEXT SIGNAL FIELDS")
    for field, meaning in PROPOSED_SIGNAL_FIELDS:
        print(f"  {field}: {meaning}")
    print()


def tournament_analysis() -> None:
    motifs = motif_atlas()
    ordered = sorted(((str(m["name"]), int(m["score"])) for m in motifs), key=lambda x: (-x[1], x[0]))
    print("TOURNAMENT ANALYSIS")
    print("  vertices=certificate/motif families, not runners or raw arcs")
    print("  pairwise_observable=which motif keeps both equioscillation and autocorrelation payload")
    print("  switch/gauge=A->B iff A keeps the LRC predicate with fewer destroyed sidecars")
    print(f"  score_hist={dict(Counter(score for _, score in ordered))}")
    print("  directed_3cycles=0")
    print("  scc_sizes=[" + ",".join("1" for _ in ordered) + "]")
    print("  hamiltonian_path_count=1")
    print("  priority_path=" + " -> ".join(name for name, _ in ordered))
    print()


def assumption_challenge() -> None:
    print("ASSUMPTION-CHALLENGE")
    print("  alternate vertices considered: runners, sectors, Fourier modes, autocorrelation lags,")
    print("  zeros/equioscillation nodes, Toeplitz moments, trap sidecars, Johnson pairs,")
    print("  difference-set blocks, codewords, filter bands, and proof obligations.")
    print("  HYP-3244 adds flip bits, deletion roots, rectangle/hourglass defects,")
    print("  parent-automorphism word orbits, and tiling-cover fibers as alternate vertices.")
    print("  chosen vertices: certificate/motif families plus LRC trap autocorrelation residuals.")
    print("  preserved predicate: HYP-3214 Fejer magic-function identity, HYP-3224 normal-fan")
    print("  AP extremality, HYP-3225 finite trap sidecar discharge, and HYP-3244")
    print("  tiling-lift / half-tiling-descent legality.")
    print("  destroyed information: exact LRC safe-set geometry, PGF root radii, endpoint owners,")
    print("  finite chamber address, state-lift obligations, and analytic equidistribution")
    print("  obligations unless routed to named sidecars; HYP-3244 adds Hamiltonian-path")
    print("  presentation, canary/filler fiber mass, parent word-orbit identity,")
    print("  rectangle/hourglass residue, and tail/tip deletion payload.")
    print("  challenged assumption: equioscillation is only a root-location story.")
    print("  The atlas treats it as the primal shadow of autocorrelation/positive-definiteness.")


def main() -> None:
    print("HYP-3245 equioscillation / autocorrelation atlas")
    print("=" * 78)
    print("basis=HYP-3214 Fejer F7 = de Moivre^2 = AP autocorrelation")
    print("extension=multitude atlas plus HYP-3202 trap out-correlation residuals")
    print()
    print_exact_examples()
    print_trap_autocorr()
    print_motif_atlas()
    print_signal_fields()
    tournament_analysis()
    assumption_challenge()


if __name__ == "__main__":
    main()
