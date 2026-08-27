#!/usr/bin/env python3
"""Exact audit for the fixed-gauge two-reversal presentation sphere.

Imports the audited literal subset-DP primitives from THM-4224.  All new
all-n formula checks below are exact integer identities.  The finite universe
is every two-reversal presentation and every unordered owner pair for
3 <= n <= 14.
"""

from __future__ import annotations

from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
DEP = ROOT / "04-computation/tournament_order11_tail_fiber_thm4224.py"
EXPECTED_DEP_SHA256 = "3c9ad5de462dcc71d4312bba12b42086d4e1e43520b5209b4a807f762d059c18"


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


need(sha256(DEP.read_bytes()).hexdigest() == EXPECTED_DEP_SHA256, "dependency hash")
spec = spec_from_file_location("thm4224_dp", DEP)
need(spec is not None and spec.loader is not None, "import specification")
dp = module_from_spec(spec)
sys.modules["thm4224_dp"] = dp
spec.loader.exec_module(dp)


def shell_formula(n: int):
    z = n - 1
    answer = set()
    for edge in dp.all_arcs(n):
        if edge not in ((0, z), (0, 1), (z - 1, z)):
            answer.add(tuple(sorted(((0, z), edge))))
    for b in range(1, z):
        for c in range(1, b + 1):
            answer.add(tuple(sorted(((0, b), (c, z)))))
    return frozenset(answer)


def type_b_parameters(reversals, z: int):
    if (0, z) in reversals:
        return None
    first = [edge for edge in reversals if edge[0] == 0]
    second = [edge for edge in reversals if edge[1] == z]
    need(len(first) == len(second) == 1, "type-B parse")
    return first[0][1], second[0][0]


def constants(b: int, c: int):
    """Return (h,K,eta,gamma,lambda) for the type-B tail formulas."""
    need(1 <= c <= b, "constant range")
    K = 2 * 3 ** (b - 1) + 2 ** (b + 1) - 2
    gamma = 2 ** max(b - c - 1, 0)
    if b == 1:
        return 1, K, 2, gamma, 1
    h = 2 ** (b - 1) + 1

    def eta_at(d: int) -> int:
        need(1 <= d <= b, "eta range")
        if d == b:
            return h
        return 2**d * 3 ** (b - d - 1) + 2 ** (b - d) - (1 if d == 1 else 0)

    eta = eta_at(c)
    lam = 2**b - 1 if c == b else eta_at(b - c)
    return h, K, eta, gamma, lam


def gap_t_n_1(n: int) -> int:
    need(n >= 5, "T(n,1) formula range")
    x = 2 ** (n - 4)
    y = 3 ** (n - 4)
    return 6 * (17 * x * y + 2 * y - 7 * x * x + 4 * x - 3)


digest = sha256()
grand_presentations = 0
grand_pair_rows = 0
grand_formula_rows = 0

print("two_reversal_sphere_audit")
print("status=PROVED_CLASSIFICATION+SYMBOLIC_TYPE_B_TAIL+FINITE_EXACT_FULL_SPHERE")
print(f"dependency_sha256={EXPECTED_DEP_SHA256}")

for n in range(3, 15):
    z = n - 1
    arcs = dp.all_arcs(n)
    brute = frozenset(
        tuple(sorted(reversals))
        for reversals in combinations(arcs, 2)
        if dp.is_strong(dp.from_reversals(n, reversals))
    )
    classified = shell_formula(n)
    need(brute == classified, f"strong classification n={n}")
    need(len(brute) == (n - 1) ** 2 - 3, f"strong count n={n}")

    product_failures = []
    plus_min_failures = []
    candidate_failures = []
    rho_failures = []
    five_failures = []
    five_equalities = []
    minimum_gap = None
    minimum_gap_rows = []
    max_b = -1
    max_p = 1
    max_rows = []
    formula_rows = 0
    pair_rows = 0

    for reversals in sorted(classified):
        tournament = dp.from_reversals(n, reversals)
        data = dp.path_data(tournament)
        profile = dp.tournament_profile(tournament, data)
        pairs = dp.all_pair_data(tournament, data, profile)
        pair_rows += len(pairs)

        if profile.gap < 0:
            five_failures.append(reversals)
        if profile.gap == 0:
            five_equalities.append(reversals)
        if minimum_gap is None or profile.gap < minimum_gap:
            minimum_gap = profile.gap
            minimum_gap_rows = [(reversals, profile)]
        elif profile.gap == minimum_gap:
            minimum_gap_rows.append((reversals, profile))

        by_pair = {(row.i, row.j): row for row in pairs}
        for row in pairs:
            tag = (reversals, (row.i, row.j), row.b, row.product, row.minimum, row.rho)
            if row.b > row.product:
                product_failures.append(tag)
            if row.b > row.product + row.minimum:
                plus_min_failures.append(tag)
            if 25 * row.b > 27 * row.product:
                candidate_failures.append(tag)
            if row.b > row.product + row.rho:
                rho_failures.append(tag)
            if row.b * max_p > max_b * row.product:
                max_b, max_p = row.b, row.product
                max_rows = [tag]
            elif row.b * max_p == max_b * row.product:
                max_rows.append(tag)
            digest.update((repr(tag) + "\n").encode())

        bc = type_b_parameters(reversals, z)
        if bc is not None:
            b, c = bc
            m = z - b - 1
            h, K, eta, gamma, lam = constants(b, c)
            need(profile.c == gamma, f"cycle formula n={n},b={b},c={c}")
            if m >= 1:
                need(profile.h == eta * 2**m + h, f"H formula n={n},b={b},c={c}")
                need(profile.end[z] == h, f"End_z formula n={n},b={b},c={c}")
                for t in range(1, m + 1):
                    owner = b + t
                    u = 2 ** (t - 1)
                    row = by_pair[(owner, z)]
                    exclusion = lam if t == m else 0
                    expected_end = eta * u
                    expected_n = K * u - exclusion
                    expected_b = (K - eta) * u - exclusion - h + gamma
                    expected_p = eta * h * u
                    expected_rho = (b + t) * (b + m) - (b + t - 1)
                    need(profile.end[owner] == expected_end, "End_tail formula")
                    need(row.n_cover == expected_n, "N_tail,z formula")
                    need(row.b == expected_b, "B_tail,z formula")
                    need(row.product == expected_p, "P_tail,z formula")
                    need(row.minimum == h, "min_tail,z formula")
                    need(row.rho == expected_rho, "rho_tail,z formula")
                    need(25 * expected_b <= 27 * expected_p, "27/25 symbolic tail")
                    formula_rows += 1
                    digest.update(
                        f"formula:{n}:{b}:{c}:{m}:{t}:{h}:{K}:{eta}:{gamma}:{lam}:"
                        f"{expected_end}:{expected_n}:{expected_b}:{expected_p}:{expected_rho}\n".encode()
                    )

    expected_pair_rows = len(classified) * n * (n - 1) // 2
    need(pair_rows == expected_pair_rows, f"pair universe n={n}")
    need(not candidate_failures, f"27/25 failure n={n}")
    need(not five_failures, f"five-copy failure n={n}")
    need(len(five_equalities) == (1 if n == 3 else 0), f"five-copy equality n={n}")

    # Every actual local failure is precisely an inherited X_m tail,z row.
    for rows, threshold in ((product_failures, 3), (plus_min_failures, 4)):
        for reversals, owner_pair, *_ in rows:
            need(reversals == ((0, 3), (3, z)), f"non-X failure n={n}")
            need(owner_pair[1] == z and owner_pair[0] >= 3 + threshold, "wrong X owner")
    for reversals, owner_pair, *_ in rho_failures:
        need(reversals == ((0, 3), (3, z)), f"non-X rho failure n={n}")
        need(owner_pair[1] == z, f"wrong X rho owners n={n}")
    expected_product = n - 7 if n >= 9 else 0
    expected_plus_min = n - 8 if n >= 10 else 0
    expected_rho = n - 11 if n >= 12 else 0
    need(len(product_failures) == expected_product, f"product failure count n={n}")
    need(len(plus_min_failures) == expected_plus_min, f"plus-min count n={n}")
    need(len(rho_failures) == expected_rho, f"rho failure count n={n}")

    if n == 3:
        need((max_b, max_p) == (0, 1), "C3 max ratio")
    elif n <= 8:
        need(max_b == max_p, f"unit max ratio n={n}")
    else:
        u = 2 ** (n - 7)
        need(max_b * (25 * u) == max_p * (27 * u - 4), f"max ratio n={n}")
        need(len(max_rows) == 1, f"max ratio multiplicity n={n}")
        need(max_rows[0][0] == ((0, 3), (3, z)), f"max ratio presentation n={n}")
        need(max_rows[0][1] == (n - 3, z), f"max ratio owners n={n}")

    if n >= 5:
        need(minimum_gap == gap_t_n_1(n), f"T(n,1) minimum gap n={n}")
        need(len(minimum_gap_rows) == (5 if n == 5 else 3), f"minimum multiplicity n={n}")
    if n >= 6:
        expected_minimizers = {
            tuple(sorted(((0, z - 2), (z - 2, z)))),
            tuple(sorted(((0, z - 1), (z - 2, z)))),
            tuple(sorted(((0, z), (z - 2, z)))),
        }
        need(
            {reversals for reversals, _ in minimum_gap_rows} == expected_minimizers,
            f"minimum presentations n={n}",
        )

    first = minimum_gap_rows[0][1]
    minimum_presentations = [reversals for reversals, _ in minimum_gap_rows]
    rho_rows = [(reversals, owner_pair) for reversals, owner_pair, *_ in rho_failures]
    print(
        f"n={n} candidates={len(tuple(combinations(arcs,2)))} strong={len(classified)} "
        f"pairs={pair_rows} formula_rows={formula_rows} "
        f"product_fail={len(product_failures)} plus_min_fail={len(plus_min_failures)} "
        f"rho_fail={len(rho_failures)} 27over25_fail=0 "
        f"Gmin={minimum_gap} Gmin_mult={len(minimum_gap_rows)} "
        f"Gmin_presentations={minimum_presentations} rho_rows={rho_rows} "
        f"Gtarget={first.target} Genergy={first.energy} maxB={max_b} maxP={max_p}"
    )
    grand_presentations += len(classified)
    grand_pair_rows += pair_rows
    grand_formula_rows += formula_rows

# Hostile filter control: the transitive tournament is not strong and has G<0.
transitive = dp.from_reversals(7, ())
transitive_profile = dp.tournament_profile(transitive)
need(not dp.is_strong(transitive) and transitive_profile.gap < 0, "filter hostile control")

print(
    f"totals_presentations={grand_presentations} totals_pairs={grand_pair_rows} "
    f"totals_typeB_tail_formula_rows={grand_formula_rows}"
)
print(f"hostile_transitive_n7_G={transitive_profile.gap}")
print(f"audit_digest={digest.hexdigest()}")
print("all_checks=PASS")
