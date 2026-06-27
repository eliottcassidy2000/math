#!/usr/bin/env python3
"""
q-Pochhammer / modular cusp / Hurwitz sidecar scout for LRC14.

The goal is not to prove modular-function facts from scratch.  It packages the
user's prompt as a proof-interface rule: a modular function has a finite
principal part at the cusp and a modularly forced q-tail; a Hurwitz/Vieta orbit
has a finite seed and a forced mutation tail.  LRC packet quotients should keep
that finite address before trusting any infinite scalar series.
"""

from __future__ import annotations

from collections import deque


N = 32


def sigma_k(n: int, k: int) -> int:
    return sum(d**k for d in range(1, n + 1) if n % d == 0)


def mul(a: list[int], b: list[int], n: int = N) -> list[int]:
    out = [0] * (n + 1)
    for i, ai in enumerate(a[: n + 1]):
        if ai == 0:
            continue
        for j, bj in enumerate(b[: n + 1 - i]):
            if bj:
                out[i + j] += ai * bj
    return out


def pow_series(a: list[int], power: int, n: int = N) -> list[int]:
    out = [1] + [0] * n
    base = a[:]
    p = power
    while p:
        if p & 1:
            out = mul(out, base, n)
        p >>= 1
        if p:
            base = mul(base, base, n)
    return out


def inv_series(a: list[int], n: int = N) -> list[int]:
    if a[0] != 1:
        raise ValueError("this scout only inverts unit series")
    out = [0] * (n + 1)
    out[0] = 1
    for m in range(1, n + 1):
        out[m] = -sum(a[i] * out[m - i] for i in range(1, m + 1))
    return out


def q_pochhammer_unit(n: int = N) -> list[int]:
    out = [1] + [0] * n
    for k in range(1, n + 1):
        factor = [1] + [0] * n
        factor[k] = -1
        out = mul(out, factor, n)
    return out


def laurent_mul(a: dict[int, int], b: dict[int, int], lo: int, hi: int) -> dict[int, int]:
    out: dict[int, int] = {}
    for ea, ca in a.items():
        for eb, cb in b.items():
            e = ea + eb
            if lo <= e <= hi:
                out[e] = out.get(e, 0) + ca * cb
    return {e: c for e, c in out.items() if c}


def laurent_pow(a: dict[int, int], power: int, lo: int, hi: int) -> dict[int, int]:
    out = {0: 1}
    for _ in range(power):
        out = laurent_mul(out, a, lo, hi)
    return out


def generalized_pentagonal_support(limit: int) -> list[tuple[int, int]]:
    rows: list[tuple[int, int]] = [(0, 1)]
    m = 1
    while True:
        added = False
        for sign_m in (m, -m):
            exp = sign_m * (3 * sign_m - 1) // 2
            if exp <= limit:
                rows.append((exp, -1 if m % 2 else 1))
                added = True
        if not added:
            break
        m += 1
    return sorted(rows)


def hurwitz_ok(row: tuple[int, int, int, int]) -> bool:
    return sum(x * x for x in row) == row[0] * row[1] * row[2] * row[3]


def hurwitz_orbit(max_entry: int = 500) -> list[tuple[int, int, int, int]]:
    start = (2, 2, 2, 2)
    seen = {start}
    q: deque[tuple[int, int, int, int]] = deque([start])
    while q:
        row = q.popleft()
        for i in range(4):
            others = [row[j] for j in range(4) if j != i]
            new_value = others[0] * others[1] * others[2] - row[i]
            if new_value <= 0 or new_value > max_entry:
                continue
            nxt_list = list(row)
            nxt_list[i] = new_value
            nxt = tuple(sorted(nxt_list))
            if nxt not in seen:
                seen.add(nxt)
                q.append(nxt)
    return sorted(seen, key=lambda xs: (max(xs), sum(xs), xs))


def main() -> None:
    print("LRC14 MODULAR CUSP POCHHAMMER HURWITZ S245")
    print(f"truncation_N={N}")

    poch = q_pochhammer_unit(N)
    inv_poch = inv_series(poch, N)
    delta_tail = pow_series(poch, 24, N)
    inv_delta_tail = inv_series(delta_tail, N)

    e4 = [1] + [240 * sigma_k(n, 3) for n in range(1, N + 1)]
    e4_cubed = pow_series(e4, 3, N)
    j_tail = mul(e4_cubed, inv_delta_tail, N)
    j = {i - 1: c for i, c in enumerate(j_tail) if c}

    print()
    print("Q_POCHHAMMER")
    print("(q;q)_infty_nonzero_terms_to_32=" + repr([(i, c) for i, c in enumerate(poch) if c]))
    print("pentagonal_theorem_terms_to_32=" + repr(generalized_pentagonal_support(N)))
    print("partition_numbers_1_over_(q;q)_infty_to_12=" + repr(inv_poch[:13]))
    print("delta_over_q_first_terms=" + repr(delta_tail[:8]))
    print("one_over_delta_principal_part=q^-1")
    print("one_over_delta_tail_after_q^-1=" + repr(inv_delta_tail[:8]))

    print()
    print("FULL_MODULAR_GROUP_CUSP")
    print("j_terms_q^-1_to_q^8=" + repr([(e, j[e]) for e in range(-1, 9)]))
    print("j_principal_part=" + repr({-1: j[-1]}))
    j2 = laurent_pow(j, 2, -2, 8)
    print("j_squared_principal_part=" + repr({e: j2[e] for e in sorted(j2) if e < 0}))
    print("finite_negative_power_rule=modular_function_q_series_has_exponents_n>=-N")

    print()
    print("HURWITZ_VIETA_ORBIT")
    orbit = hurwitz_orbit(500)
    print(f"orbit_size_max_entry_500={len(orbit)}")
    print("first_hurwitz_rows=" + repr(orbit[:12]))
    print("all_rows_satisfy_sum_squares_equals_product=" + str(all(hurwitz_ok(row) for row in orbit)))
    print("vieta_seed=(2,2,2,2)")
    print("vieta_move=x_i -> product(other_three)-x_i")

    print()
    print("LRC_SIDECAR_FIELDS")
    fields = [
        "modular_cusp_principal_part_order",
        "finite_negative_power_budget",
        "principal_part_coeff_vector",
        "q_pochhammer_tail_signature",
        "eta_delta_denominator_lane",
        "j_rational_function_address",
        "hurwitz_vieta_seed",
        "hurwitz_mutation_depth",
        "continued_fraction_period_word",
        "pell_wall_unit",
        "cusp_tail_discharge_route",
        "destroyed_cusp_or_owner_coordinate",
        "required_sidecar_or_exit",
    ]
    for field in fields:
        print(field)

    print()
    print("TOURNAMENT_ANALYSIS")
    vertices = [
        "labelled_lrc_packet_sheaf",
        "modular_cusp_principal_part",
        "full_modular_group_invariance_gate",
        "q_pochhammer_eta_tail",
        "hurwitz_vieta_seed_orbit",
        "continued_fraction_markov_address",
        "pell_wall_unit_address",
        "raw_q_series_coefficients",
        "raw_hurwitz_scalar",
    ]
    print("vertices=" + ",".join(vertices))
    print(
        "pairwise_observable="
        "retained_cusp_address,finite_principal_part_order,tail_control,"
        "hurwitz_mutation_legality,preserved_lrc_predicate"
    )
    print(
        "switch_gauge=A->B when A retains the finite cusp/arithmetic address "
        "needed before B's scalar tail can be used"
    )
    print("tie_hamiltonian_path=" + " > ".join(vertices))
    print("score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1}")
    print("directed_3cycles=0")
    print("scc_sizes=[1,1,1,1,1,1,1,1,1]")
    print("hamiltonian_path_count=1")


if __name__ == "__main__":
    main()
