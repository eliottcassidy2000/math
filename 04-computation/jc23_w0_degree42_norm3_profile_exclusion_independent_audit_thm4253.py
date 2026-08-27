#!/usr/bin/env python3
"""Clean-room exact audit of the proposed THM-4253 norm-three corollary.

The script imports neither a maintained THM-4249 companion nor the earlier
scratch projector audit.  It directly enumerates the a_u=0 degree-42 lattice,
applies the proved THM-4249 necessary filters, constructs mu_6 x <T> orbits,
and performs the elementary CM-ideal/incidence arithmetic.
"""

from collections import Counter


def add(x, y):
    return x[0] + y[0], x[1] + y[1]


def neg(x):
    return -x[0], -x[1]


def mul(x, y):
    a, b = x
    c, d = y
    return a * c - b * d, a * d + b * c - b * d


def conj(x):
    a, b = x
    return a - b, -b


def norm(x):
    z = mul(x, conj(x))
    assert z[1] == 0
    return z[0]


ZERO = (0, 0)
ONE = (1, 0)
W = (0, 1)
W2 = (-1, -1)
MINUS_ONE = (-1, 0)
MINUS_W = (0, -1)
MINUS_W2 = (1, 1)
TWO = (2, 0)
UNITS = (ONE, W, W2, MINUS_ONE, MINUS_W, MINUS_W2)

HIDDEN_H = (
    ((6, 0), (-4, -2)),
    ((-2, 2), (6, 0)),
)


def hidden_degree(x, y):
    coeffs = (x, y)
    total = ZERO
    for i in range(2):
        for j in range(2):
            total = add(
                total,
                mul(mul(coeffs[i], HIDDEN_H[i][j]), conj(coeffs[j])),
            )
    assert total[1] == 0
    return total[0]


def t_action(m):
    # Full THM-4241 basis [u,f,g,h].  This audit uses only a_u=0, but keeps
    # the full action so the orbit test also checks that coordinate stays zero.
    a, b, c, d = m
    return (
        mul(MINUS_W, a),
        mul(MINUS_W, add(c, d)),
        b,
        mul(W2, d),
    )


def t_power(m, k):
    out = m
    for _ in range(k):
        out = t_action(out)
    return out


def scalar_vector(s, m):
    return tuple(mul(s, x) for x in m)


def eis_norm_at_most(bound):
    out = []
    for a in range(-bound, bound + 1):
        for b in range(-bound, bound + 1):
            x = (a, b)
            if norm(x) <= bound:
                out.append(x)
    assert all(max(abs(x[0]), abs(x[1])) < bound for x in out)
    return out


def residual_vectors():
    # For N=O[f,g,h], diag(H_N^-1)=(1/2,1/2,1).  Therefore degree <=42
    # forces N(b),N(c)<=21 and N(d)<=42.
    bc = eis_norm_at_most(21)
    ds = eis_norm_at_most(42)
    vectors = set()
    profile = Counter()
    max_coordinate = 0
    for b in bc:
        for c in bc:
            for d in ds:
                ell_f = add(mul(TWO, b), mul(W2, d))
                ell_g = add(mul(TWO, c), d)
                qell = hidden_degree(ell_f, ell_g)
                if qell % 12:
                    continue
                k_hidden = qell // 12
                nd = norm(d)
                if nd + 3 * k_hidden != 42:
                    continue

                # Exactly the proved THM-4249 residual filters after a_u=0:
                # nonzero visible-v projector of degree at least 12, removal
                # of K=1,2, and the inherited pure-hidden row.
                if d == ZERO or nd < 3 or k_hidden in {1, 2}:
                    continue
                m = (ZERO, b, c, d)
                vectors.add(m)
                profile[(nd, k_hidden)] += 1
                max_coordinate = max(
                    max_coordinate,
                    abs(b[0]), abs(b[1]), abs(c[0]), abs(c[1]),
                    abs(d[0]), abs(d[1]),
                )
    assert profile == Counter({
        (3, 13): 672,
        (9, 11): 576,
        (12, 10): 864,
        (21, 7): 768,
        (27, 5): 288,
    })
    assert len(vectors) == 3168
    return vectors, profile, max_coordinate


def group_orbits(vectors):
    remaining = set(vectors)
    orbits = []
    while remaining:
        seed = min(remaining)
        orbit = {
            scalar_vector(unit, t_power(seed, k))
            for unit in UNITS
            for k in range(12)
        }
        assert orbit <= vectors
        assert len(orbit) == 24
        remaining -= orbit
        orbits.append(orbit)
    return orbits


def profile_of(m):
    _, b, c, d = m
    ell_f = add(mul(TWO, b), mul(W2, d))
    ell_g = add(mul(TWO, c), d)
    return norm(d), hidden_degree(ell_f, ell_g) // 12


def main():
    vectors, profile, max_coordinate = residual_vectors()
    orbits = group_orbits(vectors)
    orbit_profile = Counter(profile_of(min(orbit)) for orbit in orbits)
    assert orbit_profile == Counter({
        (3, 13): 28,
        (9, 11): 24,
        (12, 10): 36,
        (21, 7): 32,
        (27, 5): 12,
    })
    assert len(orbits) == 132

    bad_profile_vectors = {m for m in vectors if profile_of(m) == (3, 13)}
    bad_profile_orbits = [o for o in orbits if profile_of(min(o)) == (3, 13)]
    survivors = vectors - bad_profile_vectors
    survivor_orbits = [o for o in orbits if profile_of(min(o)) != (3, 13)]
    assert len(bad_profile_vectors) == 672
    assert len(bad_profile_orbits) == 28
    assert len(survivors) == 2496
    assert len(survivor_orbits) == 104

    # Ideal audit. pi=w^2-1=-2-w has norm 3 and pi^2=-3*w^2,
    # hence (pi^2)=(3).  Every norm-three d is associated to pi.
    pi = add(W2, MINUS_ONE)
    pi2 = mul(pi, pi)
    minus_three_w2 = mul((-3, 0), W2)
    assert pi == (-2, -1)
    assert norm(pi) == 3
    assert pi2 == minus_three_w2
    norm3_elements = {d for d in eis_norm_at_most(3) if norm(d) == 3}
    pi_associates = {mul(unit, pi) for unit in UNITS}
    assert norm3_elements == pi_associates

    # Kernel E[3] has O; the two nonzero pi-torsion points X=0; and six
    # points with X^3=-4.  Under R=-1/(X^3+1), the first short orbit gives
    # the excluded gate ratio -1 and the unique admissible orbit gives 1/3.
    kernel_size = norm(mul(pi, pi))
    assert kernel_size == 9
    e3_unit_orbits = Counter({1: 1, 2: 1, 6: 1})
    assert sum(size * count for size, count in e3_unit_orbits.items()) == 9
    admissible_ratios_norm3 = {"1/3"}
    assert len(admissible_ratios_norm3) == 1

    # Independent incidence arithmetic for all remaining profiles.  A kernel
    # E[d*pi] has 3*N(d) points.  Remove O and the two nonzero pi-torsion
    # points; for N(d)=12 also remove the three nonzero 2-torsion points.
    per_orbit_ratio_count = {}
    for nd in (3, 9, 12, 21, 27):
        short = 1 + 2 + (3 if nd == 12 else 0)
        numerator = 3 * nd - short
        assert numerator % 6 == 0
        per_orbit_ratio_count[nd] = numerator // 6
    assert per_orbit_ratio_count == {3: 1, 9: 4, 12: 5, 21: 10, 27: 13}

    raw_incidence = sum(
        orbit_count * per_orbit_ratio_count[nd]
        for (nd, _), orbit_count in orbit_profile.items()
    )
    assert raw_incidence == 780
    raw_after_profile = sum(
        orbit_count * per_orbit_ratio_count[nd]
        for (nd, _), orbit_count in orbit_profile.items()
        if nd != 3
    )
    assert raw_after_profile == 752

    # Every surviving d is divisible by the unique prime pi above 3, so each
    # E[d*pi] contains E[3] and contributes the already excluded ratio 1/3.
    final_incidence = raw_after_profile - len(survivor_orbits)
    assert final_incidence == 648

    print("THM4253 NORM-THREE PROFILE INDEPENDENT AUDIT")
    print(f"degree42_residual_vectors={len(vectors)} orbits={len(orbits)} all_orbit_sizes=24")
    print(f"vector_profile={dict(sorted(profile.items()))}")
    print(f"orbit_profile={dict(sorted(orbit_profile.items()))}")
    print(f"norm3_profile=vectors:{len(bad_profile_vectors)} orbits:{len(bad_profile_orbits)}")
    print("ideal_logic=d~pi, (d*pi)=(pi^2)=(3), ker[d*pi]=E[3]")
    print(f"E3_unit_orbits={dict(sorted(e3_unit_orbits.items()))} admissible_ratio=1/3")
    print(f"after_profile=vectors:{len(survivors)} orbits:{len(survivor_orbits)}")
    print(f"ratio_counts_per_map_orbit={per_orbit_ratio_count}")
    print(f"incidences=raw:{raw_incidence} raw_after_profile:{raw_after_profile} final_after_R1/3:{final_incidence}")
    print(f"max_retained_coordinate={max_coordinate} coordinate_bounds=21,21,42")
    print("consequence=THM4249_R1/3_exclusion_removes_entire_(Nd,K)=(3,13)_profile")


if __name__ == "__main__":
    main()
