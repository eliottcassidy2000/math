#!/usr/bin/env python3
"""THM-3200 exact fixed-channel quasipolynomial audit for reflected cell-90 overlaps.

Run from the repository root.  This imports only the canonical exact
digital-tent/stability engine, then independently reconstructs the eventual
quadratic polynomials, their least residue period, the step-period recurrence,
and the moving-denominator obstruction for four hostile/positive rays.
"""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = next(parent for parent in HERE.parents
            if (parent / "00-navigation").is_dir())
ENGINE = ROOT / "04-computation/lrc_small_channel_symbolic_tail_thm3171.py"
EXPECTED_ENGINE = "d73273a4cf4b88bea2890e001166d96cb07dd9b61f3a248ff1538ec44579796a"
AUDIT_MAX = 1200
CASES = (
    ("sharp_6_to_1", 3, 5, 6, 1),
    ("reverse_1_to_6", 3, 5, 1, 6),
    ("label12_12_to_1", 3, 5, 12, 1),
    ("label12_1_to_12", 3, 5, 1, 12),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def ceildiv(a, b):
    return -((-a) // b)


def divisors(n):
    return tuple(d for d in range(1, n + 1) if n % d == 0)


def polynomial_in_g(values, base, step):
    """Interpolate a quadratic from values at base+step*k, k=0,1,2."""
    v0, v1, v2 = map(F, values[:3])
    a = (v2 - 2 * v1 + v0) / 2
    b = v1 - v0 - a
    c = v0
    # Convert a*k^2+b*k+c, k=(g-base)/step, to A*g^2+B*g+C.
    A = a / (step * step)
    B = b / step - 2 * a * base / (step * step)
    C = a * base * base / (step * step) - b * base / step + c
    return A, B, C


def value(poly, g):
    a, b, c = poly
    return a * g * g + b * g + c


def collapse_period(rows, ambient):
    """Find the least divisor whose congruence classes have equal polynomials."""
    for period in divisors(ambient):
        good = True
        for residue, poly, _ in rows:
            for other, other_poly, _ in rows:
                if (residue - other) % period == 0 and poly != other_poly:
                    good = False
                    break
            if not good:
                break
        if good:
            collapsed = {}
            for residue, poly, _ in rows:
                key = residue % period
                require(key not in collapsed or collapsed[key] == poly,
                        ("period collapse", period, key))
                collapsed[key] = poly
            require(len(collapsed) == period, ("missing residues", period, collapsed))
            return period, collapsed
    raise RuntimeError(("no period", ambient))


def audit_case(engine, case):
    name, P, Q, e, f = case
    cross = Q * e - P * f
    ambient = abs(cross) or 1
    g_min = max(1, ceildiv(6, P))
    rows = []
    stable_max = 0

    # The engine's affine-branch certificate proves that each ambient residue
    # stays on one quadratic branch after the returned point.
    for residue in range(1, ambient + 1):
        h_min = max(0, ceildiv(g_min - residue, ambient))
        h = max(1, h_min)
        while not engine.ray_is_stable(P, Q, e, f, residue, ambient, h):
            h *= 2
            require(h <= 1 << 20, ("stability ceiling", case, residue))
        base = residue + ambient * h
        stable_max = max(stable_max, base)
        samples = tuple(
            engine.fast_mass(e, (base + ambient * k) * P,
                             f, (base + ambient * k) * Q)[0]
            for k in range(4)
        )
        poly = polynomial_in_g(samples, base, ambient)
        require(value(poly, base + 3 * ambient) == samples[3],
                ("fourth point", case, residue, samples, poly))
        rows.append((residue, poly, base))

    exact_period, polynomials = collapse_period(rows, ambient)

    # Exact head checks bridge the physical start to every certified tail.
    head_checks = 0
    for residue, poly, base in rows:
        g = residue
        while g < g_min:
            g += ambient
        while g < base:
            numerator = engine.fast_mass(e, g * P, f, g * Q)[0]
            require(F(numerator) == value(poly, g),
                    ("head formula", case, residue, g, numerator, poly))
            head_checks += 1
            g += ambient

    # A long direct hostile window checks the collapsed formula and
    # (E^period-1)^3 annihilator independently of the branch bookkeeping.
    numerators = {}
    formula_checks = 0
    for g in range(g_min, AUDIT_MAX + 1):
        numerator = engine.fast_mass(e, g * P, f, g * Q)[0]
        require(F(numerator) == value(polynomials[g % exact_period], g),
                ("collapsed formula", case, g, numerator,
                 polynomials[g % exact_period]))
        numerators[g] = numerator
        formula_checks += 1

    recurrence_checks = 0
    for g in range(g_min, AUDIT_MAX - 3 * exact_period + 1):
        third = (
            numerators[g + 3 * exact_period]
            - 3 * numerators[g + 2 * exact_period]
            + 3 * numerators[g + exact_period]
            - numerators[g]
        )
        require(third == 0, ("third step difference", case, g, third))
        recurrence_checks += 1

    # The exact mass is poly(g)/D(g) on each residue.  A nonzero 1/g
    # coefficient proves polynomial, rather than exponential, convergence
    # and therefore excludes an eventual constant-coefficient recurrence.
    alpha, beta = 168 * P, 168 * Q
    den_a = alpha * beta
    den_b = -(alpha * f + beta * e)
    den_c = e * f
    limits = set()
    corrections = set()
    eventually_constant = True
    for poly in polynomials.values():
        a, b, c = poly
        limit = a / den_a
        correction = (b * den_a - a * den_b) / (den_a * den_a)
        limits.add(limit)
        corrections.add(correction)
        if (a, b, c) != (limit * den_a, limit * den_b, limit * den_c):
            eventually_constant = False
    require(len(limits) == 1, ("limit mismatch", case, limits))
    require(not eventually_constant, ("unexpected constant mass", case))
    require(all(correction != 0 for correction in corrections),
            ("vanishing first correction", case, corrections))

    ordered = tuple(sorted(polynomials.items()))
    constants = tuple(poly[2] for _, poly in ordered)
    digest = sha256(repr(ordered).encode()).hexdigest()
    return {
        "name": name,
        "P": P, "Q": Q, "e": e, "f": f,
        "cross": cross,
        "ambient": ambient,
        "period": exact_period,
        "stable_max": stable_max,
        "head_checks": head_checks,
        "formula_checks": formula_checks,
        "recurrence_checks": recurrence_checks,
        "leading": next(iter(polynomials.values()))[0],
        "linear": next(iter(polynomials.values()))[1],
        "constant_count": len(set(constants)),
        "constant_min": min(constants),
        "constant_max": max(constants),
        "limit": next(iter(limits)),
        "corrections": tuple(sorted(corrections)),
        "poly_digest": digest,
        "polynomials": ordered,
    }


def main():
    actual = sha256(ENGINE.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(actual == EXPECTED_ENGINE, ("engine hash", actual, EXPECTED_ENGINE))
    spec = spec_from_file_location("fixed_channel_engine", ENGINE)
    engine = module_from_spec(spec)
    spec.loader.exec_module(engine)

    rows = tuple(audit_case(engine, case) for case in CASES)
    print("FIXED-CHANNEL CELL-90 DILATION QUASIPOLYNOMIAL AUDIT -- THM-3200")
    print(f"engine_lf_sha256={actual};audit_max={AUDIT_MAX}")
    for row in rows:
        print(
            f"case={row['name']};lane=({row['e']},{row['f']});"
            f"primitive=({row['P']},{row['Q']});cross={row['cross']};"
            f"period_bound={row['ambient']};exact_period={row['period']};degree=2;"
            f"stable_g_max={row['stable_max']};head_checks={row['head_checks']};"
            f"formula_checks={row['formula_checks']};"
            f"recurrence_checks={row['recurrence_checks']}"
        )
        print(
            f"coefficients=leading:{row['leading']},linear:{row['linear']},"
            f"constant_classes:{row['constant_count']},"
            f"constant_range:{row['constant_min']}..{row['constant_max']};"
            f"poly_sha256={row['poly_digest']}"
        )
        print(
            f"mass_limit={row['limit']};one_over_g={row['corrections']};"
            "mass_eventually_constant=NO;mass_rational_ogf=NO"
        )
        if row["period"] <= 3:
            print(f"exact_polynomials={row['polynomials']}")
    semantic = repr(rows).encode()
    print(f"semantic_sha256={sha256(semantic).hexdigest()}")


if __name__ == "__main__":
    main()
