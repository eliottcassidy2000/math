#!/usr/bin/env python3
"""Exact controls for THM-2478.

Dependency-free Fraction arithmetic only.  The finite grids are exact
integrals of rational step functions: every breakpoint is a grid endpoint
and every midpoint lies in one constant cell.
"""

from fractions import Fraction as F


P = 13
D = 5


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def cell(values, multiplier, j, cells):
    """Value of a D-cell function at multiplier*(j+1/2)/cells."""
    index = (D * multiplier * (2 * j + 1) // (2 * cells)) % D
    return F(values[index])


def exact_integral_product(factors, cells):
    return sum(
        prod(cell(values, multiplier, j, cells) for values, multiplier in factors)
        for j in range(cells)
    ) / cells


def prod(values):
    answer = F(1)
    for value in values:
        answer *= value
    return answer


def cyclic_variation(values):
    return sum(abs(F(values[(j + 1) % len(values)]) - F(values[j])) for j in range(len(values)))


def direct_handoff_control():
    # E(z) and Q(13^K z) form one positive owner-to-word block G(z).
    owner = [1, 1, 0, 0, 0]
    word = [0, 1, 1, 0, 1]
    drift = [2, -1, 3, 0, -2]
    service = [1, 0, 2, 1, 3]
    clock = 1

    handoff_cells = D * P**clock
    rho = exact_integral_product([(owner, 1), (word, P**clock)], handoff_cells)
    q = sum(map(F, drift)) / D
    mass = sum(map(F, service)) / D
    require(rho > 0 and q != 0 and mass > 0, "positive direct control")

    gate_values = [
        int(cell(owner, 1, j, handoff_cells) * cell(word, P**clock, j, handoff_cells))
        for j in range(handoff_cells)
    ]
    variation_gate = cyclic_variation(gate_values)
    variation_drift = cyclic_variation(drift)
    variation_service = cyclic_variation(service)
    constant = min(rho * (1 - rho), F(variation_gate, 12))

    rows = []
    for delay in range(1, 5):
        scale = P ** (delay - 1)
        cells = D * P ** (delay - 1 + clock)
        q_delay = exact_integral_product(
            [(owner, scale), (word, scale * P**clock), (drift, 1)], cells
        )
        mass_delay = exact_integral_product(
            [(owner, scale), (word, scale * P**clock), (service, 1)], cells
        )
        require(abs(q_delay - rho * q) <= constant * variation_drift / scale,
                "drift BV invoice")
        require(abs(mass_delay - rho * mass) <= constant * variation_service / scale,
                "service BV invoice")
        require(q_delay != 0 and mass_delay > 0, "simultaneous retention")
        rows.append((delay, q_delay, mass_delay))
    return rho, q, mass, variation_gate, constant, rows


def add_poly(a, b):
    return tuple(x + y for x, y in zip(a, b))


def scale_poly(c, a):
    return tuple(c * x for x in a)


def root_power(exponent):
    exponent %= P
    if exponent == P - 1:
        return tuple(F(-1) for _ in range(P - 1))
    answer = [F(0) for _ in range(P - 1)]
    answer[exponent] = F(1)
    return tuple(answer)


ZERO = tuple(F(0) for _ in range(P - 1))


def format_poly(poly):
    return "[" + ",".join(str(value) for value in poly) + "]"


def collision_colour_control():
    # On five base cells put one A-root and one distinct F-root.  The
    # correlation shifts are 1,...,5, so C(0)=0.  Each g_k is the DFT of
    # nonnegative Boolean shift densities before the transform.
    shifts = [1, 2, 3, 4, 5]
    drift = [2, -1, 3, 0, -2]
    q = sum(map(F, drift)) / D
    divisor = P  # a minimal owner-normalization multiplier for the control

    means = {}
    delayed = {}
    for k in range(1, P):
        values = [scale_poly(F(1, P * P), root_power(-k * shift)) for shift in shifts]
        mean = ZERO
        for value in values:
            mean = add_poly(mean, scale_poly(F(1, D), value))
        require(mean != ZERO, "every primitive collision colour survives")
        means[k] = mean

        # Exact old-drift/future-colour fibre product at three delays.
        delayed[k] = []
        errors = []
        for delay in range(3):
            multiplier = divisor * P**delay
            cells = D * multiplier
            total = ZERO
            for j in range(cells):
                old_value = cell(drift, 1, j, cells)
                future_cell = (D * multiplier * (2 * j + 1) // (2 * cells)) % D
                total = add_poly(total, scale_poly(old_value / cells, values[future_cell]))
            require(total != ZERO, "old drift survives in each future colour")
            delayed[k].append(total)

        limit = scale_poly(q, mean)
        for total in delayed[k]:
            errors.append(sum(abs(x - y) for x, y in zip(total, limit)))
        require(errors[-1] < errors[0], "future-colour convergence control")

    # C(0)=0 and one Boolean A/F pair contributes on each base cell.
    require(0 not in shifts and len(shifts) == D, "diagonal collision slice is empty")
    return len(means), means[1], delayed[1][-1]


def geometry_controls():
    # Collision-base placement: y_col=13*d*13^L*x becomes d*13^L*y on
    # x=(y+u)/13, independently of the old root u.
    y = F(3, 17)
    divisor = P
    delay = 2
    values = {
        ((P * divisor * P**delay * ((y + u) / P)) % 1)
        for u in range(P)
    }
    require(values == {(divisor * P**delay * y) % 1}, "old-root constancy")

    # Every target co-shift theta/13 is killed by the collision multiplier.
    x = F(5, 29)
    for theta in range(P):
        shifted = (P * divisor * P**delay * (x - F(theta, P))) % 1
        unshifted = (P * divisor * P**delay * x) % 1
        require(shifted == unshifted, "target neutrality")

    # Full THM-2471 stalk: with terminal clock K=1 and collision delay L=2,
    # the X,Y legs have multiplier 13^L and the ancestry Z leg has the
    # smallest multiplier 13^(L-K).  Every atomic leg is target-neutral.
    clock = 1
    for multiplier in (P**delay, P ** (delay - clock), P**delay):
        for theta in range(P):
            require(F(multiplier * theta, P).denominator == 1,
                    "full stalk target neutrality")

    # Rebase boundary for C=2*13^5.  At L=4 the deep probe is sheet-free;
    # at L=6 the sheets a=0 and a=7 have phases 0 and 1/13.
    deep = 2 * P**5
    z = F(4, 31)
    shallow_delay = 4
    phases = {
        (F(deep * (z + a), P**shallow_delay) % 1)
        for a in range(P**shallow_delay)
    }
    require(len(phases) == 1, "deep phase descends before its valuation")
    require(next(iter(phases)) == (F(deep, P**shallow_delay) * z) % 1,
            "descended deep coefficient")

    lost_delay = 6
    a0 = 0
    a1 = 7  # 2*7=1 mod 13
    phase0 = F(deep * a0, P**lost_delay) % 1
    phase1 = F(deep * a1, P**lost_delay) % 1
    require(phase0 == 0 and phase1 == F(1, P), "essential ancestry residue")
    require(phase0 < F(1, 14) and phase1 > F(1, 14), "danger/safe sheet hostile")
    return next(iter(phases)), phase0, phase1, P ** (lost_delay - 5)


def h_drift_interface_hostile():
    # A circulant owner tensor plus one off-slice nonowner atom.  It obeys
    # nonnegativity and H(t,s,t)=0, but owner drift is zero while aggregate
    # drift is positive and the untwisted nonowner transition mass is zero.
    owner = {}
    nonowner = {}
    for r in range(P):
        for s in range(P):
            for t in range(P):
                owner[r, s, t] = F(1, 4) if (r - t) % P == 1 else F(0)
                nonowner[r, s, t] = F(1, 4) if (r, s, t) == (1, 1, 0) else F(0)

    def projection(table):
        averages = {
            difference: sum(
                table[r, s, t]
                for r in range(P) for s in range(P) for t in range(P)
                if (r - t) % P == difference
            ) / (P * P)
            for difference in range(P)
        }
        return {(r, s, t): averages[(r - t) % P]
                for r in range(P) for s in range(P) for t in range(P)}

    def drift(table):
        projected = projection(table)
        return sum((table[key] - projected[key]) ** 2 for key in table) / P**3

    owner_drift = drift(owner)
    aggregate_drift = drift({key: owner[key] + nonowner[key] for key in owner})
    require(owner_drift == 0 and aggregate_drift == F(21, 742586), "sharp tensor hostile")
    require(sum(nonowner[r, 0, 0] for r in range(P)) == 0, "no untwisted nonowner edge")
    require(all(owner[t, s, t] + nonowner[t, s, t] == 0
                for s in range(P) for t in range(P)), "diagonal zero")
    return owner_drift, aggregate_drift


def main():
    rho, q, mass, variation_gate, constant, rows = direct_handoff_control()
    colours, first_mean, first_landing = collision_colour_control()
    descended, phase0, phase1, lost_modulus = geometry_controls()
    owner_drift, aggregate_drift = h_drift_interface_hostile()

    print("THM-2478 exact companion: PASS")
    print(f"direct_handoff_mean={rho}; old_drift={q}; old_service={mass}")
    print(f"direct_gate_variation={variation_gate}; BV_constant={constant}")
    print("direct_retention=" + ",".join(f"L{L}:q={qd}:M={md}" for L, qd, md in rows))
    print(f"future_collision_colours={colours}; first_mean={format_poly(first_mean)}")
    print(f"first_colour_old_drift_landing={format_poly(first_landing)}")
    print(f"deep_rebase_descended_phase={descended}; lost_phases=({phase0},{phase1}); lost_modulus={lost_modulus}")
    print(f"H_interface_owner_drift={owner_drift}; aggregate_drift={aggregate_drift}")
    print("old root/deep table retained; future owner/collision roots remain typed at a different scale")


if __name__ == "__main__":
    main()
