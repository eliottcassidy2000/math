#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2534.

The all-slope Boolean boundary transform on ``C_13`` is

    K_tau(r) = q_r (1-q_(r+tau)).

This script exhausts all 8,190 nonconstant Boolean masks and all twelve
nonzero slopes.  It checks the matched cyclic-Hilbert/Cayley flux, the
Crofton and positive-reconstruction identities, the disjoint translated
tail/head packets, sharp Boolean coercivity, and primitive-root colour
saturation.  It then audits the lossless linear equivalence between the
full directed boundary bank and a target-anchored Boolean Gram matrix,
including the deep-path gamma/beta specializations and the exact-gradient
(zero-curl) obstruction to reading the raw boundary imbalance as a cyclic
tournament.

Only integer arithmetic is used.  Primitive thirteenth-root evaluations
are represented exactly in the basis ``1,zeta,...,zeta^11``, using
``zeta^12=-(1+...+zeta^11)``.  All checks remain active under ``python -O``.
"""

from collections import Counter


P = 13
ALL = (1 << P) - 1
checks = 0


def require(condition, message):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError("FAILED: " + message)


def bit(mask, residue):
    return (mask >> (residue % P)) & 1


def vector(mask):
    return tuple(bit(mask, residue) for residue in range(P))


def shift(values, amount):
    """Return P_amount values, with (P_amount f)(r)=f(r+amount)."""

    return tuple(values[(residue + amount) % P] for residue in range(P))


def tail(values, tau):
    return tuple(
        values[residue] * (1 - values[(residue + tau) % P])
        for residue in range(P)
    )


def head(values, tau):
    return tuple(
        values[(residue - tau) % P] * (1 - values[residue])
        for residue in range(P)
    )


def hilbert(values, tau):
    """H_tau=sum_(s=1)^6(P_(-2tau s)-P_(2tau s))."""

    return tuple(
        sum(
            values[(residue - 2 * tau * s) % P]
            - values[(residue + 2 * tau * s) % P]
            for s in range(1, 7)
        )
        for residue in range(P)
    )


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def mask_of(values):
    return sum(value << residue for residue, value in enumerate(values))


def cyclotomic_value(values, colour):
    """Exact sum_r values[r] zeta^(-colour*r) in Q(zeta_13).

    The return value is its coefficient vector in the integral basis
    1,zeta,...,zeta^11.  It is zero exactly when the returned tuple is zero.
    """

    coefficients = [0] * P
    for residue, value in enumerate(values):
        coefficients[(-colour * residue) % P] += value
    top = coefficients[P - 1]
    return tuple(coefficients[index] - top for index in range(P - 1))


def score(mask, anchor, tau):
    return sum(bit(mask, anchor + j * tau) << (P - 1 - j) for j in range(P))


def selector_anchor(mask, tau):
    scores = tuple(score(mask, anchor, tau) for anchor in range(P))
    maximum = max(scores)
    anchors = tuple(anchor for anchor, value in enumerate(scores) if value == maximum)
    require(len(anchors) == 1, "unique prime-necklace anchor")
    return anchors[0]


def inverse_mod(value):
    return pow(value % P, -1, P)


def epsilon(tau):
    representative = inverse_mod(tau)
    return 1 if representative <= 6 else -1


def deep_path_masks():
    masks = [1 << j for j in range(1, P)]
    masks += [(1 << j) | (1 << (j + 1)) for j in range(1, P - 1)]
    return masks


def audit_all_boolean_masks():
    mask_slope_cases = 0
    tail_head_coordinate_cases = 0
    cayley_coordinate_cases = 0
    root_mode_cases = 0
    reconstruction_coordinates = 0
    gradient_point_cases = 0
    sharp_cases = 0
    flux_histogram = Counter()
    distinct_tail_masks = set()

    # The boundary bank alone cannot distinguish the two constant fibres;
    # their different E_0 values are exactly the sidecar used by the stated
    # reconstruction theorem.
    empty = (0,) * P
    full = (1,) * P
    require(all(not any(tail(empty, tau)) for tau in range(1, P)), "empty constant boundary kernel")
    require(all(not any(tail(full, tau)) for tau in range(1, P)), "full constant boundary kernel")

    for mask in range(1, ALL):
        values = vector(mask)
        mass = mask.bit_count()
        require(1 <= mass <= P - 1, "nonconstant mask domain")
        outgoing_degree = [0] * P
        crofton = 0

        for tau in range(1, P):
            tails = tail(values, tau)
            heads = head(values, tau)
            flux = sum(tails)
            translated = shift(tails, -tau)

            require(all(value in (0, 1) for value in tails), "Boolean tail packet")
            require(all(value in (0, 1) for value in heads), "Boolean head packet")
            require(heads == translated, "head is translated tail")
            require(all(tails[r] == heads[(r + tau) % P] for r in range(P)), "matched edge endpoints")
            require(all(tails[r] * heads[r] == 0 for r in range(P)), "tail/head disjoint at one root")
            require(1 <= flux <= 6, "every prime direction crosses a mixed mask")
            require(sum(heads) == flux, "tail/head co-support on the fibre")
            require(0 < mask_of(tails) < ALL, "tail mask is nonconstant")

            image = hilbert(values, tau)
            shifted_values = shift(values, tau)
            require(
                tuple(image[r] + image[(r + tau) % P] for r in range(P))
                == tuple(shifted_values[r] - values[r] for r in range(P)),
                "Cayley identity (I+P_tau)H_tau=P_tau-I",
            )
            require(dot(image, values) == 0, "Hilbert skew self-pair")
            require(dot(image, shifted_values) == flux, "matched Cayley flux")
            require(
                sum((shifted_values[r] - values[r]) ** 2 for r in range(P))
                == 2 * flux,
                "boundary is half the symmetric-difference energy",
            )

            require(42 * flux >= mass * (P - mass), "sharp 13/42 Boolean coercivity")
            if 42 * flux == mass * (P - mass):
                require(mass in (6, 7) and flux == 1, "sharp consecutive-block boundary")
                sharp_cases += 1

            for r, value in enumerate(tails):
                outgoing_degree[r] += value
            crofton += flux
            flux_histogram[flux] += 1
            distinct_tail_masks.add(mask_of(tails))
            mask_slope_cases += 1
            tail_head_coordinate_cases += P
            cayley_coordinate_cases += P

        require(crofton == mass * (P - mass), "Crofton complete-cut count")
        require(
            tuple(outgoing_degree)
            == tuple((P - mass) * values[r] for r in range(P)),
            "pointwise outgoing-degree formula",
        )
        reconstructed = tuple(outgoing_degree[r] // (P - mass) for r in range(P))
        require(reconstructed == values, "positive reconstruction of the Boolean mask")
        require(
            all(outgoing_degree[r] % (P - mass) == 0 for r in range(P)),
            "integral positive reconstruction",
        )
        reconstruction_coordinates += P

        # Pointwise the directed imbalance is the exact potential gradient
        # e_r-e_s.  This is the local source of the integrated zero-curl law.
        for r in range(P):
            for s in range(r + 1, P):
                forward = values[r] * (1 - values[s])
                backward = values[s] * (1 - values[r])
                require(forward - backward == values[r] - values[s], "point cut-gradient identity")
                gradient_point_cases += 1

    # Irreducibility of Phi_13 says a degree-at-most-twelve rational mask
    # has a zero primitive value iff all thirteen coefficients are equal.
    # We also evaluate every distinct tail mask exactly in Q(zeta_13).
    for tail_mask in sorted(distinct_tail_masks):
        values = vector(tail_mask)
        require(0 < tail_mask < ALL, "primitive-colour nonconstant input")
        for colour in range(1, P):
            require(
                any(cyclotomic_value(values, colour)),
                "every primitive root colour survives",
            )
            root_mode_cases += 1

    require(mask_slope_cases == (ALL - 1) * (P - 1), "all mask/slope cases")
    require(root_mode_cases == len(distinct_tail_masks) * (P - 1), "all distinct tail modes")
    require(sharp_cases == 2 * P * (P - 1), "all sharp 6/7 consecutive blocks")
    require(sum(flux_histogram.values()) == mask_slope_cases, "flux histogram total")
    return {
        "mask_slope_cases": mask_slope_cases,
        "tail_head_coordinate_cases": tail_head_coordinate_cases,
        "cayley_coordinate_cases": cayley_coordinate_cases,
        "root_mode_cases": root_mode_cases,
        "root_mode_incidence_cases": mask_slope_cases * (P - 1),
        "reconstruction_coordinates": reconstruction_coordinates,
        "gradient_point_cases": gradient_point_cases,
        "sharp_cases": sharp_cases,
        "flux_histogram": tuple(sorted(flux_histogram.items())),
        "distinct_tail_masks": len(distinct_tail_masks),
    }


def audit_anchored_gram_boundary_equivalence():
    """Audit B(r,s)=K(r,r)-K(r,s) and the anchored inverse."""

    gram = [[0] * P for _ in range(P)]
    boundary = [[0] * P for _ in range(P)]
    anchored_masks = 0
    point_entries = 0

    for mask in range(1, ALL):
        if bit(mask, 0):
            continue
        values = vector(mask)
        weight = 1 + (37 * mask + 11 * mask.bit_count() + mask * mask) % 1009
        anchored_masks += 1
        for r in range(P):
            for s in range(P):
                point_gram = values[r] * values[s]
                point_boundary = values[r] * (1 - values[s]) if r != s else 0
                if r != s:
                    require(
                        point_boundary == values[r] - point_gram,
                        "pointwise boundary/Gram polarization",
                    )
                    point_entries += 1
                gram[r][s] += weight * point_gram
                boundary[r][s] += weight * point_boundary

    require(anchored_masks == (1 << (P - 1)) - 1, "all nonempty anchored masks")
    require(all(gram[0][s] == gram[s][0] == 0 for s in range(P)), "anchored Gram row and column")

    inverse_entries = 0
    metric_entries = 0
    gradient_entries = 0
    for r in range(P):
        for s in range(P):
            if r == s:
                continue
            require(boundary[r][s] == gram[r][r] - gram[r][s], "integrated boundary from Gram")
            if r != 0:
                require(gram[r][r] == boundary[r][0], "anchor recovers Gram diagonal")
                require(
                    gram[r][s] == boundary[r][0] - boundary[r][s],
                    "full anchored Gram inverse",
                )
                inverse_entries += 1
            require(
                boundary[r][s] + boundary[s][r]
                == gram[r][r] + gram[s][s] - 2 * gram[r][s],
                "symmetric boundary is squared Gram distance",
            )
            metric_entries += 1
            require(
                boundary[r][s] - boundary[s][r] == gram[r][r] - gram[s][s],
                "skew boundary is an exact gradient",
            )
            gradient_entries += 1

    zero_curl_cases = 0
    for r in range(P):
        for s in range(P):
            if s == r:
                continue
            for t in range(P):
                if t in (r, s):
                    continue
                circulation = (
                    boundary[r][s] - boundary[s][r]
                    + boundary[s][t] - boundary[t][s]
                    + boundary[t][r] - boundary[r][t]
                )
                require(circulation == 0, "directed boundary has zero triangle curl")
                zero_curl_cases += 1

    return {
        "anchored_masks": anchored_masks,
        "point_entries": point_entries,
        "inverse_entries": inverse_entries,
        "metric_entries": metric_entries,
        "gradient_entries": gradient_entries,
        "zero_curl_cases": zero_curl_cases,
        "gram_checksum": sum((r + 1) * (s + 1) * gram[r][s] for r in range(P) for s in range(P)),
        "boundary_checksum": sum((r + 1) * (s + 1) * boundary[r][s] for r in range(P) for s in range(P)),
    }


def audit_deep_path_specialization():
    masks = deep_path_masks()
    require(len(masks) == 23, "deep path has twelve singleton and eleven pair rays")

    alpha_weights = [j * j + 1 for j in range(1, P)]
    beta_weights = [2 * j + 3 for j in range(1, P - 1)]
    weights = {
        1 << j: alpha_weights[j - 1]
        for j in range(1, P)
    }
    weights.update({
        (1 << j) | (1 << (j + 1)): beta_weights[j - 1]
        for j in range(1, P - 1)
    })

    gram = [[0] * P for _ in range(P)]
    boundary = [[0] * P for _ in range(P)]
    gamma_mass = {
        (tau, a): 0
        for tau in range(1, P)
        for a in range(P)
    }
    gamma_cases = 0
    beta_point_cases = 0
    natural_wall_cases = 0

    for mask in masks:
        values = vector(mask)
        weight = weights[mask]
        require(values[0] == 0 and sum(values) in (1, 2), "anchored one-run mask")

        natural_tails = tail(values, 1)
        require(sum(natural_tails) == 1, "positive deep direction has one terminal wall")
        natural_wall_cases += 1

        for tau in range(1, P):
            anchor = selector_anchor(mask, tau)
            eps = epsilon(tau)
            for a in range(P):
                marker = int(anchor == a)
                quadratic = values[a] * (1 - values[(a - eps) % P])
                kakeya_slice = tail(values, -eps)[a]
                require(marker == quadratic == kakeya_slice, "THM-2531 gamma is a Kakeya slice")
                gamma_mass[(tau, a)] += weight * marker
                gamma_cases += 1

        for j in range(1, P - 1):
            pair = values[j] * values[j + 1]
            require(pair == values[j] - natural_tails[j], "pointwise path beta identity")
            beta_point_cases += 1

        for r in range(P):
            for s in range(P):
                gram[r][s] += weight * values[r] * values[s]
                if r != s:
                    boundary[r][s] += weight * values[r] * (1 - values[s])

    gamma_integrated_cases = 0
    for tau in range(1, P):
        eps = epsilon(tau)
        require(sum(gamma_mass[(tau, a)] for a in range(P)) == sum(weights.values()), "gamma partitions deep mass")
        for a in range(P):
            require(
                gamma_mass[(tau, a)]
                == gram[a][a] - gram[(a - eps) % P][a]
                == boundary[a][(a - eps) % P],
                "integrated THM-2531 gamma/boundary identity",
            )
            gamma_integrated_cases += 1

    beta_integrated_cases = 0
    for j in range(1, P - 1):
        require(
            gram[j][j + 1]
            == boundary[j][0] - boundary[j][j + 1],
            "integrated path beta from target and natural needles",
        )
        require(gram[j][j + 1] == beta_weights[j - 1], "path beta coordinate")
        beta_integrated_cases += 1

    # The natural direction sees one terminal wall on both ray types.  The
    # full Crofton count distinguishes singleton (12) from pair (22), so it
    # recovers THM-2529's two-mass consumer exactly.
    alpha_total = sum(alpha_weights)
    beta_total = sum(beta_weights)
    natural_mass = sum(boundary[r][(r + 1) % P] for r in range(P))
    crofton_total = sum(boundary[r][s] for r in range(P) for s in range(P) if r != s)
    require(natural_mass == alpha_total + beta_total, "natural wall mass")
    require(crofton_total == 12 * alpha_total + 22 * beta_total, "deep Crofton two-mass law")
    require((crofton_total - 12 * natural_mass) // 10 == beta_total, "Crofton recovers pair mass")
    require((22 * natural_mass - crofton_total) // 10 == alpha_total, "Crofton recovers singleton mass")
    require(13 * natural_mass - crofton_total // 2 == 7 * alpha_total + 2 * beta_total, "THM-2529 consumer as Crofton complement")

    return {
        "path_masks": len(masks),
        "gamma_cases": gamma_cases,
        "gamma_integrated_cases": gamma_integrated_cases,
        "beta_point_cases": beta_point_cases,
        "beta_integrated_cases": beta_integrated_cases,
        "natural_wall_cases": natural_wall_cases,
        "alpha_total": alpha_total,
        "beta_total": beta_total,
        "natural_mass": natural_mass,
        "crofton_total": crofton_total,
        "deep_consumer": 7 * alpha_total + 2 * beta_total,
    }


def main():
    boolean = audit_all_boolean_masks()
    gram = audit_anchored_gram_boundary_equivalence()
    deep = audit_deep_path_specialization()

    print("THM-2534 positive Kakeya boundary transform exact referee")
    print(f"prime={P}")
    print(f"nonconstant_masks={ALL - 1}")
    print(f"mask_slope_cases={boolean['mask_slope_cases']}")
    print(f"tail_head_coordinate_cases={boolean['tail_head_coordinate_cases']}")
    print(f"cayley_coordinate_cases={boolean['cayley_coordinate_cases']}")
    print(f"distinct_tail_masks={boolean['distinct_tail_masks']}")
    print(f"primitive_root_mode_cases={boolean['root_mode_cases']}")
    print(f"primitive_root_mode_incidence_cases={boolean['root_mode_incidence_cases']}")
    print(f"reconstruction_coordinates={boolean['reconstruction_coordinates']}")
    print(f"point_gradient_cases={boolean['gradient_point_cases']}")
    print(f"sharp_13_over_42_cases={boolean['sharp_cases']}")
    print(f"flux_histogram={boolean['flux_histogram']}")
    print(f"anchored_masks={gram['anchored_masks']}")
    print(f"anchored_point_entries={gram['point_entries']}")
    print(f"anchored_inverse_entries={gram['inverse_entries']}")
    print(f"metric_entries={gram['metric_entries']}")
    print(f"gradient_entries={gram['gradient_entries']}")
    print(f"zero_curl_cases={gram['zero_curl_cases']}")
    print(f"gram_checksum={gram['gram_checksum']}")
    print(f"boundary_checksum={gram['boundary_checksum']}")
    print(f"deep_path_masks={deep['path_masks']}")
    print(f"deep_gamma_cases={deep['gamma_cases']}")
    print(f"deep_gamma_integrated_cases={deep['gamma_integrated_cases']}")
    print(f"deep_beta_point_cases={deep['beta_point_cases']}")
    print(f"deep_beta_integrated_cases={deep['beta_integrated_cases']}")
    print(f"deep_natural_wall_cases={deep['natural_wall_cases']}")
    print(f"deep_alpha_total={deep['alpha_total']}")
    print(f"deep_beta_total={deep['beta_total']}")
    print(f"deep_natural_mass={deep['natural_mass']}")
    print(f"deep_crofton_total={deep['crofton_total']}")
    print(f"deep_thm2529_consumer={deep['deep_consumer']}")
    print(f"require_checks={checks}")
    print("status=PASS")


if __name__ == "__main__":
    main()
