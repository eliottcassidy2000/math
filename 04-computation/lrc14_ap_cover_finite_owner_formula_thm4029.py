#!/usr/bin/env python3
"""Finite-owner formula for the AP seven-sector cover deficit.

For m>=13 the exact non-cover set has one component around each reduced
p/q on the circle with q<=6.  Its two radii are finite max-min expressions
using the largest time e<m in each congruence class modulo q.
"""

from fractions import Fraction as Q
from math import floor, gcd, lcm
import sys

from lrc14_ap_cover_components_thm4029 import (
    covers_theta,
    local_missing_sector_radii,
    noncover_components,
)
from lrc14_ap_cover_sequence_engine_thm4029 import farey_cover_distribution


def require(condition: bool, label):
    if not condition:
        raise RuntimeError(label)


def owners():
    return tuple(
        Q(p, q)
        for q in range(1, 7)
        for p in range(q)
        if gcd(p, q) == 1
    )


def owner_formula_deficit(m: int) -> Q:
    total_theta = Q(0)
    for x0 in owners():
        minus, plus, _, _ = local_missing_sector_radii(x0.numerator, x0.denominator, m)
        total_theta += minus + plus
    return total_theta / 7


def direct_component_deficit(m: int) -> Q:
    return sum((right - left for left, right, _ in noncover_components(m)), Q(0)) / 7


def prove_local_guards(base_m: int):
    """Exact guards needed by the track argument on every C_base component."""
    components = noncover_components(base_m)
    records = []
    for x0 in owners():
        theta0 = 7 * x0
        match = None
        for left, right, _ in components:
            lifted = theta0 + (7 if right > 7 and theta0 < left else 0)
            if left <= lifted <= right:
                match = (left, right, lifted)
                break
        require(match is not None, (base_m, x0, "owner component missing"))
        left, right, lifted = match
        minus, plus = lifted - left, right - lifted
        p, q = x0.numerator, x0.denominator
        require(q * minus < 1 and q * plus < 1, (base_m, x0, "track step can skip"))
        margins = []
        for s in range(1, q):
            A = Q(7 * p * s, q)
            frac = A - floor(A)
            margins.append((1 - frac) - s * plus)
            margins.append(frac - s * minus)
        require(all(margin >= 0 for margin in margins), (base_m, x0, "start sector lost"))
        records.append((x0, minus, plus, min(margins) if margins else Q(1)))
    return tuple(records)


def audit_m13_zero_margin_endpoints():
    """Independently audit the only two zero guard margins at the sharp base."""
    m = 13
    bps = {Q(0), Q(7)}
    for e in range(1, m):
        for r in range(7 * e + 1):
            bps.add(Q(r, e))
    bps = sorted(bps)
    records = []
    for endpoint in (Q(16, 3), Q(17, 4)):
        i = bps.index(endpoint)
        left_mid = (bps[i - 1] + endpoint) / 2
        right_mid = (endpoint + bps[i + 1]) / 2
        record = (
            endpoint,
            bps[i - 1],
            bps[i + 1],
            covers_theta(left_mid, m),
            covers_theta(endpoint, m),
            covers_theta(right_mid, m),
        )
        require(record[-3:] == (False, True, True), (endpoint, "endpoint convention changed"))
        records.append(record)
    return tuple(records)


def radius_signature(p: int, q: int, m: int):
    n = m - 1
    data = []
    occupied = set()
    for s in range(q):
        A = Q(7 * p * s, q)
        integer = floor(A)
        frac = A - integer
        residue = integer % 7
        occupied.add(residue)
        E = n - ((n - s) % q)
        if E == 0:
            E += q
        data.append((s, residue, frac, E))
    missing = set(range(7)) - occupied

    side_signatures = []
    for side in ("minus", "plus"):
        by_r = {}
        for r in missing:
            vals = []
            for s, residue, frac, E in data:
                if side == "plus":
                    d = (r - residue) % 7
                    value = (Q(d) - frac) / E
                else:
                    d = (residue - r) % 7
                    value = (Q(d - 1) + frac) / E
                vals.append((value, s))
            vmin = min(v for v, _ in vals)
            by_r[r] = (vmin, tuple(s for v, s in vals if v == vmin))
        vmax = max(v for v, _ in by_r.values())
        side_signatures.append(
            (tuple(r for r, (v, _) in by_r.items() if v == vmax),
             tuple((r, by_r[r][1]) for r in sorted(by_r)))
        )
    return tuple(side_signatures)


def _le_forever(f, g, n0: int) -> bool:
    """Prove C/(n-c) <= D/(n-d) for every real n>=n0."""
    C, c = f
    D, d = g
    # H(n)=C(n-d)-D(n-c) must be <=0.
    slope = C - D
    value0 = C * (n0 - d) - D * (n0 - c)
    return value0 <= 0 and slope <= 0


def _ge_forever(f, g, n0: int) -> bool:
    return _le_forever(g, f, n0)


def stable_radius_term(p: int, q: int, n0: int, side: str):
    """Return a proved fixed C/(n-c) max-min winner for all n>=n0 in phase."""
    data = []
    occupied = set()
    for s in range(q):
        A = Q(7 * p * s, q)
        integer = floor(A)
        frac = A - integer
        residue = integer % 7
        occupied.add(residue)
        c = (n0 - s) % q  # E_s=n-c throughout this mod-q phase
        data.append((s, residue, frac, c))
    missing = sorted(set(range(7)) - occupied)
    inner = {}
    for r in missing:
        candidates = []
        for s, residue, frac, c in data:
            if side == "plus":
                C = Q((r - residue) % 7) - frac
            else:
                C = Q(((residue - r) % 7) - 1) + frac
            candidates.append(((C, c), s))
        winner = next(
            ((f, s) for f, s in candidates if all(_le_forever(f, g, n0) for g, _ in candidates)),
            None,
        )
        if winner is None:
            raise AssertionError((p, q, n0, side, r, "inner selector crosses"))
        inner[r] = winner
    outer = next(
        ((f, r, s) for r, (f, s) in inner.items()
         if all(_ge_forever(f, g, n0) for g, _ in inner.values())),
        None,
    )
    if outer is None:
        raise AssertionError((p, q, n0, side, "outer selector crosses"))
    return outer


def prove_phase_rational_law(period: int = 60, minimum_n: int = 12):
    laws = {}
    asymptotic_constants = set()
    for phase in range(period):
        n0 = phase
        while n0 < minimum_n:
            n0 += period
        terms = []
        selectors = []
        for x0 in owners():
            p, q = x0.numerator, x0.denominator
            for side in ("minus", "plus"):
                f, r, s = stable_radius_term(p, q, n0, side)
                terms.append(f)
                selectors.append((x0, side, r, s, f))
        laws[phase] = (n0, tuple(terms), tuple(selectors))
        asymptotic_constants.add(sum((C for C, _ in terms), Q(0)) / 7)
    return laws, asymptotic_constants


def asymptotic_constant() -> Q:
    """lim m(1-a(m)), obtained by replacing every E_s(m) by m."""
    total_theta = Q(0)
    for x0 in owners():
        p, q = x0.numerator, x0.denominator
        data = []
        occupied = set()
        for s in range(q):
            A = Q(7 * p * s, q)
            integer = floor(A)
            frac = A - integer
            residue = integer % 7
            occupied.add(residue)
            data.append((residue, frac))
        missing = set(range(7)) - occupied
        plus = max(min(Q((r-residue) % 7) - frac for residue, frac in data) for r in missing)
        minus = max(min(Q(((residue-r) % 7) - 1) + frac for residue, frac in data) for r in missing)
        total_theta += plus + minus
    return total_theta / 7


def main():
    sys.stdout.reconfigure(line_buffering=True)
    print(f"owners={owners()}")
    guard_records_13 = prove_local_guards(13)
    guard_records_14 = prove_local_guards(14)
    print(f"m13_track_guard_records={guard_records_13}")
    print(f"m14_track_guard_records={guard_records_14}")
    print(f"m13_zero_margin_endpoint_audit={audit_m13_zero_margin_endpoints()}")
    print("\nINDEPENDENT COMPONENT/FIFTEEN-TRACK FORMULA CHECK")
    first_failure = direct_component_deficit(12) != owner_formula_deficit(12)
    print(
        f"m=12 minimal-tail-formula-failure={first_failure}; "
        f"direct={direct_component_deficit(12)} formula={owner_formula_deficit(12)}"
    )
    require(first_failure, "m=12 immediate-predecessor hostile")
    for m in [13, 14, 15, 20, 27, 28, 40, 60, 80, 120]:
        direct = direct_component_deficit(m)
        formula = owner_formula_deficit(m)
        print(f"m={m:3d}: direct={direct} formula={formula} equal={direct == formula}")
        require(direct == formula, (m, "direct/owner formula"))

    mass, unresolved, _ = farey_cover_distribution(200)
    cumulative = Q(0)
    mismatches = []
    for m in range(1, 201):
        cumulative += mass[m]
        if m >= 13 and 1 - cumulative != owner_formula_deficit(m):
            mismatches.append((m, 1 - cumulative, owner_formula_deficit(m)))
    print(
        f"farey_cover_time_vs_owner_formula_m13_200={not mismatches}; "
        f"mismatches={mismatches[:3]}; unresolved_at_200={unresolved}"
    )
    require(not mismatches, "Farey/owner formula tail")

    period = lcm(*range(1, 7))
    print(f"\nPHASE-SELECTOR STABILITY, candidate period={period}")
    changes = []
    for x0 in owners():
        p, q = x0.numerator, x0.denominator
        for residue in range(period):
            ms = [m for m in range(14, 20001) if (m - 1) % period == residue]
            if not ms:
                continue
            baseline = radius_signature(p, q, ms[0])
            for m in ms[1:]:
                if radius_signature(p, q, m) != baseline:
                    changes.append((x0, residue, ms[0], m))
                    break
    print(f"selector_changes_on_same_mod60_phase_through_20000={changes}")
    require(not changes, "same-phase selector stability")

    laws, phase_constants = prove_phase_rational_law(period, minimum_n=12)
    print(
        f"symbolically_proved_non_crossing_phase_laws={len(laws)}; "
        f"phase_asymptotic_constants={phase_constants}"
    )
    require(phase_constants == {Q(127, 35)}, "tail constant")

    C = asymptotic_constant()
    print(f"\nASYMPTOTIC CONSTANT C=lim m(1-a(m))={C}={float(C):.12f}")
    for m in (40, 80, 120, 500, 1000, 10000):
        d = owner_formula_deficit(m)
        print(f"m={m:5d}: deficit={d} m*deficit={float(m*d):.12f} error={float(m*d-C):+.3e}")

    print("\nSIMPLE STRICTNESS/TAIL BOUNDS")
    for m in (8, 14, 28, 100, 1000):
        lower_tail = Q(11, 7 * (m - 1))
        actual = owner_formula_deficit(m) if m >= 13 else None
        print(f"m={m}: 1-a(m) >= {lower_tail}; exact={actual}")


if __name__ == "__main__":
    main()
