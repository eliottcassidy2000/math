#!/usr/bin/env python3
"""Exact finite referee for THM-2591.

The only mathematical input from THM-2586 is its pair of proved edge-zero
sets.  Everything below is finite arithmetic in F_13.
"""

from itertools import product


P = 13
CLOCKS = tuple(range(7))
DISPLACEMENTS = tuple(range(1, 13))
MARKERS = tuple(range(1, 13))

# K_ell(s,0,0)=0 and K_ell(s,6,12)=0, respectively.
ZERO_AT_ROOT_0 = {(7, ell) for ell in (4, 5, 6)}
ZERO_AT_ROOT_6 = {(6, ell) for ell in (4, 5, 6)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def allowed_roots(s: int, ell: int) -> tuple[int, ...]:
    roots = []
    if (s, ell) not in ZERO_AT_ROOT_0:
        roots.append(0)
    if (s, ell) not in ZERO_AT_ROOT_6:
        roots.append(6)
    return tuple(roots)


def selectors(s: int):
    return product(*(allowed_roots(s, ell) for ell in CLOCKS))


def coboundary(h: tuple[int, ...]) -> tuple[int, ...]:
    # Edge ell is oriented from ell-1 to ell, matching THM-2542 (5).
    return tuple((h[ell] - h[(ell - 1) % 7]) % P for ell in CLOCKS)


def dihedral_relabels(h: tuple[int, ...]):
    for shift in CLOCKS:
        yield tuple(h[(ell - shift) % 7] for ell in CLOCKS)
        yield tuple(h[(shift - ell) % 7] for ell in CLOCKS)


print("== THM-2591: theta-zero selector versus C91 Cech holonomy ==")

expected_counts = {
    s: (16 if s in (6, 7) else 128) for s in DISPLACEMENTS
}
count_by_s = {}
selector_bank = []
for s in DISPLACEMENTS:
    bank = list(selectors(s))
    require(
        all(len(allowed_roots(s, ell)) > 0 for ell in CLOCKS),
        f"empty selector fibre at displacement {s}",
    )
    require(
        len(bank) == expected_counts[s],
        f"wrong selector count at displacement {s}",
    )
    count_by_s[s] = len(bank)
    selector_bank.extend((s, h) for h in bank)

require(sum(count_by_s.values()) == 1312, "wrong total selector count")
print("admissible selector counts by s:", count_by_s)
print("total admissible selectors:", len(selector_bank))

coboundary_checks = 0
marker_selector_checks = 0
invoice_checks = 0
for _, h in selector_bank:
    dh = coboundary(h)
    require(sum(dh) % P == 0, "selector coboundary has nonzero holonomy")
    coboundary_checks += 1

    for a in MARKERS:
        gauged = tuple((a + d) % P for d in dh)
        holonomy = sum(gauged) % P
        require(holonomy == (7 * a) % P, "gauge changed base holonomy")
        require(holonomy != 0, "nonzero marker acquired zero holonomy")
        require(any(gauged), "nonzero holonomy had all-zero edge law")

        # If an additional mixed-square correction c flattened every edge,
        # then c=-gauged and its cyclic sum is forced to be -7a.
        correction = tuple((-g) % P for g in gauged)
        require(
            sum(correction) % P == (-7 * a) % P,
            "wrong mixed-cell boundary invoice",
        )
        marker_selector_checks += 1
        invoice_checks += 1

require(coboundary_checks == 1312, "wrong coboundary check count")
require(marker_selector_checks == 1312 * 12, "wrong marker check count")
print("selector coboundary checks:", coboundary_checks)
print("marker/selector holonomy checks:", marker_selector_checks)
print("mixed-cell boundary invoices:", invoice_checks)

# Rotation/reflection of the proposed clock identification cannot alter the
# conclusion.  Reflection can reverse the class, but cannot make it zero.
dihedral_checks = 0
for _, h in selector_bank:
    for h2 in dihedral_relabels(h):
        dh2 = coboundary(h2)
        require(sum(dh2) % P == 0, "relabelled selector is not exact")
        for a in MARKERS:
            require(
                sum((a + d) % P for d in dh2) % P == (7 * a) % P,
                "dihedral relabelling changed holonomy",
            )
            dihedral_checks += 1
require(dihedral_checks == 1312 * 14 * 12, "wrong dihedral check count")
print("dihedral clock-relabel checks:", dihedral_checks)

base_holonomies = {a: (7 * a) % P for a in MARKERS}
require(
    set(base_holonomies.values()) == set(MARKERS),
    "multiplication by seven is not a unit permutation",
)
print("base holonomies a -> 7a:", base_holonomies)

minimal_cover_degrees = {}
for a in MARKERS:
    degree = next(n for n in range(1, P + 1) if n * 7 * a % P == 0)
    require(degree == 13, "wrong minimal trivializing cover degree")
    minimal_cover_degrees[a] = degree
print("minimal trivializing cover degrees:", minimal_cover_degrees)

print("SCOPE: finite exact consequence of proved THM-2586 zero sets and")
print("THM-2542 transition law; no physical mixed 2-cell or row closure.")
print("ALL CHECKS PASSED")
