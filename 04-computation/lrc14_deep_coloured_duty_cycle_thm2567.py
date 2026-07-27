#!/usr/bin/env python3
"""Exact finite referee for THM-2567's deep-coloured duty cycle."""

from __future__ import annotations


P = 13
ZERO = (0,) * (P - 1)
CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


def add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(x + y for x, y in zip(a, b))


def scale(c: int, a: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(c * x for x in a)


def zeta_pow(exponent: int) -> tuple[int, ...]:
    exponent %= P
    if exponent < P - 1:
        ans = [0] * (P - 1)
        ans[exponent] = 1
        return tuple(ans)
    return (-1,) * (P - 1)


def hostile(r: int, s: int, t: int) -> int:
    return int(s != 0 and r != t)


def fully_anchored_hostile(r: int, s: int, t: int) -> int:
    return int(r != 0 and s != 0 and r != t)


def transform(function, m: int, b: int, h: int) -> tuple[int, ...]:
    """Unnormalized three-variable transform with positive exponents."""
    ans = ZERO
    for r in range(P):
        for s in range(P):
            for t in range(P):
                value = function(r, s, t)
                if value:
                    ans = add(ans, scale(value, zeta_pow(m * r + b * s + h * t)))
    return ans


def duty(q: tuple[int, int]) -> int:
    x, y = q
    return (
        2316060
        + 210552 * (int(x == 0) + int(y == 0))
        + 12 * int((x + y) % P == 0)
        + 19128 * int(x == 0 and y == 0)
    )


POINTS = [(x, y) for x in range(P) for y in range(P)]


def add_point(q: tuple[int, int], g: tuple[int, int], coefficient: int) -> tuple[int, int]:
    return ((q[0] + coefficient * g[0]) % P, (q[1] + coefficient * g[1]) % P)


def commutator(a: dict[tuple[int, int], int], g: tuple[int, int]) -> dict[tuple[int, int], int]:
    return {
        q: sum(
            (duty(q) - duty(add_point(q, g, -j))) * a[add_point(q, g, -j)]
            for j in range(7)
        )
        for q in POINTS
    }


print("== THM-2567: deep-coloured duty-replica cycle ==")
hostile_mass = sum(hostile(r, s, t) for r in range(P) for s in range(P) for t in range(P))
require(hostile_mass == P * (P - 1) * (P - 1), "hostile mass drifted")
require(all(hostile(t, 0, t) == 0 for t in range(P)), "s=0 plane failed")
require(all(hostile(t, s, t) == 0 for s in range(P) for t in range(P)), "diagonal failed")
print(f"  rational nonnegative hostile mass: {hostile_mass}")
print("  exact zero loci: s=0 plane and r=t diagonal")

print("\n== physical target reindexing and inverse-character anchors ==")
all_a: dict[int, dict[tuple[int, int], int]] = {}
for m in range(P):
    table: dict[tuple[int, int], int] = {}
    for b in range(P):
        for q2 in range(P):
            actual = transform(hostile, m, b, q2 - m)
            if q2 != 0:
                expected = ZERO
                scalar = 0
            else:
                s_sum = P - 1 if b == 0 else -1
                scalar = P * (P - 1) * s_sum if m == 0 else -P * s_sum
                expected = scale(scalar, zeta_pow(0))
            require(actual == expected, f"hostile transform drifted at {(m,b,q2)}")
            table[(b, q2)] = scalar
    all_a[m] = table

for m in range(1, P):
    require(all_a[m][(0, 0)] == -156, "nonzero deep anchor drifted")
    require(all(all_a[m][(b, 0)] == 13 for b in range(1, P)), "off-target hostile drifted")
print("  for every m!=0: p^3 A_m(0,0)=-156")
print("  for every m,b!=0: p^3 A_m(b,0)=13")
print("  all hostile support lies on the physical q2=0 line")

for q in POINTS:
    require(sum(all_a[m][q] for m in range(P)) == 0, f"deep circulation failed at {q}")
print("  target-null circulation sum_m A_m(q)=0: 169/169 PASS")

print("\n== canonical duty replicas and exact coloured cancellation ==")
gains = [q for q in POINTS if q != (0, 0)]
duty_classes: dict[int, int] = {}
commutators: dict[tuple[int, int], list[dict[tuple[int, int], int]]] = {}
for g in gains:
    d_g = duty((0, 0)) - duty(g)
    duty_classes[d_g] = duty_classes.get(d_g, 0) + 1
    kg_by_m = [commutator(all_a[m], g) for m in range(P)]
    commutators[g] = kg_by_m
    for m in range(1, P):
        anchor = all_a[m][(0, 0)]
        for j in range(1, 7):
            require(
                kg_by_m[m][add_point((0, 0), g, j)] == -d_g * anchor,
                f"six-replica identity failed at {(g,m,j)}",
            )
    for q in POINTS:
        require(sum(kg_by_m[m][q] for m in range(P)) == 0, f"K cycle failed at {(g,q)}")

require(duty_classes == {229692: 24, 440232: 12, 440244: 132}, "duty classes drifted")
print("  gain duty classes: 229692^24, 440232^12, 440244^132")
print("  every m!=0 and every gain has six identical nonzero replicas")
print("  sum_m K_g A_m=0 for all 168 gains and 169 targets")

print("\n== stronger zero-plane and singleton controls ==")
anchored_mass = sum(
    fully_anchored_hostile(r, s, t)
    for r in range(P) for s in range(P) for t in range(P)
)
require(anchored_mass == (P - 1) * (P - 1) * (P - 1), "anchored mass drifted")
g_profile = [
    sum(
        fully_anchored_hostile((u + t) % P, s, t)
        for s in range(P) for t in range(P)
    )
    for u in range(P)
]
require(g_profile[0] == 0 and all(x > 0 for x in g_profile[1:]), "anchored profile failed")
for m in range(1, P):
    coefficient = ZERO
    for u, value in enumerate(g_profile):
        coefficient = add(coefficient, scale(value, zeta_pow(m * u)))
    require(coefficient != ZERO, f"fully anchored deep colour vanished at {m}")

singleton_controls = 0
for r in range(P):
    for s in range(P):
        for t in range(P):
            if r == t:
                continue
            for m in range(1, P):
                require(zeta_pow(m * (r - t)) != ZERO, "singleton anchor vanished")
            singleton_controls += 1
require(singleton_controls == P * P * (P - 1), "singleton control count drifted")
print(f"  r=0, s=0, and r=t hostile mass: {anchored_mass}")
print("  all 12 inverse-character anchors survive")
print(f"  diagonal-free singleton controls: {singleton_controls}")

print(f"\nexplicit checks: {CHECKS}")
print("ALL EXACT CHECKS PASSED")
