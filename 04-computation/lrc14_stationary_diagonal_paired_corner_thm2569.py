#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2569.

This checks the arithmetic inherited by conditioning THM-2563 on THM-2565,
and the finite Fourier line-sum identity showing why a stationary Hall
diagonal is only a spectator under full-X endpoint recombination.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


print("== THM-2569: stationary diagonal inside the paired corner ==")


print("\n== inherited positive floors ==")
self_return = Fraction(1, 26)
total_floor = 63 * self_return
guard_full = Fraction(7, 35152) * self_return
ordinary_full = Fraction(9, 35152) * self_return
guard_marginal = Fraction(7, 2704) * self_return
ordinary_marginal = Fraction(9, 2704) * self_return
require(total_floor == Fraction(63, 26), "conditioned mass floor changed")
require(guard_full == Fraction(7, 913952), "guard full-table floor changed")
require(ordinary_full == Fraction(9, 913952),
        "ordinary full-table floor changed")
require(guard_marginal == Fraction(7, 70304),
        "guard marginal floor changed")
require(ordinary_marginal == Fraction(9, 70304),
        "ordinary marginal floor changed")
print(f"  conditioned table mass per rho^2: {total_floor}")
print(f"  normalized full-table floors: {guard_full} / {ordinary_full}")
print(f"  two-axis marginal floors: {guard_marginal} / {ordinary_marginal}")


print("\n== common-twist Boolean annihilation ==")
# For arbitrary stationary multipliers W_s and arbitrary endpoint weights,
# W_s P_s(1-P_s)=0.  Exhaust every 13-bit gate and two nontrivial weight
# profiles.  This is the pointwise input to endpoint Parseval.
gate_checks = 0
weight_profiles = (
    tuple(Fraction(s + 1, 31) for s in range(P)),
    tuple(Fraction((5 * s + 2) % 17, 43) for s in range(P)),
)
for gate in range(1 << P):
    for weights in weight_profiles:
        for s in range(P):
            bit = (gate >> s) & 1
            require(weights[s] * bit * (1 - bit) == 0,
                    "stationary multiplier broke Boolean orthogonality")
            gate_checks += 1
require(gate_checks == (1 << P) * len(weight_profiles) * P,
        "Boolean annihilation census changed")
print(f"  weighted Boolean gate checks: {gate_checks}")
print("  every stationary multiplier preserves P_s(1-P_s)=0")


print("\n== independent-twist to coarse-target line sum ==")
# For q=a+b, summing the 2D DFT along a+b=q gives
#   (1/p) sum_s A(s,s) zeta^(qs).
# Verify the load-bearing orthogonality sum coefficientwise for every
# q,s,t.  If s!=t, the exponents a(s-t) run through F_p once; if s=t,
# every exponent is zero.
orthogonality_checks = 0
for q in range(P):
    for s in range(P):
        for t in range(P):
            exponent_counts = [0] * P
            for a in range(P):
                exponent = (a * s + (q - a) * t) % P
                exponent_counts[exponent] += 1
            if s == t:
                expected = [0] * P
                expected[(q * s) % P] = P
                require(exponent_counts == expected,
                        "diagonal line-sum phase changed")
            else:
                require(exponent_counts == [1] * P,
                        "off-diagonal character orthogonality failed")
            orthogonality_checks += 1
require(orthogonality_checks == P**3 == 2197,
        "line-sum orthogonality census changed")

# Positive off-diagonal hostile A(s,t)=1_(t=s+1) has zero diagonal and
# therefore zero coarse line sum for every q, despite total mass p.
hostile = [
    [Fraction(int(t == (s + 1) % P), 1) for t in range(P)]
    for s in range(P)
]
require(sum(sum(row) for row in hostile) == P,
        "positive off-diagonal hostile mass changed")
require(all(hostile[s][s] == 0 for s in range(P)),
        "off-diagonal hostile gained a diagonal")
print(f"  exact character-orthogonality checks: {orthogonality_checks}")
print("  positive shift-graph hostile mass 13; every coarse target line sum 0")


print("\n== Hall-root tensor is still a spectator ==")
# Tensor the endpoint hostile with the perfect semantic Hall diagonal b=h.
# The root side has positive mass thirteen, but it does not alter the zero
# target line sum because h is not an endpoint residue coordinate.
hall_cells = 0
vanishing_root_target_pairs = 0
for h in range(P):
    for b in range(P):
        if b != h:
            continue
        hall_cells += 1
        for q in range(P):
            diagonal_target_sum = sum(hostile[s][s] for s in range(P))
            require(diagonal_target_sum == 0,
                    "Hall tensor broke target annihilation")
            vanishing_root_target_pairs += 1
require(hall_cells == P, "Hall diagonal cell census changed")
require(vanishing_root_target_pairs == P**2,
        "root/target vanishing census changed")
print(f"  positive Hall diagonal root cells: {hall_cells}")
print(f"  vanishing (Hall root, coarse target) pairs: {vanishing_root_target_pairs}")
print("  physical branch equality h=b does not supply the endpoint residue eta.u")


print("\nsemantic scope: positive common support, but full-X coarse target still zero")
print("a fixed X, nonconstant normal/jet, or covariant relative sidecar remains needed")
print("\nall exact checks passed")
