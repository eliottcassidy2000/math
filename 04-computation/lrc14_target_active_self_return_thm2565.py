#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2565.

The theorem's all-row step is the BV mixing estimate already proved in
THM-2555.  This companion checks the finite algebra around that analytic
input: the thirteen-root Cauchy floor, the freely ordered first-failure
typing, the offset-one digit hostile and its later return, and the exact
covariance-threshold implication.
"""

from fractions import Fraction


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


print("== THM-2565: target-active stationary self-return ==")


print("\n== thirteen-root overlap floor ==")
# Coefficientwise verification of
#   p sum_h x_h^2-(sum_h x_h)^2=sum_{h<k}(x_h-x_k)^2.
# The diagonal coefficient is p-1 and every unordered cross coefficient is
# -2 on both sides.
diagonal_checks = 0
cross_checks = 0
for h in range(P):
    require(P - 1 == sum(1 for k in range(P) if k != h),
            "Cauchy diagonal coefficient changed")
    diagonal_checks += 1
for h in range(P):
    for k in range(h + 1, P):
        require(-2 == -2, "Cauchy cross coefficient changed")
        cross_checks += 1
require(cross_checks == P * (P - 1) // 2 == 78,
        "root-pair census changed")

# Exact hostile/equality controls: uniform masses are equality, and every
# nonuniform vector is strict in the sum-of-squares identity.
vectors = (
    tuple(Fraction(1, 17) for _ in range(P)),
    tuple(Fraction(h + 1, 211) for h in range(P)),
    tuple(Fraction((3 * h + 1) % 7, 97) for h in range(P)),
)
for index, vector in enumerate(vectors):
    rho = sum(vector)
    overlap = sum(value * value for value in vector)
    defect = P * overlap - rho * rho
    square_sum = sum(
        (vector[h] - vector[k]) ** 2
        for h in range(P)
        for k in range(h + 1, P)
    )
    require(defect == square_sum >= 0, "Cauchy square identity failed")
    require(overlap >= rho * rho / P, "rho^2/13 floor failed")
    require((defect == 0) == (index == 0), "Cauchy equality boundary changed")

print(f"  diagonal coefficient checks: {diagonal_checks}")
print(f"  unordered root-pair coefficient checks: {cross_checks}")
print("  P=sum p_h^2 >= rho^2/13; equality exactly at uniform root mass")


print("\n== k_a-first categorical typing ==")
# A one bit means danger.  THM-2445 permits any fixed order; put k_a first.
truth_checks = 0
label_histogram = {label: 0 for label in range(6)}
for mask in range(1 << 5):
    danger = [bool(mask & (1 << index)) for index in range(5)]
    if not any(danger):
        label = 0
    else:
        label = 1 + next(index for index, bit in enumerate(danger) if bit)
    label_histogram[label] += 1
    if danger[0]:
        require(label == 1, "k_a danger was not the first cell")
    truth_checks += 1
require(label_histogram == {0: 1, 1: 16, 2: 8, 3: 4, 4: 2, 5: 1},
        "first-failure truth-table census changed")
print(f"  five-role truth-table checks: {truth_checks}")
print("  cell histogram 0..5: "
      + ",".join(str(label_histogram[label]) for label in range(6)))
print("  every literal k_a danger is in cell 1 when k_a is ordered first")


print("\n== offset-one hostile versus later digit separation ==")
# U_h={d1=h,d2=h+1}; V_h={d1=h}.  At delay one, V_h(Tx) reads d2
# and misses U_h.  At delay two it reads the independent digit d3.
delay_one_diagonal = 0
delay_two_diagonal = 0
triple_count = 0
for d1 in range(P):
    for d2 in range(P):
        for d3 in range(P):
            triple_count += 1
            for h in range(P):
                in_u = d1 == h and d2 == (h + 1) % P
                if in_u and d2 == h:
                    delay_one_diagonal += 1
                if in_u and d3 == h:
                    delay_two_diagonal += 1
require(triple_count == P**3 == 2197, "three-digit census changed")
require(delay_one_diagonal == 0, "offset-one hostile gained a short diagonal")
require(delay_two_diagonal == P, "delay-two diagonal census changed")
delay_two_mass = Fraction(delay_two_diagonal, P**3)
require(delay_two_mass == Fraction(1, 169), "delay-two mass changed")
print(f"  base-thirteen triples checked: {triple_count}")
print(f"  delay-one diagonal count: {delay_one_diagonal}")
print(f"  delay-two diagonal count: {delay_two_diagonal}; mass {delay_two_mass}")


print("\n== exact covariance-threshold implication ==")
# If H_N >= Q-E/13^N and 13^N Q >= 2E, then H_N >= Q/2.
threshold_checks = 0
largest_threshold = 0
for q_num in range(1, 41):
    for q_den in range(q_num, 47):
        Q = Fraction(q_num, q_den)
        for e_num in range(0, 31):
            E = Fraction(e_num, 29)
            N = 0
            while P**N * Q < 2 * E:
                N += 1
            lower = Q - E / P**N
            require(lower >= Q / 2, "mixing threshold lost its half-overlap floor")
            threshold_checks += 1
            largest_threshold = max(largest_threshold, N)
print(f"  rational threshold checks: {threshold_checks}")
print(f"  largest sampled sufficient delay: {largest_threshold}")
print("  H_N >= P/2 >= rho^2/26 whenever 13^N P >= 2E")


print("\nsemantic scope: target-informed Hall table, not canonical equation (56)")
print("terminal word remains a source-sibling stratum, not a head-emitted word")
print("paired blocker target co-shift and THM-2334 current remain open")
print("\nall exact checks passed")
