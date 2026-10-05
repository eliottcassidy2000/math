#!/usr/bin/env python3
"""Exact controls for the accompanying three-bits / critical-flow note.

Standard library only. All checks survive python -O. No convergence assumption
enters the identities; the finite orbit census explicitly reports its universe.
"""
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

CHECKS = Counter()


def check(tag, condition, witness=None):
    CHECKS[tag] += 1
    if not condition:
        raise RuntimeError((tag, witness))


def v2(n):
    if n <= 0:
        raise ValueError(n)
    return (n & -n).bit_length() - 1


def vp(n, p):
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def odd(n):
    return n >> v2(n)


def U(n, sign=1):
    return odd(3 * n + sign)


def nu(n):
    return F(8, 3 * 4 ** n.bit_length())


def mu(n):
    m = (n + 1) // 2
    return F(1, 2 ** (2 * (m.bit_length() - 1) + 1))


def q(a):
    return F(3, 4 ** a)


def cofactor_mass(t):
    return F(8, 9) if t == 1 else nu(t) / 3


def climb_closed(n):
    """Exact resolvent of the rising-edge operator, on odd n > 1."""
    total = nu(n)
    x = n
    while (x + 1) % 3 == 0:
        x = 2 * (x + 1) // 3 - 1
        total += nu(x)
    return total


def peak(n):
    h = v2(n + 1)
    t = (n + 1) >> h
    return 2 * 3 ** (h - 1) * t - 1


def block(n):
    h = v2(n + 1)
    t = (n + 1) >> h
    return odd(3 ** h * t - 1)


def golden_b(n):
    m = n + 1
    return (3 * m - isqrt(5 * m * m) - 1) // 2


def charge(n):
    b = golden_b(n)
    return n - 2 * b, b


def zeck_digits(n):
    fs = [1, 2]
    while fs[-1] <= n:
        fs.append(fs[-1] + fs[-2])
    out = []
    for f in reversed(fs):
        if out or f <= n:
            d = int(f <= n)
            out.append(d)
            n -= d * f
    return out or [0]


def main():
    report = {"status": "FINITE-EXACT; symbolic proofs are in the companion note",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    limit = 1 << 16
    for n in range(1, limit, 2):
        m = (n + 1) // 2
        ratio = F(4, 3) if m & (m - 1) == 0 else F(1, 3)
        check("measure_comparison", nu(n) == ratio * mu(n), n)
        h = v2(n + 1)
        t = (n + 1) >> h
        check("independent_climb_cofactor", nu(n) == q(h) * cofactor_mass(t), n)
        for j in range(7):
            height = n.bit_length() + j
            joint = F(2, 4 ** height)
            check("independent_alias", joint == nu(n) * q(j + 1), (n, j))
    # Independent summation of complete finite dyadic shells plus exact tails.
    for d in range(1, 13):
        head = list(range(1, 1 << d, 2))
        check("nu_normalization", sum(map(nu, head)) + F(2, 3 * 2 ** d) == 1, d)
        for h in range(1, d + 1):
            cell = sum((nu(n) for n in head if v2(n + 1) == h), F())
            # For h<d every later shell has residue fraction 2^-h.
            # At h=d the very first later shell contains no such residue.
            tail = F(2, 3 * 2 ** d) * F(1, 2 ** h)
            if h == d:
                tail /= 2
            check("climb_law_by_shells", cell + tail == q(h), (d, h))
    report["source_laws"] = {
        "odd_source_bound_exclusive": limit, "alias_offsets": "0..6",
        "nu": "8/(3*4^bit_length(n))", "q": "3/4^a",
        "entropy_mu": "3", "entropy_nu": "log2(3)+1/3",
        "entropy_alias": "8/3-log2(3)",
        "nu_all_ones_mass": "8/9", "mu_all_ones_mass": "2/3",
        "cofactor_1_mass": "8/9", "cofactor_other_mass": "nu(t)/3"}
    # Geometric normalization, first moment and the collision tilt, with tails.
    for k in range(1, 41):
        check("q_tail", sum((q(a) for a in range(1, k + 1)), F()) + F(1, 4 ** k) == 1)
        check("q_moment", sum((a * q(a) for a in range(1, k + 1)), F())
              + F(3 * k + 4, 3 * 4 ** k) == F(4, 3))
        check("collision_tilt", F(1, 2 ** k) ** 2 / F(1, 3) == q(k))

    rays = []
    for j in range(1, 41):
        a = 2 * (4 ** j - 1) // 3
        b = a // 2
        m = (b + 1) // 2
        check("root_ray", 3 * b + 1 == 4 ** j and U(b) == 1, j)
        check("other_sheet_ray", 3 * m - 1 == 2 ** (2 * j - 1) and U(m, -1) == 1, j)
        check("ray_binary_word", bin(a)[2:] == "10" * j, j)
        check("ray_ternary_step", 2 * (4 ** (j + 1) - 1) // 3 == 4 * a + 2, j)
        if j <= 6:
            rays.append({"j": j, "height": a, "plus_root_source": b, "source_index": m})
    for d in range(1, 8):
        modulus = 3 ** d
        residues = []
        x = 0
        for j in range(modulus):
            residues.append(x)
            x = (4 * x + 2) % modulus
        check("ternary_odometer", x == 0 and len(set(residues)) == modulus, d)
    k = 40
    check("ray_nu_mass", sum((nu((4 ** j - 1) // 3) for j in range(1, k + 1)), F())
          + F(32, 45 * 16 ** k) == F(32, 45))
    check("ray_mu_mass", sum((mu((4 ** j - 1) // 3) for j in range(1, k + 1)), F())
          + F(32, 15 * 16 ** k) == F(19, 30))
    report["root_ray"] = {"rows": rays, "nu_mass": "32/45", "mu_mass": "19/30",
                          "ternary_bijections_through_depth": 7}

    # Eight-state digit reader: two charge bits and the preceding digit bit.
    states = {(0, 0, 0)}
    queue = deque(states)
    edges = []
    while queue:
        x, y, prev = queue.popleft()
        for d in (0, 1):
            if prev == d == 1:
                continue
            target = ((y + d) % 2, (x + y) % 2, d)
            edges.append([[x, y, prev], d, list(target)])
            if target not in states:
                states.add(target)
                queue.append(target)
    check("eight_states", len(states) == 8 and len(edges) == 12, states)
    carry_examples = {}
    for n in range(20000):
        x = y = prev = 0
        for d in zeck_digits(n):
            check("zeck_guard", not (prev == d == 1), n)
            x, y = y + d, x + y + 2 * d
            prev = d
        check("zeck_reader", x == n and y == 2 * n - golden_b(n), n)
        check("zeck_charge", charge(n) == (2 * y - 3 * x, 2 * x - y), n)
        eta = golden_b(3 * n + 1) - 3 * golden_b(n)
        a, b = charge(n)
        check("golden_carry", eta in (-1, 0, 1, 2)
              and charge(3 * n + 1) == (3 * a + 1 - 2 * eta, 3 * b + eta), n)
        if n % 2:
            carry_examples.setdefault(eta, n)
    check("all_four_carries", set(carry_examples) == {-1, 0, 1, 2})
    plus_parities = []
    for n in range(8):
        x, word = n, 0
        for j in range(3):
            word += (x % 2) << j
            x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
        plus_parities.append(word)
    check("fano_carry_control", plus_parities == [0, 5, 2, 3, 4, 1, 6, 7])
    report["zeckendorf"] = {"finite_universe": "0 <= n < 20000", "states": sorted(states),
                            "transitions": sorted(edges), "odd_carry_examples": carry_examples,
                            "plus_parity_map_mod_8": plus_parities,
                            "scope": "auxiliary XOR charge, not the owner-confirmed stripe"}

    # The closed ascent flow has an independently checked recurrence.
    first_fall_failure = None
    max_ratio = (F(), None)
    for n in range(3, limit, 2):
        h = v2(n + 1)
        x = n
        for j in range(h - 1):
            check("forced_climb", v2(3 * x + 1) == 1, (n, j))
            x = U(x)
        check("peak_and_block", x == peak(n) and U(x) == block(n), n)
        pred = (2 * n - 1) // 3 if n % 3 == 2 else None
        incoming = climb_closed(pred) if pred is not None else F()
        check("climb_resolvent", climb_closed(n) == nu(n) + incoming, n)
        # w_2(n)=nu(n)*4^-h pays each rise without a uniform discount.
        if h > 1:
            m = U(n)
            check("two_bit_payment", nu(n) / 4 ** h <= nu(m) / 4 ** (h - 1), n)
        else:
            m = U(n)
            if m > 1:
                ratio = climb_closed(n) / climb_closed(m)
                if ratio > 1 and first_fall_failure is None:
                    first_fall_failure = [n, m, str(ratio)]
                if ratio > max_ratio[0]:
                    max_ratio = (ratio, [n, m])
    # Sum every complete ascent from a finite source head, a second path.
    flow = Counter()
    for n in range(3, 1 << 12, 2):
        x = n
        for j in range(v2(n + 1)):
            flow[x] += nu(n)
            if j + 1 < v2(n + 1):
                x = U(x)
    for n in range(3, 1 << 12, 2):
        check("independent_climb_accumulation", flow[n] == climb_closed(n), n)
    check("finite_climb_total", sum(flow.values(), F())
          == sum((v2(n + 1) * nu(n) for n in range(3, 1 << 12, 2)), F()))
    refuels = []
    for s in range(1, 8):
        h = 2 * 3 ** s - 1
        n = (2 ** (h + 2) - 5) // 3
        m = 2 ** h - 1
        check("refuel_branch", n % 4 == 1 and U(n) == m and vp(n + 1, 3) == s, s)
        ratio = climb_closed(n) / climb_closed(m)
        check("unbounded_refuel_lower_bound", ratio > F(9, 64) * F(9, 4) ** s, s)
        refuels.append({"ternary_depth": s, "h": h, "source_bits": n.bit_length(),
                        "target_bits": m.bit_length(), "flow_ratio": str(ratio)})
    report["critical_flow"] = {
        "climb_closed_total_mass_on_n_gt_1": "2/3",
        "first_falling_edge_failure_in_head": first_fall_failure,
        "maximum_falling_edge_ratio_in_head": {"edge": max_ratio[1], "ratio": str(max_ratio[0])},
        "refuel_family": refuels,
        "status": "ascent inequalities closed exactly; full flow remains OPEN"}

    # Cycle controls for the critical-flow theorem and a fake budget proof.
    cycles = [[5, 7], [17, 25, 37, 55, 41, 61, 91]]
    outside = []
    for cycle in cycles:
        check("minus_cycle_control", [U(x, -1) for x in cycle] == cycle[1:] + cycle[:1])
        m = cycle[0]
        for a in range(1, 30):
            numerator = 2 ** a * m + 1
            if numerator % 3 == 0:
                p = numerator // 3
                if p % 2 and p > max(cycle):
                    check("outside_predecessor", U(p, -1) == m)
                    outside.append([p, m])
                    break
        else:
            raise RuntimeError("missing outside predecessor")
    for j in range(30):
        check("overdraft_false_proof", F(1, 2 ** (j + 2)) == F(1, 2) * F(1, 2 ** (j + 1)))
    report["hostile_controls"] = {
        "minus_cycles": cycles, "external_predecessors": outside,
        "budget_control": "identity base map, j -> j-1, kill j=0: contraction but no root receipt"}

    # The same valuation words at negative fixed anchors become expanding
    # positive macros. Exact source cylinders retain every intermediate guard.
    anchor_rows = []
    anchor_words = [(1, [1]), (5, [1, 2]), (17, [1, 1, 1, 2, 1, 1, 4])]
    for c, word in anchor_words:
        a = sum(word)
        r = len(word)
        p, denominator = 3 ** r, 2 ** a
        carry = c * (p - denominator)
        # Independent signed replay identifies the anchor itself.
        x = c
        for digit in word:
            check("anchor_minus_guard", v2(3 * x - 1) == digit, (c, x))
            x = U(x, -1)
        check("anchor_minus_return", x == c)
        for t in range(1, 65):
            n = 2 ** (a + 1) * t - c
            x = n
            for digit in word:
                check("anchor_plus_guard", v2(3 * x + 1) == digit, (c, n))
                x = U(x)
            check("anchor_affine", x == (p * n + carry) // denominator
                  and denominator * (x + c) == p * (n + c), (c, n))
            check("anchor_backwards", denominator * (x + c) // p - c == n)
        for n in range(3, 1 << 12, 2):
            expected = (v2(n + c) - 1) // a
            x, observed = n, 0
            while (x + c) % (2 ** (a + 1)) == 0:
                y = x
                for digit in word:
                    check("anchor_run_guard", v2(3 * y + 1) == digit)
                    y = U(y)
                x, observed = y, observed + 1
            check("anchor_repeat_depth", observed == expected, (c, n))
        for j in range(1, 21):
            k = a * j + 1
            residue = 2 ** k - c
            cylinder = nu(residue) + F(4, 3 * 4 ** k)
            check("anchor_repeat_mass", cylinder == F(1, 4 ** (a * j)), (c, j))
        expected_extra = F(1, 4 ** a - 1)
        check("anchor_resolvent_mass", sum((F(1, 4 ** (a * j)) for j in range(1, 21)), F())
              + F(1, 4 ** (20 * a)) * expected_extra == expected_extra)
        anchor_rows.append({"anchor": -c, "word": word, "P": p,
                            "Q": denominator, "B": carry,
                            "first_positive_source": 2 ** (a + 1) - c,
                            "first_endpoint": 2 * p - c,
                            "extra_resolvent_mass": str(expected_extra),
                            "total_resolvent_mass_on_n_gt_1": str(F(1, 3) + expected_extra)})
    # One ternary digit cannot authorize a macro consuming two or seven.
    check("anchor_ternary_depth_hostile", F(8 * (1 + 5), 9) - 5 == F(1, 3))
    check("anchor_ternary_depth_hostile", (2 ** 11 * (37 + 17)) % (3 ** 7) != 0)
    report["cycle_anchor_macros"] = anchor_rows

    # Explicit finite rooted universe, used only as a positive control.
    total_odd_steps = F()
    maximum = (0, 1)
    policy_counts = Counter()
    policy_lift, policy_endpoint, policy_source = Counter(), Counter(), Counter()
    stages = list(range(9))
    partial_mass = {r: F() for r in stages}
    frontier_mass = {r: F() for r in stages}
    for n in range(3, 1 << 14, 2):
        if (n + 17) % 4096 == 0:
            letter, word = "H", anchor_words[2][1]
        elif (n + 5) % 16 == 0:
            letter, word = "G", anchor_words[1][1]
        elif n % 4 == 3:
            letter, word = "R", [1]
        else:
            letter, word = "D", [v2(3 * n + 1)]
        policy_counts[letter] += 1
        policy_source[n] += nu(n)
        y = n
        for digit in word:
            policy_lift[y] += nu(n)
            check("complete_policy_guard", v2(3 * y + 1) == digit, (n, letter))
            y = U(y)
        if y != 1:
            policy_endpoint[y] += nu(n)
        check("complete_policy_direction", (y < n) if letter == "D" else (y > n), (n, letter))
        x, steps = n, 0
        while x != 1 and steps < 10000:
            x = U(x)
            steps += 1
        check("finite_root_control", x == 1, n)
        total_odd_steps += nu(n) * steps
        if steps > maximum[0]:
            maximum = steps, n
        x, cost, costs, frontiers = n, 0, [], []
        for r in stages:
            if x != 1:
                cost += v2(x + 1)
                x = block(x)
            costs.append(cost)
            frontiers.append(x != 1)
            partial_mass[r] += nu(n) * cost
            if x != 1:
                frontier_mass[r] += nu(n)
        check("block_cost_control", costs[-1] <= steps and (x != 1 or costs[-1] == steps), n)
    report["finite_root_control"] = {
        "universe": "odd 1 <= n < 2^14", "step_cap": 10000,
        "maximum_odd_steps": maximum[0], "maximum_source": maximum[1],
        "weighted_odd_steps_finite_head": str(total_odd_steps),
        "unexamined_source_mass": "1/24576",
        "scope": "no estimate of the stopping-time moment in the infinite tail"}
    report["complete_local_policy"] = {
        "finite_universe": "odd 3 <= n < 2^14", "counts": dict(sorted(policy_counts.items())),
        "scope": "universal legal block coverage; global root coverage remains OPEN",
        "exact_nu_domain_masses": {"D": "1/12", "G": "1/64", "H": "1/4194304",
                                    "R": str(F(1, 4) - F(1, 64) - F(1, 4194304))}}
    pushed_lift = Counter()
    for n, mass in policy_lift.items():
        m = U(n)
        if m != 1:
            pushed_lift[m] += mass
    for n in set(policy_lift) | set(pushed_lift) | set(policy_source) | set(policy_endpoint):
        check("policy_flow_lifting_identity", pushed_lift[n] == policy_lift[n]
              - policy_source[n] + policy_endpoint[n], n)
    report["finite_closure_stages"] = [
        {"refuels_paid": r, "mass_from_source_head": str(partial_mass[r]),
         "mass_tail_upper_bound": str((2 ** (r + 1) - 1) * F(32, 3 * 2 ** 14)),
         "unrooted_frontier_mass_from_head": str(frontier_mass[r]),
         "frontier_tail_upper_bound": "1/24576"}
        for r in stages]
    report["checks_by_group"] = dict(sorted(CHECKS.items()))
    report["checks_total"] = sum(CHECKS.values())
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
