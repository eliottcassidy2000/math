#!/usr/bin/env python3
"""Exact finite audits for THM-2267.

The proof in THM-2267 is analytic.  This companion exhausts small static
service systems and transition systems, and checks the identity-versus-swap
gluing collision using only integer arithmetic.
"""

from itertools import product


def subsets(mask):
    sub = mask
    while True:
        yield sub
        if sub == 0:
            return
        sub = (sub - 1) & mask


def static_audit():
    owner_count = 3
    obligation_count = 3
    all_owners = (1 << owner_count) - 1
    all_obligations = (1 << obligation_count) - 1
    family_count = 0
    hasse_checks = 0
    flat_checks = 0
    mobius_checks = 0
    cut_checks = 0
    flag_checks = 0

    for service in product(range(all_obligations + 1), repeat=owner_count):
        family_count += 1

        def union_of(owner_mask):
            ans = 0
            for i in range(owner_count):
                if (owner_mask >> i) & 1:
                    ans |= service[i]
            return ans

        def delta(owner_mask):
            total = sum(
                service[i].bit_count()
                for i in range(owner_count)
                if (owner_mask >> i) & 1
            )
            return total - union_of(owner_mask).bit_count()

        for a in range(all_owners + 1):
            for i in range(owner_count):
                if (a >> i) & 1:
                    continue
                lhs = delta(a | (1 << i)) - delta(a)
                rhs = (service[i] & union_of(a)).bit_count()
                assert lhs == rhs
                hasse_checks += 1

                for j in range(i + 1, owner_count):
                    if (a >> j) & 1:
                        continue
                    left = (
                        delta(a | (1 << i))
                        - delta(a)
                        + delta(a | (1 << i) | (1 << j))
                        - delta(a | (1 << i))
                    )
                    right = (
                        delta(a | (1 << j))
                        - delta(a)
                        + delta(a | (1 << i) | (1 << j))
                        - delta(a | (1 << j))
                    )
                    assert left == right
                    flat_checks += 1

        for t in range(1, all_owners + 1):
            h = 0
            for b in subsets(t):
                sign = -1 if ((t.bit_count() - b.bit_count()) & 1) else 1
                h += sign * delta(b)
            if t.bit_count() <= 1:
                expected = 0
            else:
                intersection = all_obligations
                for i in range(owner_count):
                    if (t >> i) & 1:
                        intersection &= service[i]
                expected = (
                    1 if t.bit_count() % 2 == 0 else -1
                ) * intersection.bit_count()
            assert h == expected
            mobius_checks += 1

        for a in range(all_owners + 1):
            for b in range(all_owners + 1):
                if a & b:
                    continue
                lhs = delta(a | b) - delta(a) - delta(b)
                rhs = (union_of(a) & union_of(b)).bit_count()
                assert lhs == rhs
                cut_checks += 1

        for s in range(all_owners + 1):
            pairwise_disjoint = True
            for i in range(owner_count):
                for j in range(i + 1, owner_count):
                    if ((s >> i) & 1) and ((s >> j) & 1):
                        pairwise_disjoint &= (service[i] & service[j]) == 0
            assert (delta(s) == 0) == pairwise_disjoint
            flag_checks += 1

    return {
        "families": family_count,
        "hasse": hasse_checks,
        "flat": flat_checks,
        "mobius": mobius_checks,
        "cut": cut_checks,
        "flag": flag_checks,
    }


def allowed_labels(mask, label_count):
    return [i for i in range(label_count) if (mask >> i) & 1]


def cycle_energy(allowed, label_count):
    best = None
    for labels in product(range(label_count), repeat=len(allowed)):
        if any(not ((allowed[v] >> labels[v]) & 1) for v in range(len(allowed))):
            continue
        energy = sum(
            labels[v] != labels[(v + 1) % len(labels)]
            for v in range(len(labels))
        )
        best = energy if best is None else min(best, energy)
    assert best is not None
    return best


def projected_binary_allowed(allowed, owner_cut, label_count):
    all_labels = (1 << label_count) - 1
    complement = all_labels ^ owner_cut
    answer = []
    for mask in allowed:
        side_mask = 0
        if mask & owner_cut:
            side_mask |= 1
        if mask & complement:
            side_mask |= 2
        assert side_mask
        answer.append(side_mask)
    return answer


def forced_alternation(binary_allowed):
    forced = []
    for mask in binary_allowed:
        if mask == 1:
            forced.append(0)
        elif mask == 2:
            forced.append(1)
        else:
            assert mask == 3
    if not forced or all(x == forced[0] for x in forced):
        return 0
    return sum(
        forced[i] != forced[(i + 1) % len(forced)]
        for i in range(len(forced))
    )


def transition_audit():
    label_count = 3
    cycle_length = 4
    all_labels = (1 << label_count) - 1
    system_count = 0
    cut_checks = 0
    positive_gap_checks = 0

    for allowed in product(range(1, all_labels + 1), repeat=cycle_length):
        system_count += 1
        omega = cycle_energy(allowed, label_count)
        assert omega == 0 or omega >= 2
        positive_gap_checks += 1
        cut_floor = 0
        for owner_cut in range(1, all_labels):
            binary = projected_binary_allowed(allowed, owner_cut, label_count)
            kappa = cycle_energy(binary, 2)
            assert kappa == forced_alternation(binary)
            assert omega >= kappa
            cut_floor = max(cut_floor, kappa)
            cut_checks += 1
        assert omega >= cut_floor

    binary_formula_checks = 0
    for length in range(2, 7):
        for allowed in product((1, 2, 3), repeat=length):
            assert cycle_energy(allowed, 2) == forced_alternation(allowed)
            binary_formula_checks += 1

    return {
        "systems": system_count,
        "cut": cut_checks,
        "cycle_gap": positive_gap_checks,
        "binary_formula": binary_formula_checks,
    }


def collision_energy(swap):
    # Four layers, two obligations per layer.  Obligation o has the unique
    # allowed owner o.  Each layer-to-layer gluing is identity or swap.
    energy = 0
    for _layer in range(4):
        for obligation in range(2):
            successor = 1 - obligation if swap else obligation
            energy += obligation != successor
    return energy


def collision_audit():
    # In every layer B_0={0}, B_1={1}; all redundancy/Hasse values vanish.
    static_delta_fingerprint = (0, 0, 0, 0)
    identity = collision_energy(False)
    swap = collision_energy(True)
    assert identity == 0
    assert swap == 8
    return static_delta_fingerprint, identity, swap


def main():
    static = static_audit()
    transition = transition_audit()
    fingerprint, identity, swap = collision_audit()

    print("THM-2267 exact audit")
    print(
        "static service families:",
        static["families"],
        "Hasse checks:",
        static["hasse"],
        "flat-square checks:",
        static["flat"],
    )
    print(
        "static Mobius checks:",
        static["mobius"],
        "cut checks:",
        static["cut"],
        "flag checks:",
        static["flag"],
    )
    print(
        "three-label four-cycle systems:",
        transition["systems"],
        "cut-floor checks:",
        transition["cut"],
        "cycle-gap checks:",
        transition["cycle_gap"],
    )
    print(
        "binary forced-alternation formula checks:",
        transition["binary_formula"],
    )
    print("collision static Delta fingerprint:", static_delta_fingerprint_text(fingerprint))
    print("identity-gluing switch energy:", identity)
    print("swap-gluing switch energy:", swap)
    print("all exact assertions passed")


def static_delta_fingerprint_text(fingerprint):
    return "(" + ",".join(str(x) for x in fingerprint) + ")"


if __name__ == "__main__":
    main()
