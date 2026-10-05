#!/usr/bin/env python3
"""Exact finite controls for a full-support atomic Collatz source measure.

No probabilistic sampling, third-party packages, or assumed Collatz termination.
The survivor tree concerns coefficient survival, not actual orbit non-descent.
Run from the repository root; stdout is the reproducible JSON artifact.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def weight(n):
    if n < 1 or n % 2 == 0:
        raise ValueError("positive odd source required")
    return Fraction(8, 3 * 4 ** n.bit_length())


def cylinder(r, depth):
    if not (depth >= 1 and 0 < r < 2**depth and r % 2):
        raise ValueError("canonical odd residue required")
    return weight(r) + Fraction(4, 3 * 4**depth)


def step(n, sign):
    return (3*n + sign)//2 if n % 2 else n//2


def measure_controls():
    for h in range(1, 14):
        head = sum((weight(n) for n in range(1, 2**h, 2)), Fraction())
        check(head == 1 - Fraction(2, 3*2**h), "exact head/tail split")
        # Finite complete prefix-code levels plus their explicitly summed tail.
        check(sum(Fraction(1, 2**j) for j in range(1, h+1))
              + Fraction(1, 2**h) == 1, "Kraft total")
    for depth in range(1, 11):
        total = Fraction()
        for r in range(1, 2**depth, 2):
            mass = cylinder(r, depth)
            # Independent mixture over uniform finite odd heads. For h>=depth,
            # the residue has probability 2**(1-depth) in each head.
            direct = sum((Fraction(1, 2**h) * Fraction(1, 2**(h-1))
                          for h in range(1, depth) if r < 2**h), Fraction())
            direct += Fraction(4, 4**depth)
            check(mass == direct, "cylinder from independent head mixture")
            left = cylinder(r, depth+1)
            right = cylinder(r+2**depth, depth+1)
            check(left+right == mass, "cylinder additivity")
            check(right == Fraction(1, 4**depth), "new high-one child")
            check(right <= mass/4, "new high-one conditional cost")
            total += mass
        check(total == 1, "partition mass")
    sign_rows = []
    for depth in range(1, 33):
        plus = cylinder(2**depth-1, depth)
        minus = cylinder(1, depth)
        check(plus == Fraction(4, 4**depth), "positive all-odd prefix")
        check(minus == weight(1)+Fraction(4, 3*4**depth),
              "negative all-odd prefix")
        sign_rows.append({"depth": depth, "plus": str(plus), "minus": str(minus)})
    return {"all_odd_sign_controls": sign_rows,
            "minus_cycle_minimum_atom_mass": str(weight(1)+weight(5)+weight(17))}


def survives(n, sign, depth):
    p = 1
    for t in range(1, depth+1):
        if n % 2:
            p *= 3
        n = step(n, sign)
        if p < 2**t:
            return False
    return True


def atomic_histogram_mass(hist, depth):
    count = sum(hist.values())
    numerator = 8*sum(c * 4**(depth-ell) for ell, c in hist.items()) + 4*count
    return Fraction(numerator, 3*4**depth)


def survivor_tree(max_depth=26):
    # (canonical residue r, actual endpoint s, coefficient numerator p).
    # A word maps r + 2**depth*k to s + p*k.
    states = [(0, 0, 1)]
    rows = []
    accepted_count = 0
    fringe = set()
    for depth in range(1, max_depth+1):
        old_q, q = 2**(depth-1), 2**depth
        children = []
        for r, s, p in states:
            for bit in (0, 1):
                if depth == 1 and bit == 0:
                    continue
                k = (bit-s) % 2
                rr = r + old_q*k
                pp = p * (3 if bit else 1)
                ss = (3*(s+p*k)+1)//2 if bit else (s+p*k)//2
                if pp >= q:
                    children.append((rr, ss, pp))
                else:
                    accepted_count += 1
                    threshold = max(0, (ss-rr)//(q-pp)+1)
                    # This is a finite observation, not an all-depth theorem.
                    check(threshold <= 1, "observed one-member plus fringe")
                    for j in range(threshold):
                        fringe.add(rr+q*j)
        states = children
        plus = Counter(r.bit_length() for r, _, _ in states)
        minus = Counter((q-r).bit_length() for r, _, _ in states)
        plus_mass = atomic_histogram_mass(plus, depth)
        minus_mass = atomic_histogram_mass(minus, depth)
        check(minus_mass >= weight(1)+weight(5)+weight(17),
              "negative cycle atom lower control")
        check(len({r for r, _, _ in states}) == len(states), "unique word cylinders")
        if depth <= 14:
            for sign, expected in ((1, {r for r, _, _ in states}),
                                   (-1, {q-r for r, _, _ in states})):
                native = {n for n in range(1, q, 2) if survives(n, sign, depth)}
                check(native == expected, "independent direct parity tree")
        small_plus = sorted((r for r, _, _ in states))[:12]
        small_minus = sorted((q-r for r, _, _ in states))[:12]
        rows.append({"depth": depth, "cylinders": len(states),
                     "haar_conditional_odd": str(Fraction(len(states), q//2)),
                     "plus_atomic": str(plus_mass), "minus_atomic": str(minus_mass),
                     "plus_atomic_decimal": float(plus_mass),
                     "minus_atomic_decimal": float(minus_mass),
                     "plus_height_histogram": dict(sorted(plus.items())),
                     "minus_height_histogram": dict(sorted(minus.items())),
                     "least_plus_sources": small_plus, "least_minus_sources": small_minus,
                     "accepted_ports_through_depth": accepted_count,
                     "actual_plus_fringe_through_depth": sorted(fringe)})
    return rows


def first_coefficient_exit(n, sign, cap):
    p = 1
    for t in range(1, cap+1):
        if n % 2:
            p *= 3
        n = step(n, sign)
        if p < 2**t:
            return t
    return cap+1


def finite_head(head_bits=17, cap=512):
    denominator = 3*4**head_bits
    tail = Fraction(2, 3*2**head_bits)
    rows = {}
    for sign, roots in ((1, {1}), (-1, {1, 5, 17})):
        cache = {n: 0 for n in roots}
        tau_hist = Counter()
        unknown_units = 0
        maximum = (0, 1)
        for n in range(1, 2**head_bits, 2):
            chain, seen, x = [], set(), n
            while x not in cache and x not in seen and len(chain) <= 10000:
                seen.add(x)
                chain.append(x)
                x = step(x, sign)
            units = 8*4**(head_bits-n.bit_length())
            if x not in cache:
                unknown_units += units
                continue
            tau = cache[x]
            for v in reversed(chain):
                tau += 1
                cache[v] = tau
            check(cache[n] >= 0, "finite rooted path length")
            tau_hist[cache[n]] += units
            maximum = max(maximum, (cache[n], n))
        check(sum(tau_hist.values()) + unknown_units == denominator*(1-tail),
              "finite source weights conserved")
        for n in list(range(1, 1024, 2)) + [maximum[1]]:
            x, direct_tau = n, 0
            while x not in roots and direct_tau <= 10000:
                x = step(x, sign)
                direct_tau += 1
            check(x in roots and direct_tau == cache[n], "root clock independent replay")
        intervals = []
        for t in (0, 1, 2, 4, 8, 16, 24, 32, 64, 128, 256, 512):
            lo = Fraction(sum(u for tau, u in tau_hist.items() if tau > t), denominator)
            hi = lo + Fraction(unknown_units, denominator) + tail
            intervals.append({"time": t, "unrooted_lower": str(lo),
                              "unrooted_upper": str(hi),
                              "lower_decimal": float(lo), "upper_decimal": float(hi)})
        rows[str(sign)] = {"targets": sorted(roots), "unknown_head_weight": str(Fraction(unknown_units, denominator)),
                           "max_known_head_time_and_source": maximum, "intervals": intervals}
    controls = []
    for n in (1, 3, 7, 27, 703, 1819, 6263, 1000001):
        controls.append({"n": n, "weight": str(weight(n)),
                         "plus_coefficient_exit": first_coefficient_exit(n, 1, cap)})
    for n in (1, 5, 17):
        check(first_coefficient_exit(n, -1, cap) == cap+1,
              "negative known cycle coefficient survival")
    return {"head_bits": head_bits, "odd_sources": 2**(head_bits-1),
            "explicit_tail": str(tail), "root_survival_intervals": rows,
            "coefficient_exit_controls": controls}


def partition_controls():
    out = []
    for depth in range(1, 13):
        cells = [cylinder(r, depth) for r in range(1, 2**depth, 2)]
        shell_square = sum((1 if ell == 1 else 2**(ell-2))
                           * (Fraction(8, 3*4**ell)+Fraction(4, 3*4**depth))**2
                           for ell in range(1, depth+1))
        check(sum(cells) == 1, "partition q=1")
        check(sum(x*x for x in cells) == shell_square, "q=2 independent shell sum")
        check(len(cells) == 2**(depth-1), "partition q=0")
        out.append({"depth": depth, "Z_0": len(cells), "Z_1": "1", "Z_2": str(shell_square)})
    return out


def odd_step(n):
    m = 3*n+1
    return m//(m & -m)


def sibling_controls():
    rows = []
    for k in range(4, 25):
        merges = 0
        for t in range(1, 256, 2):
            x = t*2**k-1
            small = 9*t*2**(k-3)-1
            big = odd_step(odd_step(x))
            check(big == 2*small+1, "initial sibling pair")
            for _ in range(k-4):
                big, small = odd_step(big), odd_step(small)
                check(big == 2*small+1, "forced sibling P1 chain")
            v = 3*small+1
            b = (v & -v).bit_length()-1
            check(b >= 2, "first non-one valuation")
            local_merge = b == 2
            check(local_merge == (t % 4 == pow(3, k, 4)), "one-bit merge guard")
            big, small = odd_step(big), odd_step(small)
            if local_merge:
                merges += 1
                check(big == 4*small+1, "quarter-child reached at k-3")
                check(odd_step(big) == odd_step(small), "coalescence one step later")
        check(merges == 64, "complete t mod 256 Haar control")
        cell_mass = cylinder(2**k-1, k+1)
        r = pow(3, k, 4)*2**k-1
        merge_mass = cylinder(r, k+2)
        conditional = merge_mass/cell_mass
        check(conditional == Fraction(11 if k % 2 == 0 else 1, 12),
              "atomic sibling merge probability")
        rows.append({"k": k, "merge_t_mod_4": pow(3, k, 4),
                     "haar_conditional": "1/2", "atomic_conditional": str(conditional),
                     "quarter_child_state_time": k-3, "actual_coalescence_time": k-2})
    # Two small sheet-commutator arithmetic controls.
    check(odd_step(odd_step(7)) == 17, "hostile original side")
    check(odd_step(21) == 1, "hostile multiplier side")
    # A deep-cell failed structural merge with an independently verified common
    # endpoint above 1: x=31, shadow=35; both reach 53 (at different times).
    check(9*2**(5-3)-1 == 35, "deep shadow at x31")
    original, shadow = 31, 35
    original_path = [original]
    for _ in range(40):
        original = odd_step(original)
        original_path.append(original)
    check(53 in original_path and odd_step(shadow) == 53,
          "automaton break does not exclude asynchronous common future")
    return {"deep_cells": rows, "break_but_common_future": {"x": 31, "shadow": 35, "endpoint": 53}}


def transfer_hostiles():
    rows = []
    for t in range(1, 33):
        source, target = 2**(t+1)-1, 2*3**t-1
        x = source
        for _ in range(t):
            check(x % 4 == 3, "all valuation-one climb")
            x = odd_step(x)
        check(x == target, "long climb endpoint")
        ratio = weight(source)/weight(target)
        check(ratio == 4**(target.bit_length()-source.bit_length()), "incoming weight ratio")
        check(ratio >= 1, "no strict uniform killed-transfer contraction")
        rows.append({"odd_steps": t, "source": source, "target": target,
                     "single_preimage_weight_ratio": str(ratio)})
    return rows


def climb_credit_controls():
    def height(n):
        return ((n+1) & -(n+1)).bit_length()-1

    def adjusted(n):
        return weight(n)/8**height(n)

    for n in range(3, 32768, 4):
        child = odd_step(n)
        check(height(child) == height(n)-1, "climb credit clock")
        check(adjusted(n) <= adjusted(child)/2, "all valuation-one edges paid")
    refuels = []
    for h in range(3, 64, 2):
        n, target = (2**(h+2)-5)//3, 2**h-1
        check(3*n+1 == 4*target and odd_step(n) == target, "valuation-two refuel")
        check(height(n) == 1 and height(target) == h, "new climb fuel")
        ratio = adjusted(n)/adjusted(target)
        check(ratio == Fraction(8**(h-1), 4), "unbounded refuel weight ratio")
        refuels.append({"h": h, "source": n, "target": target,
                        "adjusted_weight_ratio": str(ratio)})
    return {"valuation_one_sources_checked": 8192, "refuels": refuels}


def finite_observer_controls():
    # A finite sigma-algebra obstruction, with no nonmeasurable subset involved.
    pairs = list(combinations(range(7), 2))
    internal = [i for i, (a, b) in enumerate(pairs) if b < 4]
    outside = [i for i in range(21) if i not in internal]
    reversal = sum(1 << i for i in internal)

    def arcs(mask):
        out = [0]*7
        for i, (a, b) in enumerate(pairs):
            u, v = (a, b) if mask >> i & 1 else (b, a)
            out[u] |= 1 << v
        return out

    def lam(out):
        cover = Counter()
        for a, b, c in combinations(range(7), 3):
            subset = (1 << a) | (1 << b) | (1 << c)
            if all((out[v] & subset).bit_count() == 1 for v in (a, b, c)):
                cover.update(((a, b), (a, c), (b, c)))
        return tuple(cover[p] for p in pairs)

    def hp(out):
        dp = [[0]*7 for _ in range(128)]
        for v in range(7):
            dp[1 << v][v] = 1
        for mask in range(1, 128):
            for v in range(7):
                if not dp[mask][v]:
                    continue
                nxt = out[v] & (127 ^ mask)
                while nxt:
                    bit = nxt & -nxt
                    nxt -= bit
                    dp[mask | bit][bit.bit_length()-1] += dp[mask][v]
        return sum(dp[127])

    fixed = next(sum(((bits >> j) & 1) << idx for j, idx in enumerate(internal))
                 for bits in range(64)
                 if sorted((arcs(sum(((bits >> j) & 1) << idx for j, idx in enumerate(internal)))[v] & 15).bit_count()
                           for v in range(4)) == [1, 1, 2, 2])
    witness = None
    inspected = 0
    for bits in range(1 << len(outside)):
        inspected += 1
        mask = fixed | sum(((bits >> j) & 1) << idx for j, idx in enumerate(outside))
        a, b = arcs(mask), arcs(mask ^ reversal)
        la, lb = lam(a), lam(b)
        if la != lb:
            continue
        ha, hb = hp(a), hp(b)
        if ha == hb:
            continue
        for out, h in ((a, ha), (b, hb)):
            direct = sum(all(out[p[i]] >> p[i+1] & 1 for i in range(6))
                         for p in permutations(range(7)))
            check(h == direct, "Hamiltonian count DP versus permutations")
        check(la == lb and ha != hb, "finite labelled-observer collision")
        witness = {"mask": mask, "reversed_mask": mask ^ reversal,
                   "lambda_vector": la, "hamiltonian_paths": [ha, hb],
                   "S": [0, 1, 2, 3], "pair_order": pairs,
                   "bit_one_means": "lower vertex -> higher vertex"}
        break
    check(witness is not None, "finite observer witness found")
    cyclic = {(0, 1), (1, 2), (2, 0)}
    paths = [p for p in permutations(range(3)) if all((p[i], p[i+1]) in cyclic for i in range(2))]
    check(len(paths) == 3, "directed triangle has three Hamiltonian paths")
    # AP lonely-set witness control, entirely rational.
    ap = []
    for n in range(2, 25):
        from math import gcd
        actual = [a for a in range(n)
                  if all(min((v*a) % n, (-v*a) % n) >= 1 for v in range(1, n))]
        units = [a for a in range(1, n) if gcd(a, n) == 1]
        check(actual == units, "primitive n-gon AP witnesses")
        ap.append({"n": n, "primitive_numerators": actual})
    return {"tournament_witness": witness, "candidate_masks_inspected": inspected,
            "cyclic_triangle_paths": paths, "AP_grid_controls": ap}


def main():
    result = {"measure_controls": measure_controls(), "survivor_tree": survivor_tree(),
              "finite_head": finite_head(), "partition_controls": partition_controls(),
              "sibling_controls": sibling_controls(), "finite_observer_controls": finite_observer_controls(),
              "transfer_hostiles": transfer_hostiles(), "climb_credit_controls": climb_credit_controls()}
    result["checks"] = CHECKS
    result["status"] = "FINITE-EXACT controls passed; universal Collatz OPEN"
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
