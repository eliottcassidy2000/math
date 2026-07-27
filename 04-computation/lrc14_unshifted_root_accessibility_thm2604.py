#!/usr/bin/env python3
"""Exact referee for THM-2604's unshifted-root accessibility boundary.

For I_h=(h/13,(h+1)/13) and

    D_{k,L}={y mod 1 : ||k y|| < L/14},

the script computes the positive-measure mask

    M_{k,L}(h)=1 iff measure(I_h intersect D_{k,L})>0

by two independent exact routes: clipped danger intervals and the cleared
integer inequalities for a tooth centre j/k.  It also checks the sharp
ordinary/guard thresholds, the shifted-atlas hostile, all owner-pivot label
choices, and finite interval witnesses for the elementary cross-mixing lemma.

Everything is dependency-free and exact over ``Fraction``.  This script does
not enumerate live scalar covers and does not claim that the typed small-speed
control is one.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


P = 13


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def merge_intervals(intervals):
    """Merge a finite list of rational intervals; endpoints are measure-null."""
    ordered = sorted((a, b) for a, b in intervals if a < b)
    out = []
    for a, b in ordered:
        if not out or out[-1][1] < a:
            out.append([a, b])
        elif b > out[-1][1]:
            out[-1][1] = b
    return [(a, b) for a, b in out]


def intersect_intervals(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


@lru_cache(maxsize=None)
def danger_intervals(k, L):
    """D_{k,L} on [0,1], clipped from the centres j/k, 0<=j<=k."""
    radius = Fraction(L, 14 * k)
    intervals = []
    for j in range(k + 1):
        centre = Fraction(j, k)
        intervals.append((max(Fraction(0), centre - radius),
                          min(Fraction(1), centre + radius)))
    return merge_intervals(intervals)


def digit_interval(h):
    return [(Fraction(h, P), Fraction(h + 1, P))]


@lru_cache(maxsize=None)
def access_by_intervals(k, L, h):
    return bool(intersect_intervals(danger_intervals(k, L), digit_interval(h)))


def access_by_integer_inequalities(k, L, h):
    """Existence of j/k whose open tooth meets the open digit cell."""
    return any(
        182 * j - 13 * L < 14 * k * (h + 1)
        and 182 * j + 13 * L > 14 * k * h
        for j in range(k + 1)
    )


def translated_access_by_integer_inequalities(k, L, h, s):
    """Access for d_L(k y+s/13), with target translate s in F_13."""
    return any(
        182 * j - 14 * s - 13 * L < 14 * k * (h + 1)
        and 182 * j - 14 * s + 13 * L > 14 * k * h
        for j in range(k + 2)
    )


@lru_cache(maxsize=None)
def translated_danger_intervals(k, L, s):
    """D^s_{k,L} on [0,1], independently clipped from translated teeth."""
    radius = Fraction(L, 14 * k)
    intervals = []
    # If 0<=s<=12 and 0<=y<=1, only integer tooth labels 0,...,k+1
    # can lie within distance L/14 of ky+s/13.  The k+1 tooth is needed,
    # e.g. at (k,L,s)=(1,2,12).
    for j in range(k + 2):
        centre = (Fraction(j) - Fraction(s, P)) / k
        intervals.append((max(Fraction(0), centre - radius),
                          min(Fraction(1), centre + radius)))
    return merge_intervals(intervals)


def translated_access_by_intervals(k, L, h, s):
    return bool(intersect_intervals(
        translated_danger_intervals(k, L, s), digit_interval(h)))


@lru_cache(maxsize=None)
def access_mask(k, L):
    return tuple(h for h in range(P) if access_by_intervals(k, L, h))


def complement_intervals(intervals):
    out = []
    cursor = Fraction(0)
    for a, b in intervals:
        if cursor < a:
            out.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def shifted_boundary_points(k, L):
    """Return the full labelled target boundary atlas, separated by sign."""
    by_sign = {}
    for eps in (-1, 1):
        entries = []
        for j in range(k):
            for s in range(P):
                y = (Fraction(j) + Fraction(eps * L, 14)
                     - Fraction(s, P)) / k
                y %= 1
                entries.append((y, j, s))
        by_sign[eps] = entries
    return by_sign


def v7(n):
    out = 0
    while n % 7 == 0:
        out += 1
        n //= 7
    return out


def vp(n, prime):
    out = 0
    while n % prime == 0:
        out += 1
        n //= prime
    return out


def circle_distance(q):
    q %= 1
    return min(q, 1 - q)


def gate_root_mask_on_z_interval(k, L, z0, z1):
    """Roots r for which d_L(k(z+r)/13)=1 throughout z0<z<z1."""
    require(z0 < z1, "empty z interval")
    midpoint = (z0 + z1) / 2
    mask = []
    for r in range(P):
        # A truth value can change only at k(z+r)/13=n +/- L/14.
        qa = Fraction(k, P) * (z0 + r)
        qb = Fraction(k, P) * (z1 + r)
        for sign in (-1, 1):
            low = qa - Fraction(sign * L, 14)
            high = qb - Fraction(sign * L, 14)
            candidate = low.numerator // low.denominator + 1
            require(not (candidate < high),
                    f"gate boundary inside claimed constant cell: {(k,L,r,sign,candidate)}")
        if circle_distance(Fraction(k, P) * (midpoint + r)) < Fraction(L, 14):
            mask.append(r)
    return tuple(mask)


def canonical_wall(mask, tau):
    """THM-2531 lexicographic necklace marker and first 1-to-0 wall."""
    bits = set(mask)
    require(bits and len(bits) < P and tau % P, "wall needs a mixed mask and nonzero slope")
    words = {
        a: tuple(int((a + j * tau) % P in bits) for j in range(P))
        for a in range(P)
    }
    alpha = max(words, key=words.get)
    require(words[alpha][0] == 1, "lexicographic marker is not occupied")
    q = next(j for j in range(1, P) if (alpha + j * tau) % P not in bits)
    source = (alpha + (q - 1) * tau) % P
    head = (source + tau) % P
    return alpha, q, source, head


def safe_root_count(z, roles):
    """Number of direct roots safe for every (speed,width) role at base z."""
    return sum(
        all(circle_distance(Fraction(k, P) * (z + r)) > Fraction(L, 14)
            for k, L in roles)
        for r in range(P)
    )


def root_occupancy_wall_profile(roles):
    """Exact open-cell histogram of n(z)=#{safe direct roots}."""
    walls = {Fraction(0), Fraction(1)}
    for k, L in roles:
        radius = Fraction(L, 14)
        for r in range(P):
            # Here 0<=k(z+r)/13<=k, so this finite m range is exhaustive.
            for m in range(-1, k + 2):
                for sign in (-1, 1):
                    z = Fraction(P, k) * (m + sign * radius) - r
                    if 0 <= z <= 1:
                        walls.add(z)
    walls = sorted(walls)
    histogram = {}
    for a, b in zip(walls, walls[1:]):
        require(a < b, "duplicate occupancy wall survived")
        n = safe_root_count((a + b) / 2, roles)
        histogram[n] = histogram.get(n, Fraction(0)) + (b - a)
    require(sum(histogram.values(), Fraction(0)) == 1,
            "occupancy cells do not partition the circle")
    image_mass = sum(length for n, length in histogram.items() if n)
    physical_mass = sum(n * length for n, length in histogram.items()) / P
    return walls, histogram, image_mass, physical_mass


def role_walls(k, L):
    walls = {Fraction(0), Fraction(1)}
    radius = Fraction(L, 14)
    for r in range(P):
        for m in range(-1, k + 2):
            for sign in (-1, 1):
                z = Fraction(P, k) * (m + sign * radius) - r
                if 0 <= z <= 1:
                    walls.add(z)
    return walls


def universal_small_role_cells(role_universe):
    """One exact wall refinement and root bitsets for many finite profiles."""
    walls = sorted(set().union(*(role_walls(*role) for role in role_universe)))
    cells = tuple(zip(walls, walls[1:]))
    mids = tuple((a + b) / 2 for a, b in cells)
    lengths = tuple(b - a for a, b in cells)
    all_cells = (1 << len(cells)) - 1
    bitsets = {}
    for role in role_universe:
        k, L = role
        roots = []
        for r in range(P):
            bits = 0
            for index, z in enumerate(mids):
                if circle_distance(Fraction(k, P) * (z + r)) > Fraction(L, 14):
                    bits |= 1 << index
            roots.append(bits)
        bitsets[role] = tuple(roots)
    require(sum(lengths, Fraction(0)) == 1, "universal cells do not partition")
    return all_cells, lengths, bitsets


def simultaneous_root_bits(roles, all_cells, bitsets):
    return tuple(
        _intersect_role_root(roles, r, all_cells, bitsets)
        for r in range(P)
    )


def _intersect_role_root(roles, root, all_cells, bitsets):
    out = all_cells
    for role in roles:
        out &= bitsets[role][root]
    return out


def at_least_root_bits(root_bits, threshold, all_cells):
    """Cells on which at least ``threshold`` of the root bits are present."""
    dp = [all_cells] + [0] * threshold
    for bits in root_bits:
        for level in range(threshold, 0, -1):
            dp[level] |= dp[level - 1] & bits
    return dp[threshold]


def bit_measure(bits, lengths):
    return sum(length for index, length in enumerate(lengths) if bits >> index & 1)


def pivot_choices(maximal, big_roles=None):
    """Triples (q_star,k_a,u_0) with both distinct graft roles retained."""
    roles = set(range(6))
    allowed_k = roles if big_roles is None else roles & set(big_roles)
    return [
        (qstar, ka, u0)
        for qstar in maximal
        for ka in allowed_k
        if ka != qstar
        for u0 in range(6)
        if u0 != qstar and u0 != ka
    ]


def central_digit_interval(h):
    return (Fraction(h, P) + Fraction(1, 52),
            Fraction(h + 1, P) - Fraction(1, 52))


def grid_preimage_in_interval(y, N, interval):
    """Find x=(r+y)/13^N strictly in interval, or return None."""
    m = P**N
    a, b = interval
    lower = m * a - y
    r = lower.numerator // lower.denominator + 1
    if 0 <= r < m:
        x = (r + y) / m
        if a < x < b:
            return x
    return None


def main():
    print("== THM-2604: unshifted future-root accessibility ==")

    # Independent exhaustive agreement and the exact safe-gap geometry.
    compared = 0
    for L in (1, 2):
        for k in range(1, 201):
            intervals = danger_intervals(k, L)
            gaps = complement_intervals(intervals)
            require(gaps, f"missing safe gap at k={k}, L={L}")
            require(max(b - a for a, b in gaps) == Fraction(7 - L, 7 * k),
                    f"safe-gap length mismatch at k={k}, L={L}")
            for h in range(P):
                a = access_by_intervals(k, L, h)
                b = access_by_integer_inequalities(k, L, h)
                require(a == b, f"access routes disagree at k={k}, L={L}, h={h}")
                compared += 1
    print(f"interval-union/integer-inequality agreement: {compared} cells PASS")
    print("exact tooth safe-gap length: (7-L)/(7k) PASS")

    expected_L1 = {
        1: (0, 12),
        2: (0, 6, 12),
        3: (0, 4, 8, 12),
        4: (0, 3, 6, 9, 12),
        5: (0, 2, 5, 7, 10, 12),
        6: (0, 2, 4, 6, 8, 10, 12),
        7: (0, 1, 3, 5, 7, 9, 11, 12),
        8: (0, 1, 3, 4, 6, 8, 9, 11, 12),
        9: (0, 1, 2, 4, 5, 7, 8, 10, 11, 12),
        10: (0, 1, 2, 3, 5, 6, 7, 9, 10, 11, 12),
        11: (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
    }
    expected_L2 = {
        1: (0, 1, 11, 12),
        2: (0, 5, 6, 7, 12),
        3: (0, 3, 4, 8, 9, 12),
        4: (0, 2, 3, 6, 9, 10, 12),
        5: (0, 2, 4, 5, 7, 8, 10, 12),
        6: (0, 1, 2, 4, 6, 8, 10, 11, 12),
        7: (0, 1, 2, 3, 5, 7, 9, 10, 11, 12),
        8: (0, 1, 3, 4, 5, 6, 7, 8, 9, 11, 12),
        9: (0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12),
    }
    for k, expected in expected_L1.items():
        require(access_mask(k, 1) == expected, f"L=1 table mismatch at k={k}")
    for k, expected in expected_L2.items():
        require(access_mask(k, 2) == expected, f"L=2 table mismatch at k={k}")
    require(all(access_mask(k, 1) == tuple(range(P)) for k in range(12, 201)),
            "ordinary all-root threshold failed")
    require(all(access_mask(k, 2) == tuple(range(P)) for k in range(10, 201)),
            "guard all-root threshold failed")
    require(all(access_mask(k, 1) != tuple(range(P)) for k in range(1, 12)),
            "ordinary threshold not sharp")
    require(all(access_mask(k, 2) != tuple(range(P)) for k in range(1, 10)),
            "guard threshold not sharp")
    print("L=1 masks k=1..11:")
    for k in range(1, 12):
        print(f"  k={k:2d}: {list(access_mask(k, 1))}")
    print("L=2 masks k=1..9:")
    for k in range(1, 10):
        print(f"  k={k:2d}: {list(access_mask(k, 2))}")
    print("all-root thresholds: ordinary L=1 iff k>=12; guard L=2 iff k>=10 PASS")
    translated_checks = 0
    for L, threshold in ((1, 12), (2, 10)):
        for k in range(threshold, 201):
            for s in range(P):
                for h in range(P):
                    by_inequalities = translated_access_by_integer_inequalities(k, L, h, s)
                    by_intervals = translated_access_by_intervals(k, L, h, s)
                    require(by_inequalities == by_intervals,
                            f"translated routes disagree at {(k,L,h,s)}")
                    require(by_inequalities,
                            f"translated all-root access failed at {(k,L,h,s)}")
                    translated_checks += 1
    require(translated_access_by_integer_inequalities(1, 2, 12, 12),
            "translated k+1 endpoint tooth was omitted")
    require(translated_access_by_intervals(1, 2, 12, 12),
            "translated interval route omitted the k+1 endpoint tooth")
    print(f"uniform translated access above the thresholds: {translated_checks} (k,h,s) cells PASS")
    require(not access_by_intervals(11, 1, 6), "k=11 ordinary hostile disappeared")
    require(not access_by_intervals(9, 2, 6), "k=9 guard hostile disappeared")
    print("sharp hostiles: (k,L,h)=(11,1,6) and (9,2,6) are inaccessible PASS")

    # Shifted boundaries still pierce the hostile digit: this separates the
    # complete target atlas from the unshifted q=0 factor.
    for k, L in ((11, 1), (9, 2)):
        atlas = shifted_boundary_points(k, L)
        all_points = []
        for eps, entries in atlas.items():
            points = [y for y, _, _ in entries]
            require(len(set(points)) == P * k,
                    f"fixed-sign shifted atlas collision at k={k}, L={L}, eps={eps}")
            for h in range(P):
                count = sum(Fraction(h, P) < y < Fraction(h + 1, P) for y in points)
                require(count == k, f"wrong fixed-sign digit count at k={k},L={L},h={h}")
            all_points.extend(points)
        require(len(set(all_points)) == 2 * P * k,
                f"cross-sign shifted atlas collision at k={k}, L={L}")
        hostile_count = sum(Fraction(6, P) < y < Fraction(7, P) for y in all_points)
        require(hostile_count == 2 * k, "shifted hostile-cell count mismatch")
        print(f"shifted atlas at (k,L)=({k},{L}) has {hostile_count} boundaries in I_6 PASS")

    # Exhaust the abstract six-label pivot choice.  Label 0 is the guard;
    # labels 1..5 are ordinary.  The maximal-v_7 set is arbitrary nonempty.
    # A role is ``big`` when its unshifted mask is all-root: H>=10 for label
    # 0 (width two), or q_i>=12 for labels 1..5 (width one).
    cases = successful = fixed_qstar_checks = 0
    for max_bits in range(1, 1 << 6):
        maximal = {i for i in range(6) if max_bits >> i & 1}
        for big_bits in range(1 << 6):
            big = {i for i in range(6) if big_bits >> i & 1}
            cases += 1
            brute = bool(pivot_choices(maximal, big))
            criterion = bool(big) and not (len(big) == 1 and maximal == big)
            require(brute == criterion,
                    f"free-qstar criterion mismatch: maximal={maximal}, big={big}")
            successful += int(brute)
            for qstar in maximal:
                require(pivot_choices({qstar}),
                        f"no second graft choice for qstar={qstar}")
                fixed_big = bool(big - {qstar})
                brute_fixed = bool(pivot_choices({qstar}, big))
                require(brute_fixed == fixed_big,
                        f"fixed-qstar criterion mismatch: qstar={qstar}, big={big}")
                fixed_qstar_checks += 1
    require(cases == 63 * 64, "pivot case count mismatch")
    print(f"pivot combinatorics: {cases} maximal/big-set cases, {successful} admit an all-root k_a PASS")
    print(f"fixed-q_star tests: {fixed_qstar_checks}; a distinct k_a and u_0 always exist PASS")
    print("criterion: B_big nonempty, except B_big={b}=Max_7 with b the forced q_star")

    # A typed pivot/scalarization control, explicitly not asserted to cover.
    H = 1
    qs = (4, 2, 3, 6, 10)
    speeds = (H,) + qs
    require(all(v > 0 and gcd(v, P) == 1 for v in speeds), "typed tuple is not 13-unit")
    require(len(set(speeds)) == 6, "typed tuple is not pairwise distinct")
    require(H % 2 == 1, "guard is not odd")
    require(set(i for i, v in enumerate(speeds) if v7(v) == max(map(v7, speeds)))
            == set(range(6)), "typed tuple maximal-v7 set mismatch")
    require(H < 10 and not any(q >= 12 for q in qs),
            "typed tuple unexpectedly has an all-root role")
    print("typed non-cover control: H=1, q=(4,2,3,6,10), all six nu_7=0, no all-root role PASS")

    # Fixed-canonical-head hostile.  These are literal direct-fibre masks on
    # one open z-cell, not an assertion that the displayed tuple covers the
    # circle globally.  The strict owner interval is nested inside the cell.
    z0, z1 = Fraction(13, 56), Fraction(15, 56)
    hostile_roles = ((3, 2), (1, 1), (2, 1), (4, 1), (5, 1), (7, 1))
    hostile_masks = {
        role: gate_root_mask_on_z_interval(role[0], role[1], z0, z1)
        for role in hostile_roles
    }
    expected_hostile_masks = {
        (3, 2): (0, 4, 8, 9),
        (1, 1): (0, 12),
        (2, 1): (0, 6),
        (4, 1): (3,),
        (5, 1): (5, 10),
        (7, 1): (9, 11),
    }
    require(hostile_masks == expected_hostile_masks,
            f"fixed-head hostile masks changed: {hostile_masks}")
    occupied = tuple(sorted(set(range(P)) - set().union(*map(set, hostile_masks.values()))))
    require(occupied == (1, 2, 7), "fixed-head hostile safe mask changed")
    tau_H = (-pow(3, -1, P)) % P
    require(tau_H == 4, "guard slope convention changed")
    alpha, run, source, head = canonical_wall(occupied, tau_H)
    require((alpha, run, source, head) == (7, 1, 7, 11),
            "fixed-head canonical wall changed")
    eligible_roles = ((3, 2), (1, 1), (2, 1), (4, 1), (5, 1))
    eligible_access = set().union(*(set(access_mask(k, L)) for k, L in eligible_roles))
    require(eligible_access == set(range(P)) - {1, 11},
            "fixed-head future accessibility union changed")
    require(head not in eligible_access, "hostile head became future-accessible")
    current_failures_at_head = [role for role in hostile_roles if head in hostile_masks[role]]
    require(current_failures_at_head == [(7, 1)], "q_star is not the sole current head failure")

    blockers = (52, 169, 13**5)
    eps = Fraction(1, 100 * 13**4)
    owner0, owner1 = Fraction(1, 4) - eps, Fraction(1, 4) + eps
    require(z0 < owner0 < owner1 < z1, "owner interval is not nested in hostile cell")
    blocker_truth = tuple(
        7 in gate_root_mask_on_z_interval(c, 1, owner0, owner1)
        for c in blockers
    )
    require(tuple(vp(c, P) for c in blockers) == (1, 2, 5),
            "blocker depth profile changed")
    require(blocker_truth == (True, False, False),
            "strict exclusive-owner statuses changed")
    print("fixed canonical-head hostile on z in (13/56,15/56):")
    print("  masks H,1,2,4,5,7 =", [list(hostile_masks[role]) for role in hostile_roles])
    print("  safe mask {1,2,7}; tau_H=4 selects source 7 -> head 11 PASS")
    print("  eligible future union is F_13\\{1,11}; q_star=7 is sole current head failure PASS")
    print("  nested blockers (52,169,13^5) have strict depths (1,2,5), owner c1 only PASS")
    print("  HOSTILE SCOPE: positive local fibre, not a global scalar cover")

    # Global image-pump audit.  If A_0 is the simultaneous guard/unit-safe
    # set and n(z) counts its thirteen direct inverse branches, then
    # {n>0}=T(A_0).  A scalar 5+3 cover forces this image into the union of
    # three blocker danger sets after division by 13, of measure at most 3/7.
    hostile_walls, hostile_hist, hostile_image, hostile_mass = (
        root_occupancy_wall_profile(hostile_roles)
    )
    expected_hist = {
        3: Fraction(389, 980),
        4: Fraction(451, 980),
        5: Fraction(1, 21),
        6: Fraction(1, 42),
        7: Fraction(1, 28),
        8: Fraction(1, 140),
        9: Fraction(2, 245),
        10: Fraction(1, 49),
    }
    require(len(hostile_walls) == 46 and hostile_hist == expected_hist,
            "fixed-head occupancy wall profile changed")
    require(hostile_image == 1 and min(hostile_hist) == 3,
            "fixed-head safe set no longer projects surjectively")
    require(hostile_mass == Fraction(226, 735),
            "fixed-head physical safe mass changed")
    require(Fraction(3, 7) < hostile_image,
            "three-blocker image bound no longer contradicts projection")
    print("global image-pump audit of the fixed-head hostile:")
    print("  46 walls / 45 open cells; n(z) histogram =",
          {n: str(mass) for n, mass in sorted(hostile_hist.items())})
    print("  min n=3, mu(T(A_0))=1, mu(A_0)=226/735 > 0 PASS")
    print("  scalar cover would force T(A_0) into three 1/7 blocker sets: impossible PASS")

    # Two earlier physical selector controls die by the same global test.
    controls = (
        ("THM-2558 sample", ((33, 2), (17, 1), (44, 1), (18, 1), (28, 1), (42, 1)),
         2, Fraction(107297, 329868)),
        ("THM-2561 primitive row", ((183, 2), (95, 1), (93, 1), (114, 1),
                                     (198, 1), (304, 1)),
         1, Fraction(5785336, 17784855)),
    )
    for name, roles, expected_min, expected_mass in controls:
        walls, histogram, image_mass, physical_mass = root_occupancy_wall_profile(roles)
        require(min(histogram) == expected_min and image_mass == 1,
                f"{name} no longer has full root image")
        require(physical_mass == expected_mass, f"{name} physical mass changed")
        require(Fraction(3, 7) < image_mass, f"{name} escaped blocker image bound")
        print(f"  {name}: {len(walls)} walls, min n={min(histogram)}, "
              f"mu(T(A_0))=1, mu(A_0)={physical_mass}; globally impossible PASS")

    # Complete small-role census on one universal exact wall refinement.
    # This is finite-exact, not an analytic extrapolation beyond the listed
    # speed universe.
    ordinary_small = tuple((k, 1) for k in range(1, 13))
    guard_small = tuple((h, 2) for h in (1, 3, 5, 7, 9, 11))
    all_cells, lengths, role_bits = universal_small_role_cells(
        ordinary_small + guard_small
    )
    require(len(lengths) == 221, "universal small-role wall refinement changed")

    no_big_best = (Fraction(2), None)
    no_big_minimizers = []
    no_big_cases = 0
    for Hsmall in (1, 3, 5, 7, 9):
        for qsmall in combinations(range(1, 12), 5):
            roles = ((Hsmall, 2),) + tuple((q, 1) for q in qsmall)
            roots = simultaneous_root_bits(roles, all_cells, role_bits)
            image = bit_measure(at_least_root_bits(roots, 1, all_cells), lengths)
            no_big_cases += 1
            if image < no_big_best[0]:
                no_big_best = (image, (Hsmall, qsmall))
                no_big_minimizers = [(Hsmall, qsmall)]
            elif image == no_big_best[0]:
                no_big_minimizers.append((Hsmall, qsmall))
    require(no_big_cases == 2310, "no-big census size changed")
    require(no_big_best == (Fraction(171, 245), (1, (3, 4, 5, 7, 11)))
            and no_big_minimizers == [(1, (3, 4, 5, 7, 11))],
            f"no-big projection minimum changed: {no_big_best}")
    require(no_big_best[0] > Fraction(3, 7), "no-big universe reached blocker capacity")
    print("complete no-all-root census: 2,310 tuples on an expanded 221-cell refinement PASS")
    print("  refinement also includes q=12,H=11 cross-check walls; cell count is not invariant")
    print("  min mu(T(A_0))=171/245 at H=1, q={3,4,5,7,11} > 3/7")
    print("  consequence: every genuine scalar cover has at least one all-root role")

    # Corroborating direct-coefficient census used by earlier hostiles: all
    # six roles pairwise distinct and contained in 1..12.
    distinct_cases = 0
    distinct_image_best = (Fraction(2), None)
    distinct_physical_best = (Fraction(2), None)
    for Hsmall in (1, 3, 5, 7, 9, 11):
        available = tuple(q for q in range(1, 13) if q != Hsmall)
        for qsmall in combinations(available, 5):
            roles = ((Hsmall, 2),) + tuple((q, 1) for q in qsmall)
            roots = simultaneous_root_bits(roles, all_cells, role_bits)
            image = bit_measure(at_least_root_bits(roots, 1, all_cells), lengths)
            physical = sum(bit_measure(bits, lengths) for bits in roots) / P
            distinct_image_best = min(distinct_image_best, (image, (Hsmall, qsmall)))
            distinct_physical_best = min(
                distinct_physical_best, (physical, (Hsmall, qsmall))
            )
            distinct_cases += 1
    require(distinct_cases == 2772, "pairwise-distinct small census size changed")
    require(distinct_image_best ==
            (Fraction(171, 245), (1, (3, 4, 5, 7, 11))),
            "pairwise-distinct image minimum changed")
    require(distinct_physical_best ==
            (Fraction(19841, 97020), (11, (1, 5, 7, 8, 9))),
            "pairwise-distinct physical minimum changed")
    print("pairwise-distinct 1..12 census: 2,772 tuples PASS")
    print("  image minimum 171/245; physical A_0 minimum 19841/97020")

    # If the sole all-root role is the forced q_star, an ordinary q_star can
    # delete at most two direct roots.  The five remaining small roles leave
    # >=3 roots on at least 43/77 of z, already above blocker capacity.  A
    # guard q_star can delete four, and the analogous >=5-root floor is only
    # 137/441<3/7: this is the sharp unresolved magnitude exception.
    ordinary_forced_best = (Fraction(2), None)
    ordinary_forced_cases = 0
    for Hsmall in (1, 3, 5, 7, 9):
        for qsmall in combinations(range(1, 12), 4):
            roles = ((Hsmall, 2),) + tuple((q, 1) for q in qsmall)
            roots = simultaneous_root_bits(roles, all_cells, role_bits)
            mass = bit_measure(at_least_root_bits(roots, 3, all_cells), lengths)
            ordinary_forced_best = min(
                ordinary_forced_best, (mass, (Hsmall, qsmall))
            )
            ordinary_forced_cases += 1
    require(ordinary_forced_cases == 1650, "ordinary forced-qstar census changed")
    require(ordinary_forced_best ==
            (Fraction(43, 77), (1, (3, 4, 5, 11))),
            "ordinary forced-qstar capacity floor changed")
    require(ordinary_forced_best[0] > Fraction(3, 7),
            "ordinary forced q_star escaped image contradiction")

    ordinary_extreme_roles = ((1, 2), (3, 1), (4, 1), (5, 1), (11, 1))
    ordinary_extreme_walls, ordinary_extreme_hist, _, _ = root_occupancy_wall_profile(
        ordinary_extreme_roles
    )
    require(len(ordinary_extreme_walls) == 48,
            "ordinary forced-qstar extremal wall count changed")
    require(ordinary_extreme_hist == {
        1: Fraction(12, 385),
        2: Fraction(158, 385),
        3: Fraction(1, 44),
        4: Fraction(3, 28),
        5: Fraction(9, 77),
        6: Fraction(61, 231),
        7: Fraction(1, 84),
        8: Fraction(1, 140),
        9: Fraction(6, 385),
        10: Fraction(1, 77),
    }, "ordinary forced-qstar extremal histogram changed")

    guard_forced_best = (Fraction(2), None)
    for qsmall in combinations(range(1, 12), 5):
        roles = tuple((q, 1) for q in qsmall)
        roots = simultaneous_root_bits(roles, all_cells, role_bits)
        mass = bit_measure(at_least_root_bits(roots, 5, all_cells), lengths)
        guard_forced_best = min(guard_forced_best, (mass, qsmall))
    require(guard_forced_best ==
            (Fraction(137, 441), (1, 5, 7, 8, 9)),
            "guard forced-qstar capacity floor changed")
    require(guard_forced_best[0] < Fraction(3, 7),
            "guard cardinality control no longer exposes its gap")

    # Cardinality alone is insufficient for a forced guard, but its root
    # danger mask has shape: for u=H mod 13 it is a length-three or
    # length-four consecutive arc in the u^{-1} root order.  A four-block
    # therefore contains every possible guard danger mask.  Exhaust exact
    # compatibility of the five-unit safe set with those blocks.
    guard_shape_worst = (Fraction(-1), None)
    guard_shape_maximizers = []
    guard_shape_cases = 0
    guard_remaining_min_roots = P
    for qsmall in combinations(range(1, 12), 5):
        roles = tuple((q, 1) for q in qsmall)
        roots = simultaneous_root_bits(roles, all_cells, role_bits)
        root_masks = []
        for cell in range(len(lengths)):
            mask = sum(1 << r for r in range(P) if roots[r] >> cell & 1)
            root_masks.append(mask)
            guard_remaining_min_roots = min(
                guard_remaining_min_roots, mask.bit_count()
            )
        for u in range(1, P):
            step = pow(u, -1, P)
            four_blocks = tuple(
                sum(1 << ((start + j * step) % P) for j in range(4))
                for start in range(P)
            )
            compatible = sum(
                length for mask, length in zip(root_masks, lengths)
                if any(mask & ~block == 0 for block in four_blocks)
            )
            key = (qsmall, u)
            if compatible > guard_shape_worst[0]:
                guard_shape_worst = (compatible, key)
                guard_shape_maximizers = [key]
            elif compatible == guard_shape_worst[0]:
                guard_shape_maximizers.append(key)
            guard_shape_cases += 1
    require(guard_shape_cases == 5544, "guard shape census size changed")
    require(guard_remaining_min_roots == 3,
            "five-small-unit safe mask changed minimum cardinality")
    require(guard_shape_worst ==
            (Fraction(904, 2695), ((1, 4, 5, 7, 11), 4)),
            f"guard shape maximum changed: {guard_shape_worst}")
    require(guard_shape_maximizers == [
        ((1, 4, 5, 7, 11), 4),
        ((1, 4, 5, 7, 11), 9),
    ], f"guard shape maximizers changed: {guard_shape_maximizers}")
    guard_survival = 1 - guard_shape_worst[0]
    require(guard_survival == Fraction(1791, 2695)
            and guard_survival > Fraction(3, 7),
            "guard shape survival floor failed")

    guard_extreme_roles = tuple((q, 1) for q in (1, 4, 5, 7, 11))
    guard_extreme_walls, guard_extreme_hist, _, _ = root_occupancy_wall_profile(
        guard_extreme_roles
    )
    require(len(guard_extreme_walls) == 58 and guard_extreme_hist == {
        3: Fraction(6, 49),
        4: Fraction(3581, 10780),
        5: Fraction(209, 980),
        6: Fraction(74, 539),
        7: Fraction(4, 77),
        8: Fraction(3, 28),
        9: Fraction(1, 140),
        10: Fraction(2, 245),
        11: Fraction(4, 539),
        12: Fraction(1, 77),
    }, "guard forced-qstar extremal histogram changed")
    guard_extreme_roots = simultaneous_root_bits(
        guard_extreme_roles, all_cells, role_bits
    )
    guard_extreme_blocks = tuple(
        sum(1 << ((start + j * pow(4, -1, P)) % P) for j in range(4))
        for start in range(P)
    )
    guard_compatible_by_size = {}
    for cell, length in enumerate(lengths):
        mask = sum(1 << r for r in range(P) if guard_extreme_roots[r] >> cell & 1)
        if any(mask & ~block == 0 for block in guard_extreme_blocks):
            size = mask.bit_count()
            guard_compatible_by_size[size] = (
                guard_compatible_by_size.get(size, Fraction(0)) + length
            )
    require(guard_compatible_by_size == {
        3: Fraction(6, 49),
        4: Fraction(82, 385),
    }, "guard compatible-size split changed")
    print("forced sole all-root role census:")
    print("  ordinary q_star: remaining >=3-root mass >=43/77 >3/7, impossible PASS")
    print("  guard cardinality alone gives 137/441 <3/7 (hostile control)")
    print("  guard four-block shape: 5,544 cases, compatible mass <=904/2695")
    print("  hence guard survival >=1791/2695 >3/7, impossible PASS")
    print("uniform conclusion: every scalar cover admits a pivot-eligible all-root k_a PASS")

    # Concrete exact controls for the all-large-N interval proof.  A is the
    # middle half of I_h; choose y in an open component of B_h.  Mesh 13^-N
    # is already smaller than |A|=1/26 for N>=2.
    mixing_checks = 0
    for L, top_k in ((1, 12), (2, 10)):
        for k in range(1, top_k + 1):
            for h in access_mask(k, L):
                pieces = intersect_intervals(danger_intervals(k, L), digit_interval(h))
                require(pieces, "accessible cell has no interval component")
                a, b = max(pieces, key=lambda ab: ab[1] - ab[0])
                y = (a + b) / 2
                A = central_digit_interval(h)
                for N in range(2, 6):
                    x = grid_preimage_in_interval(y, N, A)
                    require(x is not None, f"preimage grid missed A at {(k,L,h,N)}")
                    require((P**N * x) % 1 == y, "preimage identity failed")
                    require(A[0] < x < A[1] and a < y < b,
                            "mixing witness is not interior")
                    mixing_checks += 1
    print(f"exact preimage-grid cross-mixing controls: {mixing_checks} interior witnesses PASS")
    print("SCOPE: exact accessibility, image-pump, and finite small-role discharge;")
    print("all scalar rows admit an eligible all-root role; canonical-head cross-mixing follows;")
    print("no target current,")
    print("row exclusion, or LRC(14) conclusion. All exact checks passed.")


if __name__ == "__main__":
    main()
