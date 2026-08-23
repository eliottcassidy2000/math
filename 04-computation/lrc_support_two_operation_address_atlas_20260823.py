#!/usr/bin/env python3
"""Exact operation-address atlas for THM-3743/THM-3793/THM-3823.

The all-scale decoder theorem is elementary; the finite computations audit
the complete support-two atlas a+b<=356, its
inert/cube-free subatlas, cube-sum collision fibres, and the induced
Stern--Brocot forest.  Gates remain active under ``python -O``.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from math import gcd, isqrt
import sys


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CAP = 356


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor(n: int) -> dict[int, int]:
    require(n >= 1, "factor domain")
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def phi(n: int) -> int:
    answer = n
    for p in factor(n):
        answer -= answer // p
    return answer


def admissible_shell(s: int) -> bool:
    fs = factor(s)
    return s >= 3 and bool(fs) and all(p % 3 == 2 and e <= 2 for p, e in fs.items())


def parent(pair: tuple[int, int]) -> tuple[int, int] | None:
    a, b = pair
    require(1 <= a < b and gcd(a, b) == 1, "parent input")
    if (a, b) == (1, 2):
        return None
    if b > 2 * a:
        return a, b - a
    if b < 2 * a:
        return b - a, a
    raise RuntimeError("only the root can lie on b=2a")


def children(pair: tuple[int, int]) -> tuple[tuple[int, int], tuple[int, int]]:
    a, b = pair
    return (a, a + b), (b, a + b)


def path_and_heap(pair: tuple[int, int]) -> tuple[str, int]:
    reverse_letters: list[str] = []
    cursor = pair
    while cursor != (1, 2):
        a, b = cursor
        if b > 2 * a:
            reverse_letters.append("L")
            cursor = (a, b - a)
        elif b < 2 * a:
            reverse_letters.append("R")
            cursor = (b - a, a)
        else:
            raise RuntimeError("bad Euclidean cone")
    word = "".join(reversed(reverse_letters))
    heap = 1
    for letter in word:
        heap = 2 * heap + (letter == "R")
    return word, heap


def decode_from_factorization(fm: dict[int, int]) -> tuple[int, int, int, int] | None:
    """Decode (g,s,X,Y) on the inert-scale/cube-free-shell carrier."""
    if 3 in fm:
        return None
    g = 1
    s = 1
    q = 1
    for p, exponent in fm.items():
        if p % 3 == 2:
            quotient, remainder = divmod(exponent, 3)
            g *= p ** quotient
            s *= p ** remainder
        elif p % 3 == 1:
            q *= p ** exponent
        else:
            return None
    numerator = 4 * q - s * s
    if numerator <= 0 or numerator % 3:
        return None
    delta2 = numerator // 3
    delta = isqrt(delta2)
    if delta * delta != delta2 or not (0 < delta < s) or (s - delta) % 2:
        return None
    x = (s - delta) // 2
    y = (s + delta) // 2
    if gcd(x, y) != 1:
        return None
    return g, s, x, y


def merged_factors(left: dict[int, int], right: dict[int, int], multiplier: int) -> dict[int, int]:
    out = dict(left)
    for p, exponent in right.items():
        out[p] = out.get(p, 0) + multiplier * exponent
    return out


def main() -> None:
    checks = 0

    full_nodes = [
        (a, s - a)
        for s in range(3, CAP + 1)
        for a in range(1, (s + 1) // 2)
        if a < s - a and gcd(a, s) == 1
    ]
    full_set = set(full_nodes)
    require(len(full_nodes) == len(full_set) == 19_314, "full atlas count")

    shell_sizes = Counter(a + b for a, b in full_nodes)
    for s in range(3, CAP + 1):
        require(shell_sizes[s] == phi(s) // 2, f"totient shell {s}")
        checks += 1

    # Dense shell-first rank, and its exact inverse.
    dense_rank: dict[tuple[int, int], int] = {}
    dense_inverse: dict[int, tuple[int, int]] = {}
    for rank, pair in enumerate(full_nodes, 1):
        dense_rank[pair] = rank
        dense_inverse[rank] = pair
    require(set(dense_inverse) == set(range(1, 19_315)), "dense full rank interval")
    checks += 2 * len(full_nodes)

    # The full atlas is a prefix-closed finite Stern--Brocot subtree.
    depth_histogram: Counter[int] = Counter()
    heaps: set[int] = set()
    max_depth = -1
    deepest: list[tuple[int, int]] = []
    child_inside = Counter()
    for pair in full_nodes:
        p = parent(pair)
        if p is not None:
            require(p in full_set, f"prefix closure at {pair}")
        word, heap = path_and_heap(pair)
        require(heap not in heaps, f"heap collision at {pair}")
        heaps.add(heap)
        require(len(word) <= sum(pair) - 3, f"depth bound at {pair}")
        if p is not None:
            pword, pheap = path_and_heap(p)
            require(word[:-1] == pword and heap // 2 == pheap, "heap parent law")
        left, right = children(pair)
        for letter, child in (("L", left), ("R", right)):
            s, delta = pair[0] + pair[1], pair[1] - pair[0]
            child_s, child_delta = child[0] + child[1], child[1] - child[0]
            if letter == "L":
                require(2 * child_s == 3 * s - delta, "L affine shell law")
                require(2 * child_delta == s + delta, "L affine gap law")
                require(path_and_heap(child)[1] == 2 * heap, "L heap append")
            else:
                require(2 * child_s == 3 * s + delta, "R affine shell law")
                require(2 * child_delta == s - delta, "R affine gap law")
                require(path_and_heap(child)[1] == 2 * heap + 1, "R heap append")
            if child in full_set:
                child_inside[letter] += 1
        depth = len(word)
        depth_histogram[depth] += 1
        if depth > max_depth:
            max_depth = depth
            deepest = [pair]
        elif depth == max_depth:
            deepest.append(pair)
        checks += 12
    require(max_depth == 353 and deepest == [(1, 355)], "sharp depth boundary")

    # Reserve seven low bits for one of the 78 unordered labelled placements.
    # This keeps the full ratio-tree operations affine on a single natural
    # number; the actual common speed scale remains an explicit sidecar.
    tree_label_checks = 0
    for pair in full_nodes:
        _, heap = path_and_heap(pair)
        for label_index in range(78):
            address = 128 * heap + label_index
            recovered_label = address % 128
            recovered_heap = (address - recovered_label) // 128
            require((recovered_heap, recovered_label) == (heap, label_index), "tree-label decoder")
            require(2 * address - label_index == 128 * (2 * heap) + label_index,
                    "tree L address law")
            require(2 * address - label_index + 128 == 128 * (2 * heap + 1) + label_index,
                    "tree R address law")
            tree_label_checks += 3

    # Exact cube-sum fibres on the complete primitive ratio universe.
    cube_fibres: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for pair in full_nodes:
        cube_fibres[pair[0] ** 3 + pair[1] ** 3].append(pair)
    collision_values = {m: pairs for m, pairs in cube_fibres.items() if len(pairs) > 1}
    collision_histogram = Counter(len(pairs) for pairs in cube_fibres.values())
    require(cube_fibres[1729] == [(1, 12), (9, 10)], "Ramanujan hostile fibre")
    require(max(collision_histogram) >= 2, "positive collision control")
    collision_inert_profiles = Counter()
    for m, pairs in collision_values.items():
        inert_part = 1
        for p, exponent in factor(m).items():
            if p % 3 == 2:
                inert_part *= p ** exponent
        require(all(sum(pair) % inert_part == 0 for pair in pairs), "shared inert divisor")
        collision_inert_profiles["inert_free" if inert_part == 1 else "nontrivial_inert"] += 1

    admissible_sums = [s for s in range(3, CAP + 1) if admissible_shell(s)]
    selected_nodes = [pair for pair in full_nodes if admissible_shell(sum(pair))]
    selected_set = set(selected_nodes)
    require(len(admissible_sums) == 94, "admissible shell count")
    require(len(selected_nodes) == 5_855, "selected ratio count")
    require(sum(phi(s) // 2 for s in admissible_sums) == 5_855, "selected totient mass")

    selected_dense_rank = {pair: rank for rank, pair in enumerate(selected_nodes, 1)}
    require(set(selected_dense_rank.values()) == set(range(1, 5_856)), "dense selected interval")

    # Prime-colour decoder on every selected primitive value.
    decoded_values: set[int] = set()
    split_cofactor_factor_count = 0
    for pair in selected_nodes:
        a, b = pair
        s = a + b
        q = a * a - a * b + b * b
        m = a ** 3 + b ** 3
        fq = factor(q)
        require(all(p % 3 == 1 for p in fq), f"non-split cofactor {pair}")
        split_cofactor_factor_count += len(fq)
        decoded = decode_from_factorization(factor(m))
        require(decoded == (1, s, a, b), f"primitive decoder {pair}: {decoded}")
        require(m not in decoded_values, "selected cube collision")
        decoded_values.add(m)
        checks += 5 + len(fq)
    require(len(decoded_values) == 5_855, "selected singleton image")

    # The unused prime 3 carries the 78 unordered runner-pair labels.  This is
    # sparse rather than dense, but common inert dilation remains multiplication
    # by a cube on the single natural-number address.
    labelled_addresses: set[int] = set()
    label_decode_checks = 0
    for m in decoded_values:
        for label_index in range(78):
            address = 3 ** label_index * m
            require(address not in labelled_addresses, "labelled natural-address collision")
            labelled_addresses.add(address)
            stripped = address
            recovered_label = 0
            while stripped % 3 == 0:
                recovered_label += 1
                stripped //= 3
            require(recovered_label == label_index and stripped == m, "3-adic label decoder")
            label_decode_checks += 2
    require(len(labelled_addresses) == 456_690, "labelled address census")

    # All finite inert dilations through 100: exponents split as 3*scale + shell digit.
    inert_scales = [
        g for g in range(1, 101)
        if g == 1 or all(p % 3 == 2 for p in factor(g))
    ]
    scaled_decoder_checks = 0
    for a, b in selected_nodes:
        fm = factor(a ** 3 + b ** 3)
        for g in inert_scales:
            scaled_factors = merged_factors(fm, factor(g), 3)
            decoded = decode_from_factorization(scaled_factors)
            require(decoded == (g, a + b, a, b), f"scaled decoder {(g, a, b)}")
            # Inert dilation is exactly a cube action and leaves the ratio state fixed.
            require((g * a) ** 3 + (g * b) ** 3 == g ** 3 * (a ** 3 + b ** 3), "cube action")
            scaled_decoder_checks += 2

    # Selected nodes form an induced forest, not a subtree or operation-closed carrier.
    selected_children: dict[tuple[int, int], list[tuple[int, int]]] = defaultdict(list)
    selected_parent_edges = 0
    branch_landings_unbounded = Counter()
    branch_landings_inside = Counter()
    unbounded_child_profile = Counter()
    for pair in selected_nodes:
        p = parent(pair)
        if p in selected_set:
            selected_children[p].append(pair)
            selected_parent_edges += 1
        left, right = children(pair)
        unbounded_count = 0
        for letter, child in (("L", left), ("R", right)):
            if admissible_shell(sum(child)):
                unbounded_count += 1
                branch_landings_unbounded[letter] += 1
                if sum(child) <= CAP:
                    require(child in selected_set, "selected child membership")
                    branch_landings_inside[letter] += 1
        unbounded_child_profile[unbounded_count] += 1

    roots = [pair for pair in selected_nodes if parent(pair) not in selected_set]
    require(len(roots) + selected_parent_edges == len(selected_nodes), "forest Euler identity")
    component_sizes: list[tuple[int, tuple[int, int], int]] = []
    for root in roots:
        queue = deque([(root, 0)])
        size = 0
        max_internal_depth = 0
        while queue:
            node, internal_depth = queue.popleft()
            size += 1
            max_internal_depth = max(max_internal_depth, internal_depth)
            queue.extend((child, internal_depth + 1) for child in selected_children[node])
        component_sizes.append((size, root, max_internal_depth))
    require(sum(size for size, _, _ in component_sizes) == len(selected_nodes), "component partition")
    component_sizes.sort(reverse=True)

    # Small canonical hostiles for the respective quotients.
    require(children((1, 4)) == ((1, 5), (4, 5)), "first branch hostile pairs")
    require(not admissible_shell(6) and not admissible_shell(9), "first branch hostile shells")
    require(decode_from_factorization(factor(1729)) is None, "1729 must stay outside selected decoder")
    split_scaled = 7 ** 3 * (1 ** 3 + 4 ** 3)
    require(split_scaled == 22_295, "split-scale hostile value")
    require(decode_from_factorization(factor(split_scaled)) is None, "split scale must not fake inert scale")
    require(decode_from_factorization(factor(515_375)) is None, "exponent-three hostile must fail")
    ap_ratio_hits = 0
    ap_all_scale_hits = 0
    for u in range(1, 14):
        for v in range(u + 1, 14):
            g = gcd(u, v)
            a, b = u // g, v // g
            if not admissible_shell(a + b):
                continue
            ap_ratio_hits += 1
            if g == 1 or all(p % 3 == 2 for p in factor(g)):
                require(decode_from_factorization(factor(u ** 3 + v ** 3)) == (g, a + b, a, b),
                        "AP all-scale decoder")
                ap_all_scale_hits += 1
    require(ap_ratio_hits == 30 and ap_all_scale_hits == 27, "AP hostile census")
    checks += 8 + ap_all_scale_hits

    leaf_profile = Counter(len(selected_children[pair]) for pair in selected_nodes)
    top_components = component_sizes[:10]
    first_collisions = sorted(collision_values.items())[:10]

    print("GRAVER RATIO / PRIME-COLOUR ADDRESS PROBE")
    print(f"universe=coprime 1<=a<b, a+b<={CAP}")
    print(f"full_nodes={len(full_nodes)}; shell_formula=(1/2)sum_{{s=3}}^{CAP}phi(s)")
    print(f"dense_full_address=1..{len(full_nodes)}")
    print(f"stern_brocot_root=(1,2); prefix_closed=yes; max_depth={max_depth}; deepest={deepest}")
    print(f"tree_label_address=128*heap+k; count={78*len(full_nodes)}; checks={tree_label_checks}")
    print("tree_label_operations: L(A)=2A-k, R(A)=2A-k+128; scale remains a sidecar")
    print(f"full_child_edges_inside=L:{child_inside['L']},R:{child_inside['R']}")
    print(f"cube_distinct_values={len(cube_fibres)}; fibre_histogram={dict(sorted(collision_histogram.items()))}")
    print(f"cube_collision_values={len(collision_values)}; collision_excess={len(full_nodes)-len(cube_fibres)}")
    print(f"collision_inert_profiles={dict(sorted(collision_inert_profiles.items()))}")
    print(f"first_cube_collisions={first_collisions}")
    print(f"admissible_shells={len(admissible_sums)}; selected_nodes={len(selected_nodes)}")
    print(f"dense_selected_address=1..{len(selected_nodes)}")
    print(f"selected_cube_values={len(decoded_values)}; split_cofactor_prime_factor_checks={split_cofactor_factor_count}")
    print(f"labelled_prime3_addresses={len(labelled_addresses)}; label_decode_checks={label_decode_checks}")
    print("labelled_address=3^k*g^3*(a^3+b^3), 0<=k<78; v_3 recovers the labelled pair")
    print(f"inert_scales_through_100={inert_scales}; scaled_decoder_checks={scaled_decoder_checks}")
    print("decoder: v_p(m)=3*v_p(g)+v_p(s) on p=2 mod 3; quotient=scale, remainder=shell")
    print("decoder: q=split-prime part; delta^2=(4q-s^2)/3; (a,b)=((s-delta)/2,(s+delta)/2)")
    print(f"selected_parent_edges={selected_parent_edges}; selected_forest_roots={len(roots)}")
    print(f"selected_child_profile={dict(sorted(leaf_profile.items()))}")
    print(f"selected_unbounded_child_profile={dict(sorted(unbounded_child_profile.items()))}")
    print(f"selected_branch_landings_inside=L:{branch_landings_inside['L']},R:{branch_landings_inside['R']}")
    print(f"selected_branch_landings_unbounded=L:{branch_landings_unbounded['L']},R:{branch_landings_unbounded['R']}")
    print(f"largest_selected_components={top_components}")
    print("hostile_branch=(1,4) shell5 -> shells6,9; neither admissible")
    print("hostile_scalar=1729=(1,12)=(9,10); cube value is not a full-atlas address")
    print("hostile_scale=7*(1,4): 22295; split scale is entangled with the Eisenstein cofactor")
    print("hostile_exponent=515375: primitive shell 125=5^3; valuation remainder erases the shell")
    print("hostile_AP=(1,...,13) is safe at t=1/14 but has 30 selected ratio placements, 27 with inert scale")
    print(f"active_checks={checks + scaled_decoder_checks + label_decode_checks + tree_label_checks}")
    print("RESULT PASS")


if __name__ == "__main__":
    main()
