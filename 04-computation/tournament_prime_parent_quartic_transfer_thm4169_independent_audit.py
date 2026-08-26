#!/usr/bin/env python3
"""Self-contained exact Python audit for THM-4169.

The frozen order-ten parent Q is extended by a vertex x with x -> i exactly
when bit i of z is one.  One endpoint DP on Q constructs Boolean-polynomial
coordinates for

    H(Q+x_z), all 55 child capacities, C, D, D+2C, and D-2C.

An independent literal endpoint DP is then run on every one of the 1,024
children and every coordinate is compared exactly.  Small labelled universes
are exhausted as positive controls.  Arithmetic is integral and every
load-bearing check uses ``require`` rather than ``assert``, so ``python -O``
executes the same audit.
"""

from __future__ import annotations

from collections import Counter
from hashlib import sha256
import gc


FROZEN_PARENT_BITS = (
    "111111101111111111111101111110110111110111111"
)


class AuditFailure(RuntimeError):
    """Raised when an exact audit obligation fails."""


def require(condition: bool, message: str) -> None:
    """Keep validation live under both ordinary and optimized Python."""
    if not condition:
        raise AuditFailure(message)


def adj_from_bits(bits: str) -> list[int]:
    """Decode upper-triangle bits; bit 1 at (i,j) means i -> j."""
    require(isinstance(bits, str), "tournament code must be a string")
    require(all(ch in "01" for ch in bits), "tournament code is not binary")
    edge_count = len(bits)
    n = 1
    while n * (n - 1) // 2 < edge_count:
        n += 1
    require(
        n * (n - 1) // 2 == edge_count,
        f"bit length {edge_count} is not triangular",
    )
    out = [0] * n
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if bits[k] == "1":
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return out


def bits_from_code(q: int, code: int) -> str:
    edge_count = q * (q - 1) // 2
    return "".join("1" if (code >> k) & 1 else "0" for k in range(edge_count))


def extend(out: list[int], pattern: int) -> list[int]:
    """Add x=n, with x -> i iff bit i of pattern is one."""
    q = len(out)
    require(0 <= pattern < (1 << q), f"attachment pattern {pattern} out of range")
    child = list(out) + [0]
    for i in range(q):
        if (pattern >> i) & 1:
            child[q] |= 1 << i
        else:
            child[i] |= 1 << q
    return child


def is_strong(out: list[int]) -> bool:
    n = len(out)
    if n == 0:
        return False
    full = (1 << n) - 1

    def reaches_all(reverse: bool) -> bool:
        seen = 1
        frontier = 1
        while frontier:
            ubit = frontier & -frontier
            frontier ^= ubit
            u = ubit.bit_length() - 1
            if reverse:
                neighbors = 0
                for v in range(n):
                    if (out[v] >> u) & 1:
                        neighbors |= 1 << v
            else:
                neighbors = out[u]
            new = neighbors & full & ~seen
            seen |= new
            frontier |= new
        return seen == full

    return reaches_all(False) and reaches_all(True)


def module_witness(out: list[int]) -> int | None:
    """Return a proper nontrivial homogeneous set, if one exists."""
    n = len(out)
    full = (1 << n) - 1
    for mask in range(1, full):
        if mask.bit_count() < 2:
            continue
        outside = full ^ mask
        vertices = outside
        homogeneous = True
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            v = bit.bit_length() - 1
            relation = out[v] & mask
            if relation and relation != mask:
                homogeneous = False
                break
        if homogeneous:
            return mask
    return None


def is_prime(out: list[int]) -> bool:
    return len(out) >= 3 and module_witness(out) is None


def endpoint_dp(out: list[int]) -> tuple[list[list[int]], list[list[int]]]:
    """Count Hamilton paths of every induced subset by either endpoint."""
    n = len(out)
    size = 1 << n
    incoming = [0] * n
    for u in range(n):
        vertices = out[u]
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            incoming[bit.bit_length() - 1] |= 1 << u

    ending = [[0] * n for _ in range(size)]
    starting = [[0] * n for _ in range(size)]
    for v in range(n):
        ending[1 << v][v] = 1
        starting[1 << v][v] = 1

    for mask in range(1, size):
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit

            predecessors = rest & incoming[v]
            total = 0
            while predecessors:
                ubit = predecessors & -predecessors
                predecessors ^= ubit
                total += ending[rest][ubit.bit_length() - 1]
            ending[mask][v] += total

            successors = rest & out[v]
            total = 0
            while successors:
                ubit = successors & -successors
                successors ^= ubit
                total += starting[rest][ubit.bit_length() - 1]
            starting[mask][v] += total
    return ending, starting


def boundary_tables(
    out: list[int], ending: list[list[int]], starting: list[list[int]]
) -> tuple[list[list[int]], list[list[int]]]:
    """Path counts that can sit immediately before/after each boundary."""
    n = len(out)
    size = 1 << n
    incoming = [0] * n
    for u in range(n):
        vertices = out[u]
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            incoming[bit.bit_length() - 1] |= 1 << u

    before = [[0] * n for _ in range(size)]
    after = [[0] * n for _ in range(size)]
    before[0] = [1] * n
    after[0] = [1] * n
    for mask in range(1, size):
        end_row = ending[mask]
        start_row = starting[mask]
        before_row = before[mask]
        after_row = after[mask]
        for boundary in range(n):
            vertices = mask & incoming[boundary]
            total = 0
            while vertices:
                bit = vertices & -vertices
                vertices ^= bit
                total += end_row[bit.bit_length() - 1]
            before_row[boundary] = total

            vertices = mask & out[boundary]
            total = 0
            while vertices:
                bit = vertices & -vertices
                vertices ^= bit
                total += start_row[bit.bit_length() - 1]
            after_row[boundary] = total
    return before, after


def capacities_from_tables(
    out: list[int], before: list[list[int]], after: list[list[int]]
) -> list[list[int]]:
    """Literal exposed-endpoint convolution for every unordered pair."""
    n = len(out)
    full = (1 << n) - 1
    capacities = [[0] * n for _ in range(n)]
    for x in range(n):
        for y in range(x + 1, n):
            remainder = full ^ (1 << x) ^ (1 << y)
            left = remainder
            value = 0
            while True:
                right = remainder ^ left
                value += before[left][x] * after[right][y]
                value += before[left][y] * after[right][x]
                if left == 0:
                    break
                left = (left - 1) & remainder
            capacities[x][y] = value
            capacities[y][x] = value
    return capacities


def literal_caps(out: list[int]) -> tuple[int, list[list[int]]]:
    """Fresh literal H and capacities for one tournament."""
    ending, starting = endpoint_dp(out)
    before, after = boundary_tables(out, ending, starting)
    hamilton_paths = sum(ending[(1 << len(out)) - 1])
    return hamilton_paths, capacities_from_tables(out, before, after)


# Boolean multilinear polynomials are sparse dictionaries mask -> coefficient.


def padd(a: dict[int, int], b: dict[int, int]) -> dict[int, int]:
    result = dict(a)
    for monomial, coefficient in b.items():
        result[monomial] = result.get(monomial, 0) + coefficient
        if result[monomial] == 0:
            del result[monomial]
    return result


def pscale(a: dict[int, int], scalar: int) -> dict[int, int]:
    if scalar == 0:
        return {}
    return {
        monomial: scalar * coefficient
        for monomial, coefficient in a.items()
        if scalar * coefficient != 0
    }


def pmul(
    a: dict[int, int],
    b: dict[int, int],
    max_degree: int | None = None,
    context: str = "polynomial product",
) -> dict[int, int]:
    """Multiply in Z[z]/(z_i^2-z_i), rejecting rather than truncating."""
    result: dict[int, int] = {}
    for left_mask, left_coefficient in a.items():
        for right_mask, right_coefficient in b.items():
            monomial = left_mask | right_mask
            if max_degree is not None:
                require(
                    monomial.bit_count() <= max_degree,
                    f"{context} exceeds Boolean degree {max_degree}: mask={monomial}",
                )
            result[monomial] = (
                result.get(monomial, 0)
                + left_coefficient * right_coefficient
            )
    return {m: c for m, c in result.items() if c != 0}


def zvar(i: int) -> dict[int, int]:
    return {1 << i: 1}


def one_minus_z(i: int) -> dict[int, int]:
    return {0: 1, 1 << i: -1}


def poly_degree(poly: dict[int, int]) -> int:
    return max((mask.bit_count() for mask in poly), default=0)


def peval(poly: dict[int, int], pattern: int) -> int:
    return sum(
        coefficient
        for monomial, coefficient in poly.items()
        if monomial & ~pattern == 0
    )


def zeta_values(poly: dict[int, int], q: int) -> list[int]:
    values = [poly.get(mask, 0) for mask in range(1 << q)]
    for i in range(q):
        bit = 1 << i
        for mask in range(1 << q):
            if mask & bit:
                values[mask] += values[mask ^ bit]
    return values


def capacity_polynomials(
    out: list[int],
) -> tuple[dict[int, int], list[list[dict[int, int]]]]:
    """Build H and every capacity of Q+x_z from one endpoint DP of Q."""
    q = len(out)
    size = 1 << q
    full = size - 1
    ending, starting = endpoint_dp(out)
    before, after = boundary_tables(out, ending, starting)

    # Hamilton paths on Q[mask]+x with an old endpoint v.
    mixed_end = [[{} for _ in range(q)] for _ in range(size)]
    mixed_start = [[{} for _ in range(q)] for _ in range(size)]

    def end_at_x(mask: int) -> dict[int, int]:
        if mask == 0:
            return {0: 1}
        answer = {0: sum(ending[mask])}
        for u in range(q):
            if (mask >> u) & 1 and ending[mask][u] != 0:
                answer = padd(answer, {1 << u: -ending[mask][u]})
        return answer

    def start_at_x(mask: int) -> dict[int, int]:
        if mask == 0:
            return {0: 1}
        return {
            1 << u: starting[mask][u]
            for u in range(q)
            if (mask >> u) & 1 and starting[mask][u] != 0
        }

    paths_end_x = [end_at_x(mask) for mask in range(size)]
    paths_start_x = [start_at_x(mask) for mask in range(size)]

    for mask in range(1, size):
        vertices = mask
        while vertices:
            bit = vertices & -vertices
            vertices ^= bit
            v = bit.bit_length() - 1
            rest = mask ^ bit
            old_end: dict[int, int] = {}
            old_start: dict[int, int] = {}
            for u in range(q):
                if not ((rest >> u) & 1):
                    continue
                if (out[u] >> v) & 1:
                    old_end = padd(old_end, mixed_end[rest][u])
                if (out[v] >> u) & 1:
                    old_start = padd(old_start, mixed_start[rest][u])
            mixed_end[mask][v] = padd(
                old_end,
                pmul(
                    zvar(v), paths_end_x[rest], 2,
                    "mixed path ending at old vertex",
                ),
            )
            mixed_start[mask][v] = padd(
                old_start,
                pmul(
                    one_minus_z(v), paths_start_x[rest], 2,
                    "mixed path starting at old vertex",
                ),
            )

    def before_mix(mask: int, boundary: int) -> dict[int, int]:
        answer: dict[int, int] = {}
        for u in range(q):
            if (mask >> u) & 1 and (out[u] >> boundary) & 1:
                answer = padd(answer, mixed_end[mask][u])
        return padd(
            answer,
            pmul(
                zvar(boundary), paths_end_x[mask], 2,
                "mixed-before boundary",
            ),
        )

    def after_mix(mask: int, boundary: int) -> dict[int, int]:
        answer: dict[int, int] = {}
        for u in range(q):
            if (mask >> u) & 1 and (out[boundary] >> u) & 1:
                answer = padd(answer, mixed_start[mask][u])
        return padd(
            answer,
            pmul(
                one_minus_z(boundary), paths_start_x[mask], 2,
                "mixed-after boundary",
            ),
        )

    mixed_before: list[list[dict[int, int] | None]] = [
        [None] * q for _ in range(size)
    ]
    mixed_after: list[list[dict[int, int] | None]] = [
        [None] * q for _ in range(size)
    ]
    for mask in range(size):
        for boundary in range(q):
            if not ((mask >> boundary) & 1):
                mixed_before[mask][boundary] = before_mix(mask, boundary)
                mixed_after[mask][boundary] = after_mix(mask, boundary)

    capacities: list[list[dict[int, int]]] = [
        [{} for _ in range(q + 1)] for _ in range(q + 1)
    ]

    # For two old endpoints, x occurs on exactly one side of the partition.
    for a in range(q):
        for b in range(a + 1, q):
            remainder = full ^ (1 << a) ^ (1 << b)
            left = remainder
            answer: dict[int, int] = {}
            while True:
                right = remainder ^ left
                bm_left_a = mixed_before[left][a]
                bm_left_b = mixed_before[left][b]
                am_right_a = mixed_after[right][a]
                am_right_b = mixed_after[right][b]
                require(
                    bm_left_a is not None and bm_left_b is not None,
                    "missing mixed-before table entry",
                )
                require(
                    am_right_a is not None and am_right_b is not None,
                    "missing mixed-after table entry",
                )
                answer = padd(answer, pscale(bm_left_a, after[right][b]))
                answer = padd(answer, pscale(bm_left_b, after[right][a]))
                answer = padd(answer, pscale(am_right_b, before[left][a]))
                answer = padd(answer, pscale(am_right_a, before[left][b]))
                if left == 0:
                    break
                left = (left - 1) & remainder
            require(
                poly_degree(answer) <= 2,
                f"old capacity ({a},{b}) exceeds degree two",
            )
            capacities[a][b] = answer
            capacities[b][a] = answer

    # For {a,x}, all remaining vertices are old.  The mutual arc does not
    # occur, so the polynomial must be affine and independent of z_a.
    x = q
    for a in range(q):
        remainder = full ^ (1 << a)
        left = remainder
        answer: dict[int, int] = {}
        while True:
            right = remainder ^ left
            answer = padd(answer, pscale(paths_end_x[left], after[right][a]))
            answer = padd(answer, pscale(paths_start_x[right], before[left][a]))
            if left == 0:
                break
            left = (left - 1) & remainder
        require(
            all(
                monomial.bit_count() <= 1 and not (monomial & (1 << a))
                for monomial in answer
            ),
            f"new capacity ({a},x) is not affine/independent of z_{a}",
        )
        capacities[a][x] = answer
        capacities[x][a] = answer

    require(
        all(
            poly_degree(capacities[i][j]) <= 2
            for i in range(q + 1)
            for j in range(i + 1, q + 1)
        ),
        "a child capacity exceeds degree two",
    )

    # THM-4114 parent cut formula:
    # H(Q+x_z)=H(Q)+sum_{i->j} c^Q_ij z_i(1-z_j).
    parent_capacities = capacities_from_tables(out, before, after)
    hpoly: dict[int, int] = {0: sum(ending[full])}
    for i in range(q):
        for j in range(q):
            if i != j and (out[i] >> j) & 1:
                contribution = pmul(
                    zvar(i), one_minus_z(j), 2, "Hamilton-path cut formula"
                )
                hpoly = padd(
                    hpoly, pscale(contribution, parent_capacities[i][j])
                )
    require(poly_degree(hpoly) <= 2, "H polynomial exceeds degree two")
    return hpoly, capacities


def packet_from_caps(
    out: list[int], hamilton_paths: int, capacities: list[list[int]]
) -> dict[str, int]:
    n = len(out)
    degree = [sum(capacities[i]) for i in range(n)]
    signed = [
        sum(
            capacities[i][j] if (out[i] >> j) & 1 else -capacities[i][j]
            for j in range(n)
            if j != i
        )
        for i in range(n)
    ]
    c_value = sum(d * s for d, s in zip(degree, signed))
    edges = [
        (i, j, capacities[i][j])
        for i in range(n)
        for j in range(i + 1, n)
    ]
    d_value = 0
    for k, (i, j, capacity) in enumerate(edges):
        for u, v, other_capacity in edges[k + 1 :]:
            if len({i, j, u, v}) == 4:
                d_value += capacity * other_capacity

    total_capacity = sum(capacity for _, _, capacity in edges)
    d_identity_numerator = (
        total_capacity * total_capacity
        - sum(value * value for value in degree)
        + sum(capacity * capacity for _, _, capacity in edges)
    )
    require(
        d_identity_numerator % 2 == 0,
        "closed D identity has odd numerator",
    )
    d_identity = d_identity_numerator // 2
    require(
        d_value == d_identity,
        f"disjoint-pair D={d_value} differs from closed identity {d_identity}",
    )
    return {
        "H": hamilton_paths,
        "C": c_value,
        "D": d_value,
        "D+2C": d_value + 2 * c_value,
        "D-2C": d_value - 2 * c_value,
    }


def packet_polynomials(
    out: list[int],
) -> tuple[
    dict[int, int],
    list[list[dict[int, int]]],
    dict[str, dict[int, int]],
]:
    hpoly, capacities = capacity_polynomials(out)
    q = len(out)
    n = q + 1
    x = q
    degree: list[dict[int, int]] = [{} for _ in range(n)]
    signed: list[dict[int, int]] = [{} for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            capacity = capacities[i][j]
            degree[i] = padd(degree[i], capacity)
            degree[j] = padd(degree[j], capacity)
            if j == x:
                orientation = padd({0: 1}, pscale(zvar(i), -2))
                signed_term = pmul(
                    orientation, capacity, 2, "signed new-edge capacity"
                )
                signed[i] = padd(signed[i], signed_term)
                signed[x] = padd(signed[x], pscale(signed_term, -1))
            elif (out[i] >> j) & 1:
                signed[i] = padd(signed[i], capacity)
                signed[j] = padd(signed[j], pscale(capacity, -1))
            else:
                signed[i] = padd(signed[i], pscale(capacity, -1))
                signed[j] = padd(signed[j], capacity)

    cpoly: dict[int, int] = {}
    for i in range(n):
        cpoly = padd(
            cpoly,
            pmul(degree[i], signed[i], 4, f"C contribution at vertex {i}"),
        )

    dpoly: dict[int, int] = {}
    edges = [
        (i, j, capacities[i][j])
        for i in range(n)
        for j in range(i + 1, n)
    ]
    for k, (i, j, capacity) in enumerate(edges):
        for u, v, other_capacity in edges[k + 1 :]:
            if len({i, j, u, v}) == 4:
                dpoly = padd(
                    dpoly,
                    pmul(
                        capacity,
                        other_capacity,
                        4,
                        f"D contribution from ({i},{j}) and ({u},{v})",
                    ),
                )

    packet = {
        "C": cpoly,
        "D": dpoly,
        "D+2C": padd(dpoly, pscale(cpoly, 2)),
        "D-2C": padd(dpoly, pscale(cpoly, -2)),
    }
    for name, poly in packet.items():
        require(poly_degree(poly) <= 4, f"{name} exceeds degree four")
    return hpoly, capacities, packet


def digest_poly(poly: dict[int, int]) -> str:
    raw = ";".join(f"{mask}:{poly[mask]}" for mask in sorted(poly))
    return sha256(raw.encode("ascii")).hexdigest()


def verify_parent_transfer(
    out: list[int], patterns: range, collect_value_digest: bool
) -> dict[str, object]:
    """Compare symbolic coordinates with fresh literal child computations."""
    q = len(out)
    hpoly, capacity_poly, packet_poly = packet_polynomials(out)
    h_values = zeta_values(hpoly, q)
    capacity_values = {
        (i, j): zeta_values(capacity_poly[i][j], q)
        for i in range(q + 1)
        for j in range(i + 1, q + 1)
    }
    packet_values = {
        name: zeta_values(poly, q) for name, poly in packet_poly.items()
    }

    value_hasher = sha256()
    pattern_count = 0
    capacity_check_count = 0
    packet_check_count = 0
    selected_packets: dict[int, dict[str, int]] = {}

    for pattern in patterns:
        child = extend(out, pattern)
        literal_h, literal_capacity = literal_caps(child)
        require(
            h_values[pattern] == literal_h,
            f"H mismatch q={q} z={pattern}: "
            f"symbolic={h_values[pattern]} literal={literal_h}",
        )

        for i in range(q + 1):
            for j in range(i + 1, q + 1):
                symbolic = capacity_values[(i, j)][pattern]
                literal = literal_capacity[i][j]
                require(
                    symbolic == literal,
                    f"capacity mismatch q={q} z={pattern} edge=({i},{j}): "
                    f"symbolic={symbolic} literal={literal}",
                )
                capacity_check_count += 1

        literal_packet = packet_from_caps(child, literal_h, literal_capacity)
        for name in ("C", "D", "D+2C", "D-2C"):
            symbolic = packet_values[name][pattern]
            literal = literal_packet[name]
            require(
                symbolic == literal,
                f"{name} mismatch q={q} z={pattern}: "
                f"symbolic={symbolic} literal={literal}",
            )
            packet_check_count += 1

        if pattern in (0, 308, (1 << q) - 1):
            selected_packets[pattern] = literal_packet

        if collect_value_digest:
            fields = [str(pattern), str(literal_h)]
            fields.extend(
                str(literal_capacity[i][j])
                for i in range(q + 1)
                for j in range(i + 1, q + 1)
            )
            fields.extend(
                str(literal_packet[name])
                for name in ("C", "D", "D+2C", "D-2C")
            )
            value_hasher.update((",".join(fields) + "\n").encode("ascii"))
        pattern_count += 1

    old_degrees = [
        poly_degree(capacity_poly[i][j])
        for i in range(q)
        for j in range(i + 1, q)
    ]
    new_degrees = [poly_degree(capacity_poly[i][q]) for i in range(q)]
    return {
        "patterns": pattern_count,
        "H_checks": pattern_count,
        "capacity_checks": capacity_check_count,
        "packet_checks": packet_check_count,
        "H_degree": poly_degree(hpoly),
        "old_capacity_max_degree": max(old_degrees, default=0),
        "new_capacity_max_degree": max(new_degrees, default=0),
        "capacity_degree_hist": tuple(
            sorted(Counter(old_degrees + new_degrees).items())
        ),
        "packet_degrees": tuple(
            (name, poly_degree(packet_poly[name]))
            for name in ("C", "D", "D+2C", "D-2C")
        ),
        "terms": tuple(
            [("H", len(hpoly))]
            + [
                (name, len(packet_poly[name]))
                for name in ("C", "D", "D+2C", "D-2C")
            ]
        ),
        "poly_hashes": tuple(
            [("H", digest_poly(hpoly))]
            + [
                (name, digest_poly(packet_poly[name]))
                for name in ("C", "D", "D+2C", "D-2C")
            ]
        ),
        "values_sha256": value_hasher.hexdigest() if collect_value_digest else None,
        "selected_packets": selected_packets,
    }


def exhaustive_transfer_controls(max_q: int = 4) -> list[tuple[int, int, int]]:
    rows: list[tuple[int, int, int]] = []
    for q in range(2, max_q + 1):
        edge_count = q * (q - 1) // 2
        parents = 1 << edge_count
        patterns = 0
        for code in range(parents):
            out = adj_from_bits(bits_from_code(q, code))
            result = verify_parent_transfer(out, range(1 << q), False)
            require(
                result["patterns"] == 1 << q,
                f"small-universe pattern count failed at q={q}",
            )
            patterns += 1 << q
        rows.append((q, parents, patterns))
    return rows


def augmentation_partition_controls(max_q: int = 5) -> tuple[int, int, tuple]:
    """Exhaust prime parents and verify the source/clone/prime partition."""
    total_prime_parents = 0
    total_patterns = 0
    rows = []
    for q in range(3, max_q + 1):
        edge_count = q * (q - 1) // 2
        prime_parents = 0
        checked_patterns = 0
        for code in range(1 << edge_count):
            out = adj_from_bits(bits_from_code(q, code))
            if not is_prime(out):
                continue
            require(is_strong(out), f"prime parent at q={q} is not strong")
            prime_parents += 1
            expected_uniform = {0, (1 << q) - 1}
            expected_clones = set()
            for v in range(q):
                base = out[v]
                expected_clones.add(base)
                expected_clones.add(base | (1 << v))
            require(
                len(expected_clones) == 2 * q,
                f"clone patterns collide for prime parent q={q} code={code}",
            )
            require(
                expected_uniform.isdisjoint(expected_clones),
                f"uniform and clone patterns overlap q={q} code={code}",
            )

            actual_uniform = set()
            actual_clones = set()
            actual_prime = set()
            for pattern in range(1 << q):
                child = extend(out, pattern)
                if not is_strong(child):
                    actual_uniform.add(pattern)
                elif is_prime(child):
                    actual_prime.add(pattern)
                else:
                    actual_clones.add(pattern)
            require(
                actual_uniform == expected_uniform,
                f"nonstrong partition mismatch q={q} code={code}",
            )
            require(
                actual_clones == expected_clones,
                f"clone partition mismatch q={q} code={code}",
            )
            require(
                len(actual_prime) == (1 << q) - 2 - 2 * q,
                f"prime attachment count mismatch q={q} code={code}",
            )
            checked_patterns += 1 << q
        rows.append((q, prime_parents, checked_patterns))
        total_prime_parents += prime_parents
        total_patterns += checked_patterns
    return total_prime_parents, total_patterns, tuple(rows)


def main() -> None:
    print("THM-4169 independent Python audit")
    print("arithmetic=exact-integer validation=explicit-require optimized_checks=live")

    partition_parents, partition_patterns, partition_rows = (
        augmentation_partition_controls(5)
    )
    require(partition_parents == 266, "unexpected prime-parent control total")
    require(partition_patterns == 8464, "unexpected augmentation pattern total")
    print(
        "augmentation_controls",
        f"rows={partition_rows}",
        f"prime_parents={partition_parents}",
        f"patterns={partition_patterns}",
        "status=PASS",
    )

    transfer_rows = exhaustive_transfer_controls(4)
    require(
        transfer_rows == [(2, 2, 8), (3, 8, 64), (4, 64, 1024)],
        f"unexpected small transfer universe {transfer_rows}",
    )
    print("small_transfer_controls", f"rows={tuple(transfer_rows)}", "status=PASS")

    parent = adj_from_bits(FROZEN_PARENT_BITS)
    require(len(parent) == 10, "frozen parent does not have order ten")
    require(is_prime(parent), "frozen order-ten parent is not prime")

    # None of the audit structures is cyclic; suppressing cyclic-GC polling
    # makes the 1,024 fresh child DPs faster without changing semantics.
    gc_was_enabled = gc.isenabled()
    if gc_was_enabled:
        gc.disable()
    try:
        result = verify_parent_transfer(parent, range(1 << 10), True)
    finally:
        if gc_was_enabled:
            gc.enable()

    require(result["patterns"] == 1024, "q10 audit did not cover all patterns")
    require(result["H_checks"] == 1024, "q10 H check count mismatch")
    require(result["capacity_checks"] == 1024 * 55, "q10 capacity check count mismatch")
    require(result["packet_checks"] == 1024 * 4, "q10 packet check count mismatch")
    require(result["H_degree"] == 2, "frozen H polynomial degree changed")
    require(
        result["old_capacity_max_degree"] == 2,
        "frozen old-capacity maximum degree changed",
    )
    require(
        result["new_capacity_max_degree"] == 1,
        "frozen new-capacity maximum degree changed",
    )
    require(
        result["packet_degrees"]
        == (("C", 4), ("D", 4), ("D+2C", 4), ("D-2C", 4)),
        "frozen packet degrees changed",
    )
    require(
        result["terms"]
        == (("H", 56), ("C", 386), ("D", 386), ("D+2C", 386), ("D-2C", 386)),
        "frozen polynomial term counts changed",
    )

    selected = result["selected_packets"]
    source_packet = selected[1023]
    require(
        source_packet["D"] == 2_735_733_720
        and min(source_packet["D+2C"], source_packet["D-2C"])
        == -1_082_674_184,
        "frozen source hostile changed",
    )
    sink_packet = selected[0]
    require(
        min(sink_packet["D+2C"], sink_packet["D-2C"])
        == 209_423_320,
        "frozen sink control changed",
    )

    print("parent", FROZEN_PARENT_BITS, "order=10", "prime=yes")
    print(
        "full_transfer",
        f"patterns={result['patterns']}",
        f"H_checks={result['H_checks']}",
        f"capacity_checks={result['capacity_checks']}",
        f"packet_checks={result['packet_checks']}",
        "status=PASS",
    )
    print(
        "degrees",
        f"H={result['H_degree']}",
        f"old_capacity_max={result['old_capacity_max_degree']}",
        f"new_capacity_max={result['new_capacity_max_degree']}",
        f"capacity_hist={result['capacity_degree_hist']}",
        f"packet={result['packet_degrees']}",
    )
    print("terms", result["terms"])
    print("polynomial_sha256", result["poly_hashes"])
    print("all_values_sha256", result["values_sha256"])
    print(
        "boundary_controls",
        f"z0_margin={min(sink_packet['D+2C'], sink_packet['D-2C'])}",
        f"z1023_twice_abs_C={2 * abs(source_packet['C'])}",
        f"z1023_D={source_packet['D']}",
        f"z1023_margin={min(source_packet['D+2C'], source_packet['D-2C'])}",
        "status=PASS",
    )
    print("overall status=PASS")


if __name__ == "__main__":
    main()
