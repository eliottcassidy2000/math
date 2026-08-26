#!/usr/bin/env python3
"""Exact replay for THM-4172.

This file has three deliberately separate lanes:

* literal tagged odd-cycle packings verify the capacity deletion transform;
* direct polynomial evaluation verifies support-degree diagonalization and
  the all-depth normalized holonomy law; and
* subset-DP capacities evaluate the named order-11/order-12 controls.

It imports no project computation module.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def decode(code, n):
    adj = [0] * n
    for bit, (i, j) in enumerate(combinations(range(n), 2)):
        if code >> bit & 1:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
    return tuple(adj)


def restrict_pair(adj, z, removed):
    keep = tuple(v for v in range(len(adj)) if v not in removed)
    where = {v: i for i, v in enumerate(keep)}
    card = [0] * len(keep)
    tensor = [[0] * len(keep) for _ in keep]
    for i, j in combinations(keep, 2):
        a, b = where[i], where[j]
        if adj[i] >> j & 1:
            card[a] |= 1 << b
        else:
            card[b] |= 1 << a
        tensor[a][b] = tensor[b][a] = z[i][j]
    return tuple(card), tuple(tuple(row) for row in tensor), keep, where


def capacity(adj):
    """Return H and c_ij=Q(i,j)+Q(j,i) using start/end subset DPs."""
    n = len(adj)
    size = 1 << n
    full = size - 1
    starts = [[0] * n for _ in range(size)]
    ends = [[0] * n for _ in range(size)]
    for u in range(n):
        starts[1 << u][u] = 1
        ends[1 << u][u] = 1
    for mask in range(1, size):
        if mask & (mask - 1) == 0:
            continue
        for u in range(n):
            if not mask >> u & 1:
                continue
            rest = mask ^ (1 << u)
            for v in range(n):
                if not rest >> v & 1:
                    continue
                if adj[u] >> v & 1:
                    starts[mask][u] += starts[rest][v]
                if adj[v] >> u & 1:
                    ends[mask][u] += ends[rest][v]
    q = [[0] * n for _ in range(n)]
    for left in range(1, full):
        right = full ^ left
        for i in range(n):
            if not left >> i & 1 or not ends[left][i]:
                continue
            for j in range(n):
                if right >> j & 1:
                    q[i][j] += ends[left][i] * starts[right][j]
    c = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        c[i][j] = c[j][i] = q[i][j] + q[j][i]
    return sum(starts[full]), tuple(tuple(row) for row in c)


def packet(adj, z):
    degrees = [sum(row) for row in z]
    fields = [
        sum(
            z[i][j] if adj[i] >> j & 1 else -z[i][j]
            for j in range(len(adj))
            if i != j
        )
        for i in range(len(adj))
    ]
    c_value = sum(h * d for h, d in zip(fields, degrees))
    edges = tuple(combinations(range(len(adj)), 2))
    d_value = sum(
        z[i][j] * z[k][l]
        for index, (i, j) in enumerate(edges)
        for k, l in edges[index + 1 :]
        if len({i, j, k, l}) == 4
    )
    return c_value, d_value


def tau(adj, z):
    c_value, d_value = packet(adj, z)
    need(d_value > 0, "positive normalized denominator")
    divisor = 2 if len(adj) % 2 == 0 else 4
    return Fraction((len(adj) - 3) * c_value, divisor * d_value)


def add_ear(adj, outgoing):
    n = len(adj)
    child = list(adj) + [0]
    for u in range(n):
        if outgoing >> u & 1:
            child[n] |= 1 << u
        else:
            child[u] |= 1 << n
    return tuple(child)


def directed_odd_cycles(adj):
    answer = []
    for length in range(3, len(adj) + 1, 2):
        for support in combinations(range(len(adj)), length):
            first = support[0]
            for tail in permutations(support[1:]):
                word = (first,) + tail
                if all(
                    adj[word[k]] >> word[(k + 1) % length] & 1
                    for k in range(length)
                ):
                    answer.append((sum(1 << u for u in support), word))
    return tuple(answer)


def odd_cycle_packings(adj):
    cycles = directed_odd_cycles(adj)
    answer = []

    def visit(start, used, chosen):
        answer.append(tuple(chosen))
        for index in range(start, len(cycles)):
            mask, word = cycles[index]
            if not used & mask:
                chosen.append((mask, word))
                visit(index + 1, used | mask, chosen)
                chosen.pop()

    visit(0, 0, [])
    return tuple(answer)


def tagged_support_histogram(adj, edge):
    """Literal W_s for the two ambient-state tags of THM-4167."""
    n = len(adj)
    x = n
    histogram = [0] * (n - 1)
    endpoint_mask = (1 << edge[0]) | (1 << edge[1])
    for outgoing_endpoint in edge:
        child = add_ear(adj, 1 << outgoing_endpoint)
        for packing in odd_cycle_packings(child):
            x_cycle = None
            used = 0
            for mask, word in packing:
                used |= mask
                if mask >> x & 1:
                    x_cycle = word
            if x_cycle is None:
                continue
            index = x_cycle.index(x)
            neighbours = frozenset(
                (x_cycle[index - 1], x_cycle[(index + 1) % len(x_cycle)])
            )
            if neighbours != frozenset(edge):
                continue
            outside = (used & ((1 << n) - 1) & ~endpoint_mask).bit_count()
            histogram[outside] += 1 << len(packing)
    return tuple(histogram)


def deletion_layers(adj, edge):
    outside = tuple(v for v in range(len(adj)) if v not in edge)
    zero_tensor = tuple(tuple(0 for _ in adj) for _ in adj)
    layers = []
    for r in range(len(outside) + 1):
        total = 0
        for removed in combinations(outside, r):
            card, _, keep, where = restrict_pair(adj, zero_tensor, frozenset(removed))
            _, c_card = capacity(card)
            total += c_card[where[edge[0]]][where[edge[1]]]
        layers.append(total)
    return tuple(layers)


def invert_layers(layers):
    """Coefficient extraction from A(x-1)=sum_s W_s x^(N-s)."""
    n_outside = len(layers) - 1
    answer = []
    for s in range(n_outside + 1):
        degree = n_outside - s
        answer.append(
            sum(
                (-1) ** (r - degree) * comb(r, degree) * layers[r]
                for r in range(degree, n_outside + 1)
            )
        )
    return tuple(answer)


def audit_capacity_tomography():
    rows = layer_gates = 0
    for n in range(2, 5):
        for code in range(1 << (n * (n - 1) // 2)):
            adj = decode(code, n)
            _, c_parent = capacity(adj)
            for edge in combinations(range(n), 2):
                histogram = tagged_support_histogram(adj, edge)
                layers = deletion_layers(adj, edge)
                need(sum(histogram) == c_parent[edge[0]][edge[1]], "tagged atom sum")
                predicted = tuple(
                    sum(comb(n - 2 - s, r) * weight for s, weight in enumerate(histogram))
                    for r in range(n - 1)
                )
                need(layers == predicted, "capacity binomial deletion transform")
                need(invert_layers(layers) == histogram, "capacity support tomography")
                need(all(weight >= 0 for weight in invert_layers(layers)), "positive support weights")
                rows += 1
                layer_gates += len(layers)
    need((rows, layer_gates) == (410, 1202), "small-universe gate count")
    return rows, layer_gates


def binomial_or_zero(a, b):
    return comb(a, b) if 0 <= b <= a else 0


def audit_multideletion_algebra():
    rows = restrictions = 0
    for n in range(5, 10):
        adj = decode(((1 << (n * (n - 1) // 2)) - 1) // 3, n)
        z = [[0] * n for _ in range(n)]
        for index, (i, j) in enumerate(combinations(range(n), 2)):
            z[i][j] = z[j][i] = 1 + (29 * index + 11 * n) % 31
        z = tuple(tuple(row) for row in z)
        parent_c, parent_d = packet(adj, z)
        parent_tau = tau(adj, z)
        for r in range(1, n - 3):
            deck = []
            for removed in combinations(range(n), r):
                child, child_z, _, _ = restrict_pair(adj, z, frozenset(removed))
                child_c, child_d = packet(child, child_z)
                deck.append((child_d, tau(child, child_z), child_c))
            sum_c = sum(row[2] for row in deck)
            sum_d = sum(row[0] for row in deck)
            need(sum_c == binomial_or_zero(n - 3, r) * parent_c, "C support degree three")
            need(sum_d == binomial_or_zero(n - 4, r) * parent_d, "D support degree four")
            weighted = sum(Fraction(d, sum_d) * child_tau for d, child_tau, _ in deck)
            child_order = n - r
            kappa_parent = 2 if n % 2 == 0 else 4
            kappa_child = 2 if child_order % 2 == 0 else 4
            need(parent_tau == Fraction(kappa_child, kappa_parent) * weighted,
                 "all-depth parity holonomy")
            if r % 2 == 0:
                need(parent_tau == weighted, "same-parity barycenter")
            rows += 1
            restrictions += len(deck)
    need(rows == 15, "multideletion row count")
    return rows, restrictions


def multideck(adj, z, r):
    rows = []
    for removed in combinations(range(len(adj)), r):
        child, child_z, _, _ = restrict_pair(adj, z, frozenset(removed))
        c_value, d_value = packet(child, child_z)
        rows.append((removed, tau(child, child_z), c_value, d_value))
    total_d = sum(row[3] for row in rows)
    mean = sum(Fraction(row[3], total_d) * row[1] for row in rows)
    return mean, rows


def named_controls():
    prime11 = decode(3169369058263173, 11)
    h11, z11 = capacity(prime11)
    need(h11 == 23685, "prime order-eleven H")
    parent11 = tau(prime11, z11)
    mean11, pairs11 = multideck(prime11, z11, 2)
    need(parent11 == mean11 == Fraction(1055017002, 11090656697),
         "prime order-eleven two-deletion mean")
    need(sum(abs(row[1]) >= 1 for row in pairs11) == 0,
         "prime order-eleven pair centrality")
    worst11 = max(pairs11, key=lambda row: abs(row[1]))
    need((worst11[0], worst11[1]) == ((7, 10), Fraction(1629373665, 9484374388)),
         "prime order-eleven worst pair")

    hostile12 = (3070, 3644, 3704, 3824, 4064, 4032,
                 3970, 3846, 3598, 1024, 2049, 512)
    h12, z12 = capacity(hostile12)
    need(h12 == 27759, "order-twelve hostile H")
    parent12 = tau(hostile12, z12)
    mean_one, one_cards = multideck(hostile12, z12, 1)
    mean_two, pairs12 = multideck(hostile12, z12, 2)
    need(parent12 == 2 * mean_one == mean_two
         == Fraction(-53092739331, 40435524866), "hostile multideletion identities")
    need(sum(abs(row[1]) >= 1 for row in one_cards) == 0,
         "hostile one-deletion false negative")
    need(sum(abs(row[1]) >= 1 for row in pairs12) == 65,
         "hostile two-deletion exposure count")
    central = tuple(row for row in pairs12 if abs(row[1]) < 1)
    need(len(central) == 1 and central[0][:2]
         == ((9, 11), Fraction(-4461429, 4143186457)), "unique central hostile pair")
    worst12 = max(pairs12, key=lambda row: abs(row[1]))
    need(worst12[:2] == ((9, 10), Fraction(-770665315, 458646486)),
         "hostile worst pair")
    return {
        "prime11": (parent11, len(pairs11), worst11[:2]),
        "hostile12": (parent12, len(one_cards), 65, central[0][:2], worst12[:2]),
    }


def main():
    print("capacity_tomography", audit_capacity_tomography())
    print("multideletion_algebra", audit_multideletion_algebra())
    print("named_controls", named_controls())
    print("THM4172_AUDIT_PASS")


if __name__ == "__main__":
    main()
