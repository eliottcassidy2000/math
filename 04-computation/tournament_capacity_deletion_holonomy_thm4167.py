#!/usr/bin/env python3
"""Exact audit for THM-4167 capacity deletion and parity holonomy.

The implementation uses subset dynamic programming for capacities, literal
permutations for marked exposed words, and literal directed-cycle packings for
the small OCF audit.  It imports no project computation module.
"""

from fractions import Fraction
from itertools import combinations, permutations


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def decode(code, n):
    adj = [0] * n
    bit = 0
    for i, j in combinations(range(n), 2):
        if code >> bit & 1:
            adj[i] |= 1 << j
        else:
            adj[j] |= 1 << i
        bit += 1
    return tuple(adj)


def encode(adj):
    code = 0
    for bit, (i, j) in enumerate(combinations(range(len(adj)), 2)):
        code |= ((adj[i] >> j) & 1) << bit
    return code


def delete_vertex(adj, v):
    keep = tuple(u for u in range(len(adj)) if u != v)
    where = {u: i for i, u in enumerate(keep)}
    card = [0] * len(keep)
    for u, w in combinations(keep, 2):
        if adj[u] >> w & 1:
            card[where[u]] |= 1 << where[w]
        else:
            card[where[w]] |= 1 << where[u]
    return tuple(card), keep, where


def is_strong(adj):
    full = (1 << len(adj)) - 1
    for root in range(len(adj)):
        seen = todo = 1 << root
        while todo:
            bit = todo & -todo
            todo ^= bit
            u = bit.bit_length() - 1
            fresh = adj[u] & ~seen
            seen |= fresh
            todo |= fresh
        if seen != full:
            return False
    return True


def is_prime(adj):
    full = (1 << len(adj)) - 1
    for module in range(1, full):
        if not 2 <= module.bit_count() < len(adj):
            continue
        outside = full ^ module
        ok = True
        while outside and ok:
            bit = outside & -outside
            outside ^= bit
            row = adj[bit.bit_length() - 1] & module
            ok = row == 0 or row == module
        if ok:
            return False
    return True


def capacity(adj):
    """Return H and c_ij=Q_ij+Q_ji by independent start/end subset DPs."""
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
            ub = 1 << u
            if not mask & ub:
                continue
            rest = mask ^ ub
            todo = rest
            while todo:
                vb = todo & -todo
                todo ^= vb
                v = vb.bit_length() - 1
                if adj[u] & vb:
                    starts[mask][u] += starts[rest][v]
                if adj[v] & ub:
                    ends[mask][u] += ends[rest][v]
    q = [[0] * n for _ in range(n)]
    for left in range(1, full):
        right = full ^ left
        todo_left = left
        while todo_left:
            ib = todo_left & -todo_left
            todo_left ^= ib
            i = ib.bit_length() - 1
            if not ends[left][i]:
                continue
            todo_right = right
            while todo_right:
                jb = todo_right & -todo_right
                todo_right ^= jb
                j = jb.bit_length() - 1
                q[i][j] += ends[left][i] * starts[right][j]
    c = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        c[i][j] = c[j][i] = q[i][j] + q[j][i]
    return sum(starts[full]), tuple(tuple(row) for row in c)


def packet(adj, z):
    degree = [sum(row) for row in z]
    field = [
        sum(z[i][j] if adj[i] >> j & 1 else -z[i][j]
            for j in range(len(adj)) if i != j)
        for i in range(len(adj))
    ]
    cval = sum(field[i] * degree[i] for i in range(len(adj)))
    edges = tuple(combinations(range(len(adj)), 2))
    dval = sum(
        z[i][j] * z[k][l]
        for a, (i, j) in enumerate(edges)
        for k, l in edges[a + 1:]
        if len({i, j, k, l}) == 4
    )
    return cval, dval


def restrict_tensor(adj, z, v):
    card, keep, where = delete_vertex(adj, v)
    out = [[0] * len(keep) for _ in keep]
    for i, j in combinations(keep, 2):
        out[where[i]][where[j]] = out[where[j]][where[i]] = z[i][j]
    return card, tuple(tuple(row) for row in out)


def add_ear(adj, outgoing):
    n = len(adj)
    out = list(adj) + [0]
    for u in range(n):
        if outgoing >> u & 1:
            out[n] |= 1 << u
        else:
            out[u] |= 1 << n
    return tuple(out)


def directed_odd_cycles(adj):
    cycles = []
    for length in range(3, len(adj) + 1, 2):
        for support in combinations(range(len(adj)), length):
            first = support[0]
            for tail in permutations(support[1:]):
                word = (first,) + tail
                if all(adj[word[k]] >> word[(k + 1) % length] & 1
                       for k in range(length)):
                    cycles.append((sum(1 << u for u in support), word))
    return tuple(cycles)


def ocf_packings(adj):
    cycles = directed_odd_cycles(adj)
    result = []

    def visit(start, used, chosen):
        result.append(tuple(chosen))
        for index in range(start, len(cycles)):
            mask, word = cycles[index]
            if not used & mask:
                chosen.append((mask, word))
                visit(index + 1, used | mask, chosen)
                chosen.pop()

    visit(0, 0, [])
    return tuple(result)


def tagged_ocf_atom(adj, edge):
    n = len(adj)
    x = n
    total = moment = 0
    for outgoing in edge:
        for packing in ocf_packings(add_ear(adj, 1 << outgoing)):
            xcycle = None
            used = 0
            for mask, word in packing:
                used |= mask
                if mask >> x & 1:
                    xcycle = word
            if xcycle is None:
                continue
            index = xcycle.index(x)
            neighbours = frozenset((xcycle[index - 1],
                                    xcycle[(index + 1) % len(xcycle)]))
            if neighbours != frozenset(edge):
                continue
            weight = 1 << len(packing)
            outside = (used & ((1 << n) - 1)
                       & ~((1 << edge[0]) | (1 << edge[1]))).bit_count()
            total += weight
            moment += outside * weight
    return total, moment


def exposed_words(adj, edge):
    marked_edge = frozenset(edge)
    words = []
    for word in permutations(range(len(adj))):
        marked = [k for k in range(len(adj) - 1)
                  if frozenset((word[k], word[k + 1])) == marked_edge]
        if len(marked) != 1:
            continue
        if all(
            frozenset((word[k], word[k + 1])) == marked_edge
            or adj[word[k]] >> word[k + 1] & 1
            for k in range(len(adj) - 1)
        ):
            words.append(word)
    return tuple(words)


def marked_identity(adj, v, edge):
    card, keep, where = delete_vertex(adj, v)
    local_edge = (where[edge[0]], where[edge[1]])
    card_words = [tuple(keep[u] for u in word)
                  for word in exposed_words(card, local_edge)]
    redundancy = theft = extensions = 0
    for word in card_words:
        signature = [int(bool(adj[v] >> u & 1)) for u in word]
        drops = sum(signature[k] == 1 and signature[k + 1] == 0
                    for k in range(len(word) - 1))
        mark = next(k for k in range(len(word) - 1)
                    if frozenset((word[k], word[k + 1])) == frozenset(edge))
        stolen = int(signature[mark] == 0 and signature[mark + 1] == 1)
        predicted = 1 + drops - stolen
        literal = int(bool(adj[v] >> word[0] & 1))
        literal += int(bool(adj[word[-1]] >> v & 1))
        literal += sum(
            int(bool(adj[word[k]] >> v & 1)
                and bool(adj[v] >> word[k + 1] & 1))
            for k in range(len(word) - 1) if k != mark
        )
        need(predicted == literal, "marked insertion count")
        redundancy += drops
        theft += stolen
        extensions += literal
    card_set = set(card_words)
    full_words = exposed_words(adj, edge)
    orphans = sum(tuple(u for u in word if u != v) not in card_set
                  for word in full_words)
    need(len(full_words) == extensions + orphans, "deletion fibres")
    need(len(full_words) - len(card_words)
         == redundancy + orphans - theft, "marked identity")


def audit_ocf_and_marked_words():
    ocf_gates = marked_gates = 0
    for n in range(2, 5):
        for code in range(1 << (n * (n - 1) // 2)):
            adj = decode(code, n)
            _, c = capacity(adj)
            card_sum = [[0] * n for _ in range(n)]
            for v in range(n):
                card, keep, where = delete_vertex(adj, v)
                _, cc = capacity(card)
                for i, j in combinations(keep, 2):
                    card_sum[i][j] += cc[where[i]][where[j]]
                    card_sum[j][i] = card_sum[i][j]
            for edge in combinations(range(n), 2):
                total, moment = tagged_ocf_atom(adj, edge)
                need(total == c[edge[0]][edge[1]], "tagged OCF capacity")
                need(moment == (n - 2) * total
                     - card_sum[edge[0]][edge[1]], "OCF support moment")
                ocf_gates += 1
    for n in range(3, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            adj = decode(code, n)
            for v in range(n):
                for edge in combinations(tuple(u for u in range(n) if u != v), 2):
                    marked_identity(adj, v, edge)
                    marked_gates += 1
    return ocf_gates, marked_gates


def audit_n6_monotonicity():
    cache = {code: capacity(decode(code, 5))[1] for code in range(1 << 10)}
    gates = zeros = 0
    minimum = None
    witness = None
    for code in range(1 << 15):
        adj = decode(code, 6)
        _, c = capacity(adj)
        for v in range(6):
            card, keep, where = delete_vertex(adj, v)
            cc = cache[encode(card)]
            for i, j in combinations(keep, 2):
                delta = c[i][j] - cc[where[i]][where[j]]
                gates += 1
                zeros += delta == 0
                if minimum is None or delta < minimum:
                    minimum = delta
                    witness = (code, v, i, j, c[i][j], cc[where[i]][where[j]])
    need(minimum == 0, "capacity deletion monotonicity")
    return gates, zeros, minimum, witness


def audit_restriction_algebra():
    rows = 0
    for n in range(5, 9):
        adj = decode((1 << (n * (n - 1) // 2)) // 3, n)
        z = [[0] * n for _ in range(n)]
        for index, (i, j) in enumerate(combinations(range(n), 2)):
            z[i][j] = z[j][i] = 1 + (17 * index + 5 * n) % 23
        z = tuple(tuple(row) for row in z)
        cval, dval = packet(adj, z)
        deck = [packet(*restrict_tensor(adj, z, v)) for v in range(n)]
        need(sum(row[0] for row in deck) == (n - 3) * cval,
             "C restriction deck")
        need(sum(row[1] for row in deck) == (n - 4) * dval,
             "D restriction deck")
        need(dval > 0 and all(row[1] > 0 for row in deck), "test denominators")
        if n % 2:
            parent_tau = Fraction((n - 3) * cval, 4 * dval)
            child_tau = [Fraction((n - 4) * cv, 2 * dv) for cv, dv in deck]
            factor = Fraction(1, 2)
        else:
            parent_tau = Fraction((n - 3) * cval, 2 * dval)
            child_tau = [Fraction((n - 4) * cv, 4 * dv) for cv, dv in deck]
            factor = Fraction(2, 1)
        weighted = sum(Fraction(deck[v][1], sum(row[1] for row in deck))
                       * child_tau[v] for v in range(n))
        need(parent_tau == factor * weighted, "parity holonomy")
        rows += 1
    return rows


def normalized_deck(adj, z):
    n = len(adj)
    cval, dval = packet(adj, z)
    deck = [packet(*restrict_tensor(adj, z, v)) for v in range(n)]
    need(dval > 0 and all(row[1] > 0 for row in deck), "capacity denominators")
    if n % 2:
        parent = Fraction((n - 3) * cval, 4 * dval)
        children = [Fraction((n - 4) * cv, 2 * dv) for cv, dv in deck]
        factor = Fraction(1, 2)
    else:
        parent = Fraction((n - 3) * cval, 2 * dval)
        children = [Fraction((n - 4) * cv, 4 * dv) for cv, dv in deck]
        factor = Fraction(2, 1)
    total_d = sum(row[1] for row in deck)
    mean = sum(Fraction(deck[v][1], total_d) * children[v]
               for v in range(n))
    need(parent == factor * mean, "named parity identity")
    return parent, mean, max(abs(x) for x in children), sum(abs(x) >= 1 for x in children)


def named_hostiles():
    prime11 = decode(3169369058263173, 11)
    need(is_prime(prime11) and is_strong(prime11), "prime-11 hostile type")
    h11, c11 = capacity(prime11)
    p11 = packet(prime11, c11)
    card10, _, _ = delete_vertex(prime11, 10)
    h10, c10 = capacity(card10)
    p10 = packet(card10, c10)
    need(not is_strong(card10), "prime-11 hostile card")
    need(7 * abs(p10[0]) > 2 * p10[1], "noncentral actual card")

    order12 = (3070, 3644, 3704, 3824, 4064, 4032,
               3970, 3846, 3598, 1024, 2049, 512)
    h12, c12 = capacity(order12)
    p12 = packet(order12, c12)
    deck12 = normalized_deck(order12, c12)
    need((h12, p12) == (27759, (-94387092144, 323484198928)),
         "THM-4133 packet")
    need(deck12[:3] == (
        Fraction(-53092739331, 40435524866),
        Fraction(-53092739331, 80871049732),
        Fraction(9073595176, 12026131621)),
        "THM-4133 restriction ratios")
    return {
        "prime11": (h11, p11, Fraction(2 * abs(p11[0]), p11[1])),
        "bad_card10": (h10, p10, (7 * abs(p10[0]), 2 * p10[1])),
        "order12": (h12, p12, deck12),
    }


def main():
    print("ocf_marked_gates", audit_ocf_and_marked_words())
    print("n6_monotonicity", audit_n6_monotonicity())
    print("restriction_algebra_rows", audit_restriction_algebra())
    print("named_hostiles", named_hostiles())
    print("THM4167_AUDIT_PASS")


if __name__ == "__main__":
    main()
