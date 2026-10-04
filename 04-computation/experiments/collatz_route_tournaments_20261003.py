"""Exact finite audits of intrinsically decodable Collatz route tournaments.

Run: python3 04-computation/experiments/collatz_route_tournaments_20261003.py
No packages, external data, or inherited filters.  These are certificates for
the already-generated basin of 1, not a proof that every integer enters it.
Graphs use integer adjacency bit masks and the ordinary map n/2 or 3n+1,
stopped at 1.  SCC decoding never receives the construction boundaries.
"""

from itertools import product


def v2(n):
    if n <= 0:
        raise ValueError("valuation requires a positive integer")
    return (n & -n).bit_length() - 1


def minimum_exponent(y):
    if y <= 0 or y % 2 == 0 or y % 3 == 0:
        raise ValueError("an odd predecessor needs an odd target prime to 3")
    return 2 if y % 3 == 1 else 1


def decode_odd_word(word):
    """Decode compressed inverse choices, excluding padding around odd 1."""
    y = 1
    for j in reversed(tuple(word)):
        if not isinstance(j, int) or j < 0:
            raise ValueError("inverse choices are nonnegative integers")
        k = minimum_exponent(y) + 2 * j
        numerator = (y << k) - 1
        assert numerator % 3 == 0
        y = numerator // 3
        assert y > 0 and y % 2 == 1
        if y == 1:
            raise ValueError("the reduced route stops at its first arrival at 1")
    return y


def encode_odd_word(n, max_steps=10000):
    """Return a certified first-hit word; the explicit bound can fail."""
    if n <= 0 or n % 2 == 0:
        raise ValueError("source must be positive and odd")
    word = []
    while n != 1:
        if len(word) >= max_steps:
            raise ValueError("no arrival at 1 within the supplied step bound")
        z = 3 * n + 1
        k = v2(z)
        y = z >> k
        k0 = minimum_exponent(y)
        assert k >= k0 and (k - k0) % 2 == 0
        word.append((k - k0) // 2)
        n = y
    return tuple(word)


def route_state(n):
    e = v2(n)
    return e, encode_odd_word(n >> e)


def regular_cyclic(m):
    if m < 3 or m % 2 == 0:
        raise ValueError("cyclic strong blocks have odd size at least 3")
    return tuple(sum(1 << ((i + d) % m) for d in range(1, (m + 1) // 2))
                 for i in range(m))


def ordinal_sum(graphs):
    n = sum(map(len, graphs))
    rows, offset = [], 0
    for graph in graphs:
        later = ((1 << n) - 1) ^ ((1 << (offset + len(graph))) - 1)
        rows.extend((row << offset) | later for row in graph)
        offset += len(graph)
    return tuple(rows)


def clone(graph, size=2):
    """Lexicographic substitution by a transitive size-vertex tournament."""
    if size < 1:
        raise ValueError("clone size must be positive")
    rows = []
    fibre = (1 << size) - 1
    for u, row in enumerate(graph):
        external = 0
        for v in bits(row):
            external |= fibre << (v * size)
        for r in range(size):
            internal = (fibre ^ ((1 << (r + 1)) - 1)) << (u * size)
            rows.append(external | internal)
    return tuple(rows)


def route_graph(word, exponent=0):
    decode_odd_word(word)  # The graph family contains only arithmetic certificates.
    base = ordinal_sum([regular_cyclic(2 * j + 3) for j in word] + [(0,)])
    return clone(base, 1 << exponent)


def bits(mask):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def check_tournament(graph):
    n = len(graph)
    for u, row in enumerate(graph):
        assert 0 <= row < (1 << n) and not ((row >> u) & 1)
        for v in range(u):
            assert ((row >> v) & 1) + ((graph[v] >> u) & 1) == 1


def reach(graph, start):
    seen, pending = 0, 1 << start
    while pending:
        bit = pending & -pending
        pending ^= bit
        if seen & bit:
            continue
        seen |= bit
        pending |= graph[bit.bit_length() - 1] & ~seen
    return seen


def ordered_sccs(graph):
    """Compute SCCs from adjacency, then their intrinsic tournament order."""
    reverse = [0] * len(graph)
    for u, row in enumerate(graph):
        for v in bits(row):
            reverse[v] |= 1 << u
    remaining = (1 << len(graph)) - 1
    components = []
    while remaining:
        v = (remaining & -remaining).bit_length() - 1
        component = reach(graph, v) & reach(reverse, v)
        components.append(tuple(bits(component)))
        remaining &= ~component
    unordered = tuple(components)
    components.sort(key=lambda c: -sum((graph[c[0]] >> d[0]) & 1
                                      for d in unordered))
    for i, component in enumerate(components):
        for later in components[i + 1:]:
            assert all((graph[u] >> v) & 1 for u in component for v in later)
    return tuple(components)


def decode_graph(graph):
    """Read a generated carrier by SCCs, without any supplied block labels.

    This decodes routes, but is not a complete recognition test for arbitrary
    tournaments isomorphic to the chosen cyclic-block family.
    """
    components = ordered_sccs(graph)
    trailing = 0
    for component in reversed(components):
        if len(component) != 1:
            break
        trailing += 1
    if not trailing or trailing & (trailing - 1):
        raise ValueError("home must be a nonempty power-of-two singleton suffix")
    word = []
    for component in components[:-trailing]:
        if len(component) % trailing:
            raise ValueError("SCC size is not divisible by the clone scale")
        q = len(component) // trailing
        if q < 3 or q % 2 == 0:
            raise ValueError("odd strong-block size is invalid")
        word.append((q - 3) // 2)
    word = tuple(word)
    return trailing.bit_length() - 1, word, trailing * decode_odd_word(word)


def induced(graph, vertices):
    return tuple(sum(((graph[u] >> v) & 1) << j
                     for j, v in enumerate(vertices)) for u in vertices)


def twin_fibres(graph):
    """Connected components of the identical-external-neighborhood relation."""
    parent = list(range(len(graph)))

    def root(v):
        while parent[v] != v:
            parent[v] = parent[parent[v]]
            v = parent[v]
        return v

    for u in range(len(graph)):
        for v in range(u):
            outside = ~((1 << u) | (1 << v))
            if not ((graph[u] ^ graph[v]) & outside):
                parent[root(u)] = root(v)
    groups = {}
    for v in range(len(graph)):
        groups.setdefault(root(v), []).append(v)
    ordered = []
    for group in sorted(groups.values(), key=min):
        mask = sum(1 << v for v in group)
        group.sort(key=lambda v: -(graph[v] & mask).bit_count())
        assert [(graph[v] & mask).bit_count() for v in group] == list(
            range(len(group) - 1, -1, -1))
        ordered.append(tuple(group))
    return tuple(ordered)


def halve_graph(graph):
    e, _, _ = decode_graph(graph)
    if e == 0:
        raise ValueError("halving requires a positive clone exponent")
    fibres = twin_fibres(graph)
    assert all(len(fibre) == 1 << e for fibre in fibres)
    pairs = [fibre[i:i + 2] for fibre in fibres for i in range(0, len(fibre), 2)]
    for i, pair in enumerate(pairs):
        for other in pairs[i + 1:]:
            assert len({(graph[u] >> v) & 1 for u in pair for v in other}) == 1
    return induced(graph, [pair[0] for pair in pairs])


def ordinary_step_graph(graph, vertex_bound=256):
    """Stopped ordinary Collatz rewrite; reject oversized audit constructions."""
    e, word, n = decode_graph(graph)
    if n == 1:
        return graph  # First-hit certificates use a stopped root.
    if e:
        return halve_graph(graph)
    first = set(ordered_sccs(graph)[0])
    tail = induced(graph, [v for v in range(len(graph)) if v not in first])
    _, tail_word, y = decode_graph(tail)
    assert tail_word == word[1:]
    k = minimum_exponent(y) + 2 * word[0]
    if len(tail) * (1 << k) > vertex_bound:
        raise ValueError("graph rewrite exceeds this finite audit's vertex bound")
    return clone(tail, 1 << k)


def main():
    valid, rejected = [], 0
    for length in range(1, 6):
        for word in product(range(5), repeat=length):
            try:
                n = decode_odd_word(word)
            except ValueError:
                rejected += 1
                continue
            assert encode_odd_word(n) == word
            assert decode_odd_word((word[0] + 1,) + word[1:]) == 4 * n + 1
            valid.append((word, n))
    assert decode_odd_word(()) == 1 and encode_odd_word(1) == ()
    print("Universe: j in {0,1,2,3,4}, lengths 1..5; no inherited filters")
    print(f"Inverse grammar: {len(valid)} valid, {rejected} rejected, plus the empty root")
    print(f"Forward replay and safe first-block +2 / source 4n+1: {len(valid)} each")
    macro_positions = phase_checks = 0
    for word, n in valid:
        values, valuations, x = [n], [], n
        while x != 1:
            z = 3 * x + 1
            k = v2(z)
            valuations.append(k)
            x = z >> k
            values.append(x)
        A = S = 0
        for p, old_k in enumerate(valuations):
            period = 3 ** p
            R = 4 ** period
            numerator = (R - 1) * ((1 << A) + 3 * S)
            denominator = 3 ** (p + 1)
            assert numerator % denominator == 0
            K = numerator // denominator
            changed = word[:p] + (word[p] + period,) + word[p + 1:]
            new_n = decode_odd_word(changed)
            assert new_n == R * n + K and new_n > n
            new_valuations, y = [], new_n
            while y != 1:
                z = 3 * y + 1
                k = v2(z)
                new_valuations.append(k)
                y = z >> k
            expected = valuations[:]
            expected[p] += 2 * period
            assert new_valuations == expected
            assert encode_odd_word(new_n) == changed
            # Every smaller positive increment fails to retain this exact head.
            h = values[p]
            for delta in range(1, period + 1):
                h_delta = 4 ** delta * h + (4 ** delta - 1) // 3
                assert (((1 << A) * h_delta - S) % period == 0) == (delta == period)
                phase_checks += 1
            macro_positions += 1
            S = 3 * S + (1 << A)
            A += old_k
    print(f"General internal macro: {macro_positions} word/position pairs; increment j_p by 3^p")
    print(f"Exact affine recurrence and all valuations pass; {phase_checks} minimal-period controls")
    for m in range(3, 32, 2):
        a = (m - 1) // 2
        embedding = [i if i <= a else i + 1 for i in range(m)]
        assert induced(regular_cyclic(m + 2), embedding) == regular_cyclic(m)
    print("Cyclic +2 is an induced extension: every odd block size 3..31")

    multisets, witness = {}, None
    for word, n in valid:
        key = tuple(sorted(word))
        if key in multisets and multisets[key][0] != word:
            witness = multisets[key], (word, n)
            break
        multisets[key] = word, n
    assert witness and witness[0][1] != witness[1][1]
    print("Same multiset, different certified order/source:", witness)
    assert decode_odd_word((0, 1)) == 3
    try:
        decode_odd_word((0, 2))
    except ValueError as exc:
        print("Internal hostile (0,1)->(0,2): 3->5->1 becomes tail 21;", str(exc))
    else:
        raise AssertionError("internal block mutation was incorrectly accepted")
    assert decode_odd_word((0, 3)) == 113
    assert v2(3 * 3 + 1) == 1 and v2(3 * 113 + 1) == 2
    print("Minimality scope: (0,1)->(0,3) is valid at delta=2<3, but head valuation changes 1->2")

    # A phase-guarded internal macro really is lawful, unlike a bare internal +2.
    macro_graphs, macro_sources = [], []
    for t in range(4):
        numerator = 4 ** (9 * t + 6) - 19
        assert numerator % 27 == 0
        n = numerator // 27
        word = (0, 0, 9 * t + 4)
        assert decode_odd_word(word) == n and encode_odd_word(n) == word
        values, valuations, x = [n], [], n
        while x != 1:
            z = 3 * x + 1
            k = v2(z)
            valuations.append(k)
            x = z >> k
            values.append(x)
        assert valuations == [1, 1, 18 * t + 10]
        assert values[1] == (3 * n + 1) // 2
        assert values[2] == (4 ** (9 * t + 5) - 1) // 3
        graph = route_graph(word)
        assert decode_graph(graph) == (0, word, n)
        assert [len(c) for c in ordered_sccs(graph)] == [3, 3, 18 * t + 11, 1]
        assert len(graph) == 18 * t + 18
        assert sum(k + 1 for k in valuations) == 18 * t + 15
        macro_graphs.append(graph)
        macro_sources.append(n)
    for t in range(3):
        assert macro_sources[t + 1] == 262144 * macro_sources[t] + 184471
        old, new = macro_graphs[t:t + 2]
        m = 18 * t + 11
        embedding = list(range(m))
        for _ in range(9):
            a = (m - 1) // 2
            embedding = [i if i <= a else i + 1 for i in embedding]
            m += 2
        whole = list(range(6)) + [6 + i for i in embedding] + [len(new) - 1]
        assert induced(new, whole) == old
    print("Guarded internal +18 macro: t=0..3, word (0,0,9t+4), valuations (1,1,18t+10)")
    print("Macro sources:", macro_sources)
    print("Macro: recurrence 262144n+184471; SCC sizes (3,3,18t+11,1); vertices 18t+18")
    print("Macro: 3 induced old-edge-preserving extensions; ordinary remaining steps 18t+15")

    samples = [(), (1,), (0, 1), (2,), (2, 1), (1, 2, 1)]
    carriers = halved = 0
    for word in samples:
        source = decode_odd_word(word)
        for e in (0, 1, 2):
            graph = route_graph(word, e)
            check_tournament(graph)
            assert decode_graph(graph) == (e, word, (1 << e) * source)
            # A reversed vertex labelling is a hostile control for boundary-free decoding.
            relabelled = induced(graph, list(reversed(range(len(graph)))))
            assert decode_graph(relabelled) == (e, word, (1 << e) * source)
            assert len(graph) == (1 << e) * (1 + sum(2 * j + 3 for j in word))
            if e == 0 and word:
                larger = route_graph((word[0] + 1,) + word[1:])
                m = 2 * word[0] + 3
                a = (m - 1) // 2
                embedding = [i if i <= a else i + 1 for i in range(m)]
                embedding.extend(range(m + 2, len(graph) + 2))
                assert induced(larger, embedding) == graph
            carriers += 1
            if e:
                assert halve_graph(graph) == route_graph(word, e - 1)
                assert decode_graph(halve_graph(relabelled))[2] == (1 << (e - 1)) * source
                halved += 1
    print(f"Actual tournaments: {carriers}; SCC decoding also under reversed vertex labels")
    print(f"Intrinsic twin-fibre halving: {halved}, also under reversed labels")

    odd_steps = 0
    for word, n in valid:
        tail_word = word[1:]
        y = decode_odd_word(tail_word)
        k = minimum_exponent(y) + 2 * word[0]
        size = (1 << k) * (1 + sum(2 * j + 3 for j in tail_word))
        if size > 128:
            continue
        rewritten = ordinary_step_graph(route_graph(word), vertex_bound=128)
        assert rewritten == route_graph(tail_word, k)
        assert decode_graph(rewritten)[2] == 3 * n + 1
        original = route_graph(word)
        relabelled = induced(original, list(reversed(range(len(original)))))
        assert decode_graph(ordinary_step_graph(relabelled, vertex_bound=128))[2] == 3 * n + 1
        odd_steps += 1
        if odd_steps == 30:
            break
    assert ordinary_step_graph((0,)) == (0,)
    assert decode_graph(route_graph((), 2))[2] == 4
    assert decode_graph(ordinary_step_graph(route_graph((), 2)))[2] == 2
    assert decode_graph(ordinary_step_graph(route_graph((), 1)))[2] == 1
    print(f"Ordinary odd-step rewrites: {odd_steps}, also reversed labels; outputs <=128 vertices")
    print("Root convention: stop at 1; the explicit even chain 4->2->1 also passes")
    print("Size: 2^e*(1+sum(2j+3)); n=3 has 9 vertices, n=12 has 36")
    print("Compression: retain ordered SCC indices and binary clone exponent; no expanded graph needed")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
