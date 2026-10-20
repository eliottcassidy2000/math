"""Exact controls for paired exceptional ports; no producer imports.

The universe and each sufficient guard are explicit. Runtime checks use
need(), so optimized Python replays precisely the same mathematical tests.
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import random

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def opposite(arcs):
    return frozenset((v, u) for u, v in arcs)


def reverse_path(arcs, word):
    out = set(arcs)
    for u, v in zip(word, word[1:]):
        if (u, v) in out:
            out.remove((u, v)); out.add((v, u))
        else:
            need((v, u) in out, "complete tournament")
            out.remove((v, u)); out.add((u, v))
    return frozenset(out)


def kmove(arcs, word):
    return opposite(reverse_path(arcs, word))


def directed(arcs, word):
    return all((u, v) in arcs for u, v in zip(word, word[1:]))


def blue(arcs, word):
    rho = dict(zip(word, word[::-1]))
    return all((rho[v], rho[u]) in arcs for u, v in arcs)


def degree(arcs, vertex):
    return sum(u == vertex for u, _ in arcs)


def rows(arcs, n):
    return tuple(sum(1 << v for u, v in arcs if u == i) for i in range(n))


def presentations(n, constraints=None, samples=None, seed=0):
    """Every reflected pair-orbit gets one bit; fixed path arcs remain forward."""
    fixed = {(i, i + 1): 1 for i in range(n - 1)}
    for pair, bit in (constraints or {}).items():
        if pair in fixed and fixed[pair] != bit:
            return
        fixed[pair] = bit
    unseen = set(combinations(range(n), 2)); base = set(); free = []
    while unseen:
        u, v = min(unseen)
        orbit = {(u, v), (n - 1 - v, n - 1 - u)}
        unseen -= orbit
        choices = {fixed[p] for p in orbit if p in fixed}
        if len(choices) > 1:
            return
        if choices:
            bit = choices.pop()
            base.update((a, b) if bit else (b, a) for a, b in orbit)
        else:
            free.append(orbit)
    assignments = product((0, 1), repeat=len(free))
    if samples is not None:
        rng = random.Random(seed)
        assignments = [tuple(rng.randrange(2) for _ in free) for _ in range(samples)]
    for bits in assignments:
        arcs = set(base)
        for orbit, bit in zip(free, bits):
            arcs.update((a, b) if bit else (b, a) for a, b in orbit)
        yield frozenset(arcs)


def layer_constraints(n, h):
    fixed = {(0, j): int(j not in {h, n - 2, n - 1}) for j in range(1, n)}
    fixed[2, n - 2] = 0  # c->x, the A branch.
    return fixed


def check_k(arcs, words):
    n = len(words[0]); target = tuple(range(n)); state = arcs; excesses = []
    need(directed(arcs, target) and blue(arcs, target), "source standard-blue")
    for word in words:
        need(len(word) == n and set(word) == set(target), "Hamiltonian vertex set")
        need(directed(state, word), "actual K word directed")
        need(not blue(state, word), "actual K word black")
        excesses.append(degree(state, word[0]) + degree(state, word[-1]) - n + 1)
        state = kmove(state, word)
    need(state == kmove(arcs, target), "literal labelled K target")
    return excesses


def check_r(arcs, words):
    state = arcs; sums = []
    for word in words:
        need(directed(state, word), "actual R word directed")
        need(not blue(state, word), "actual R word black")
        sums.append(degree(state, word[0]) + degree(state, word[-1]))
        state = reverse_path(state, word)
    need(state == reverse_path(arcs, tuple(range(len(words[0])))), "literal R target")
    check_k(arcs, [w if i % 2 == 0 else w[::-1] for i, w in enumerate(words)])
    return sums


def family(middle, r):
    A = (0, 1, 2); M = tuple(range(3, r + 3)); Z = (r + 3, r + 4, r + 5)
    return frozenset(set(combinations(A, 2)) | set(combinations(Z, 2)) |
                     {(u + 3, v + 3) for u, v in middle} |
                     {(u, v) for u in A for v in M} |
                     {(u, v) for u in M for v in Z} |
                     {(u, v) for u in Z for v in A})


def end_pair_words(n):
    W = tuple(range(3, n - 3)); h, c, z = n - 3, n - 2, n - 1
    H = (2,) + W[::-1] + (c, h, z)
    return [(0, 2) + W + (h, c, z, 1), (0, 1) + H, (1,) + H[::-1] + (0,)]


def check_family():
    counts = Counter(); larger = []; example = None
    for r in range(2, 7):
        for middle in presentations(r):
            arcs = family(middle, r); n = r + 6
            sums = check_r(arcs, end_pair_words(n))
            need(sums == [2 * r + 3, r + 4, 2 * r + 2], "family sums")
            need((degree(arcs, 0), degree(arcs, n - 1)) == (n - 4, 3), "endpoint layer")
            cycle = tuple(range(3)) + tuple(range(r + 3, r + 6)) + tuple(range(3, r + 3)) + (0,)
            need(directed(kmove(arcs, tuple(range(n))), cycle), "strong target cycle")
            counts[r] += 1
            if r == 4 and all((i, j) in middle for i, j in combinations(range(r), 2)):
                example = dict(rows=rows(arcs, n), R_words=end_pair_words(n), endpoint_sums=sums)
    for r in (7, 8, 9, 12, 17, 24, 31):
        for middle in presentations(r, samples=3, seed=72000 + r):
            check_r(family(middle, r), end_pair_words(r + 6))
        larger.append(dict(middle_order=r, presentations=3))
    small = family(frozenset(), 1)
    words = [(0, 2, 3, 4, 5, 6, 1), (3, 0, 2, 4, 5, 6, 1), (3, 0, 1, 2, 4, 5, 6)]
    excesses = check_k(small, words)
    state = kmove(small, words[0])
    need((3, 4) in state and (4, 1) not in state, "small balanced anti hostile")
    return dict(exhaustive_blue_middle_counts=dict(counts), large_controls=larger,
                ten_vertex_343=example,
                r1=dict(rows=rows(small, 7), K_words=words, excesses=excesses))


def a_words(n, length):
    D = tuple(range(3, 3 + length)); U = tuple(range(3 + length, n - 1)); z = n - 1
    return [(0, 2) + D + U + (z, 1), D + (0, 2) + U + (z, 1),
            D + (0, 1, 2) + U + (z,)]


def check_nonmodule_and_a():
    results = {}; digests = {}
    for name in ("end_pair", "A1_full", "A2_singleton", "A2_odd_prefix"):
        records = []; checks = 0
        for n in ((8, 9, 10, 14, 25) if name == "end_pair" else (11, 12, 15, 20, 31)):
            lengths = (None,) if name == "end_pair" else ((n - 5,) if name == "A1_full" else ((1,) if name == "A2_singleton" else (3, 5, 7)))
            for length in lengths:
                if name == "A2_odd_prefix" and n < 2 * length + 6:
                    continue
                for h in ((n - 3,) if name == "end_pair" else range(2, n - 2)):
                    if name == "A1_full" and h in {2, 3, 4, n - 8, n - 3}:
                        continue
                    if name == "A2_singleton" and h in {2, 3}:
                        continue
                    if name == "A2_odd_prefix" and h in {2, 4, length - 1, length + 2}:
                        continue
                    fixed = layer_constraints(n, h)
                    if name == "end_pair":
                        fixed[1, n - 4] = 1; fixed[2, n - 4] = 1
                    elif name in {"A2_singleton", "A2_odd_prefix"}:
                        fixed[2, 3 + length] = 0
                    count = 0
                    for arcs in presentations(n, fixed, samples=8, seed=100000 * n + 100 * h + (length or 0)):
                        need((degree(arcs, 0), degree(arcs, n - 1)) == (n - 4, 3), "guarded endpoint layer")
                        if name == "end_pair":
                            check_r(arcs, end_pair_words(n))
                        else:
                            check_k(arcs, a_words(n, length))
                        records.append(rows(arcs, n)); count += 1
                    need(count == 8, "nonempty compatible guarded universe")
                    checks += count
        results[name] = checks
        digests[name] = hashlib.sha256(json.dumps(records).encode()).hexdigest()
    return dict(guarded_presentations_checked=results, row_bank_sha256=digests,
                note="Eight deterministic reflected-chord assignments per listed parameter; not a tournament census.")


def check_reflection_order():
    counts = 0
    for length in (1, 3, 5, 7, 9, 11):
        for n in range(2 * length + 6, 2 * length + 16):
            word = a_words(n, length)[2]
            reversal = dict(zip(word, word[::-1]))
            g = [reversal[n - 1 - v] for v in range(n)]
            unseen = set(range(n)); sizes = []
            while unseen:
                v = min(unseen); orbit = []; u = v
                while u not in orbit:
                    orbit.append(u); unseen.remove(u); u = g[u]
                sizes.append(len(orbit))
            need(any(s % 2 == 0 for s in sizes), "reflection product has even order")
            counts += 1
    return dict(disjoint_reflection_commutators=counts)


def check_naive_lift():
    n = 5; middle = (1, 2, 3)
    macro = frozenset(set(combinations(middle, 2)) | {(0, v) for v in middle} |
                      {(v, 4) for v in middle} | {(4, 0)})
    words = [(1, 2, 3, 4, 0), (4, 2, 3, 1, 0), (0, 4, 2, 3, 1)]
    check_k(macro, words)
    blocks = {0: (0, 1, 2), 1: (3,), 2: (4,), 3: (5,), 4: (6, 7, 8)}
    arcs = frozenset({(x, y) for u, v in macro for x in blocks[u] for y in blocks[v]} |
                     set(combinations(blocks[0], 2)) | set(combinations(blocks[4], 2)))
    expanded = [tuple(x for v in w for x in blocks[v]) for w in words]
    state = arcs; legal = []; bad = []
    for word in expanded:
        legal.append(directed(state, word))
        bad.append([(u, v) for u, v in zip(word, word[1:]) if (u, v) not in state])
        state = kmove(state, word)
    need(legal == [True, True, False], "TT3 naive prefix lift hostile")
    need(state != kmove(arcs, tuple(range(9))), "naive expanded target wrong")
    check_r(arcs, end_pair_words(9))
    return dict(expanded_K_words=expanded, directed=legal, wrong_arcs=bad,
                literal_target_matches=False, repaired_by_theorem_M=True)


def final_endpair_words(n, x_to_y):
    y, c, z = n - 3, n - 2, n - 1; W = tuple(range(3, y))
    Q1 = (0, 2) + W + (y, c, z, 1)
    if x_to_y:
        return [Q1, (z,) + W + (y, 2, 1, 0, c), (1, z) + W + (y, 2, 0, c)]
    return [Q1, W + (0, 2, y, c, z, 1), W + (0, 1, 2, y, c, z)]


def check_final_endpair():
    receipt = []
    for x_to_y in (0, 1):
        orders = (11, 12, 13, 17, 24, 31) if not x_to_y else (9, 10, 11, 12, 13, 17, 24, 31)
        for n in orders:
            fixed = layer_constraints(n, n - 3); fixed[2, n - 3] = x_to_y
            sample = None if n in (9, 10) else 64
            count = 0; middles = Counter()
            for arcs in presentations(n, fixed, samples=sample, seed=n * 100 + x_to_y):
                words = final_endpair_words(n, x_to_y); check_k(arcs, words)
                S = kmove(arcs, words[0]); B = kmove(arcs, tuple(range(n)))
                if x_to_y:
                    need((0, 5) not in S and (2, 3) in S, "new Q2 fixed-arc contradiction")
                    need((0, 3) not in B and (2, n - 1) in B, "new Q3 fixed-arc contradiction")
                else:
                    for state, word, vertex in [(S, words[1], 4), (B, words[2], 3)]:
                        rho = dict(zip(word, word[::-1]))
                        need(rho[vertex] == n - 1, "variable cut witness paired with z")
                        if n == 11:
                            need(rho[0] == 0 and degree(state, 0) == 4, "n11 fixed score mismatch")
                        else:
                            need(rho[0] == 8 and (vertex, 0) in state and (8, n - 1) not in state,
                                 "variable cut exact-neighborhood contradiction")
                middles["open" if (1, n - 2) in arcs else "closed"] += 1; count += 1
            receipt.append(dict(n=n, x_to_y=bool(x_to_y), checked=count,
                                universe="complete fixed-blue presentations" if sample is None else "64 deterministic mirrored chord assignments",
                                middle_directions=dict(middles)))
    return receipt


def check_sidecar_hostile():
    n = 8; fixed = layer_constraints(n, n - 3)
    fixed[2, 4] = 0; fixed[2, 5] = 0
    arcs = next(presentations(n, fixed)); right = tuple(range(2, n))
    changed = reverse_path(arcs, right)
    need(directed(arcs, tuple(range(n))) and blue(arcs, tuple(range(n))), "sidecar hostile actual blue source")
    need(all((2, v) not in changed for v in right if v != 2), "prescribed start becomes sink")
    return dict(rows=rows(arcs, n), reversed_right_block=right, start=2,
                start_outdegree_inside_reversed_block=0,
                conclusion="This particular endpoint-preserving right-block path does not exist.")


def check_A3():
    count = 0
    for n in range(11, 201):
        word = a_words(n, n - 5)[2]; rho = dict(zip(word, word[::-1]))
        g = {v: rho[n - 1 - v] for v in range(n)}
        formula = {0: 3, 1: 4, **{i: i + 6 for i in range(2, n - 8)},
                   n - 8: 0, n - 7: 1, n - 6: 2, n - 5: n - 2,
                   n - 4: n - 1, n - 3: 5, n - 2: 6, n - 1: 7}
        need(g == formula, "A3 explicit permutation")
        orbit = []; v = 0
        while v not in orbit:
            orbit.append(v); v = g[v]
        residue = n % 6
        expected = ((n, {1, n - 2, n - 1}) if residue in (1, 3) else
                    (n // 2, {1}) if residue in (0, 4) else
                    ((n + 1) // 3, {n - 2}) if residue == 2 else
                    ((n + 1) // 6, set()))
        need((len(orbit), set(orbit) & {1, n - 2, n - 1}) == expected, "A3 orbit-capacity table")
        count += 1
    witness = (32758, 9348, 93336, 130355, 73762, 91334, 126366, 65816,
               75318, 118014, 67728, 78258, 92150, 18048, 101778, 81335, 8195)
    n = 17; arcs = frozenset((u, v) for u in range(n) for v in range(n) if witness[u] & (1 << v))
    need(len(arcs) == n * (n - 1) // 2 and all(((u, v) in arcs) != ((v, u) in arcs) for u, v in combinations(range(n), 2)), "A3 witness tournament")
    need(all((u, u) not in arcs for u in range(n)), "no witness loops")
    P = tuple(range(n)); B = kmove(arcs, P); word = a_words(n, n - 5)[2]
    need(directed(arcs, P) and blue(arcs, P), "A3 witness source blue")
    need((degree(arcs, 0), degree(arcs, n - 1)) == (13, 3), "A3 witness endpoint degrees")
    need({u for u, v in arcs if v == 0} == {3, 15, 16} and (15, 2) in arcs, "A3 witness exceptional ports and chord")
    need(directed(B, word) and blue(B, word), "A3 order17 final word genuinely blue")
    alternative = [tuple(range(13)) + (16, 13, 14, 15),
                   tuple(range(1, 13)) + (16, 13, 14, 0, 15),
                   tuple(range(1, 15)) + (0, 15, 16)]
    alternative_excesses = check_k(arcs, alternative)
    return dict(permutation_orders_checked=[11, 200], order_count=count,
                source_of_order17_witness="root lane", witness_rows=witness,
                blue_final_K_word=word, movable_cut_K_words=alternative,
                movable_cut_excesses=alternative_excesses,
                conclusion="Fixed-word hostile with independently replayed three-black alternative.")


def main():
    receipt = dict(status="PASS", scope="Paired-port sufficient repairs, not full endpoint-three coverage",
                   family=check_family(), nonmodule_and_A=check_nonmodule_and_a(),
                   reflection_order=check_reflection_order(), naive_lift=check_naive_lift(),
                   final_endpair=check_final_endpair(), missing_sidecar=check_sidecar_hostile(),
                   A3=check_A3())
    (HERE / "lift_check.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print("PASS: TT3 family, nonmodule end-pair repair, whole/odd-prefix A repairs, and naive lift hostile")
    print(json.dumps(receipt["nonmodule_and_A"]["guarded_presentations_checked"], sort_keys=True))


if __name__ == "__main__":
    main()
