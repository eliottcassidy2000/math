"""Independent exact controls for mirrored blow-ups and port-aware K lifts.

No producer modules, graph-isomorphism package, or large tournament census.
Every guard survives python -O.
"""
from collections import Counter
from itertools import combinations, product
import hashlib
import json
from pathlib import Path
import random


HERE = Path(__file__).resolve().parent


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def opposite(arcs):
    return frozenset((v, u) for u, v in arcs)


def path_edges(word):
    return {frozenset((u, v)) for u, v in zip(word, word[1:])}


def toggle(arcs, word):
    out = set(arcs)
    for u, v in zip(word, word[1:]):
        if (u, v) in out:
            out.remove((u, v)); out.add((v, u))
        else:
            need((v, u) in out, "a complete tournament pair")
            out.remove((v, u)); out.add((u, v))
    return frozenset(out)


def kmove(arcs, word):
    return opposite(toggle(arcs, word))


def directed(arcs, word):
    return all((u, v) in arcs for u, v in zip(word, word[1:]))


def blue(arcs, word):
    reversal = dict(zip(word, word[::-1]))
    return all((reversal[v], reversal[u]) in arcs for u, v in arcs)


def scores(arcs, n):
    out = [0] * n
    for u, _ in arcs:
        out[u] += 1
    return out


def standard_blue(n):
    unseen = set(combinations(range(n), 2))
    fixed, free = [], []
    while unseen:
        u, v = min(unseen)
        orbit = {(u, v), (n - 1 - v, n - 1 - u)}
        unseen -= orbit
        (fixed if any(b == a + 1 for a, b in orbit) else free).append(orbit)
    base = set().union(*fixed)
    for bits in product((0, 1), repeat=len(free)):
        arcs = set(base)
        for orbit, bit in zip(free, bits):
            arcs.update((u, v) if bit else (v, u) for u, v in orbit)
        yield frozenset(arcs)


def family(middle, r):
    left, mid, right = (0, 1), tuple(range(2, r + 2)), (r + 2, r + 3)
    return frozenset({left, right} |
                     {(u, v) for u in left for v in mid} |
                     {(u, v) for u in mid for v in right} |
                     {(u, v) for u in right for v in left} |
                     {(u + 2, v + 2) for u, v in middle})


def check_r_ledger(arcs, words, target_word):
    state, kstate = arcs, arcs
    endpoint_sums = []
    for i, word in enumerate(words):
        need(directed(state, word), "R word directed")
        need(not blue(state, word), "R word black")
        d = scores(state, len(word))
        endpoint_sums.append(d[word[0]] + d[word[-1]])
        state = toggle(state, word)
        kword = word if i % 2 == 0 else word[::-1]
        need(directed(kstate, kword) and not blue(kstate, kword), "alternating K word black")
        kstate = kmove(kstate, kword)
        need(kstate == (opposite(state) if i % 2 == 0 else state), "K/R complement phase")
    need(state == toggle(arcs, target_word), "literal R endpoint")
    need(kstate == kmove(arcs, target_word), "literal K endpoint")
    return endpoint_sums


def family_ledger(r):
    mid = tuple(range(2, r + 2)); z0, z1 = r + 2, r + 3
    repaired = (mid[0], z0) + mid[:0:-1] + (z1,)
    return [(0,) + mid + (z0, z1, 1), (0, 1) + repaired,
            (1,) + repaired[::-1] + (0,)]


def check_family():
    cases = Counter()
    for r in range(3, 7):
        for middle in standard_blue(r):
            arcs = family(middle, r); n = r + 4; target = tuple(range(n))
            need(blue(arcs, target), "marked family path blue")
            words = family_ledger(r)
            sums = check_r_ledger(arcs, words, target)
            need(sums == [2 * r + 1, r + 2, 2 * r], "family endpoint sums")
            before_last = toggle(toggle(arcs, words[0]), words[1])
            need((r + 3, 0) in before_last and (2, 1) in before_last,
                 "explicit final anti-involution hostile")
            cases[r] += 1
    # Large middles have arbitrary mirrored chord data, not just transitive interiors.
    random_source = random.Random(602020)
    larger = []
    for r in (7, 8, 9, 12, 17, 24, 31):
        unseen = set(combinations(range(r), 2)); arcs = set()
        while unseen:
            u, v = min(unseen); orbit = {(u, v), (r - 1 - v, r - 1 - u)}
            unseen -= orbit
            bit = 1 if any(b == a + 1 for a, b in orbit) else random_source.randrange(2)
            arcs.update((a, b) if bit else (b, a) for a, b in orbit)
        whole = family(frozenset(arcs), r)
        check_r_ledger(whole, family_ledger(r), tuple(range(r + 4)))
        larger.append(r)
    # r=2 has a different literal three-black certificate.
    r2 = family(frozenset({(0, 1)}), 2)
    small_words = [(0, 2, 5, 1, 3, 4), (4, 0, 1, 5, 2, 3), (3, 1, 2, 0, 4, 5)]
    small_sums = check_r_ledger(r2, small_words, tuple(range(6)))
    # r=1 is a two-black detour in the ordinary-isomorphism quotient.
    r1 = family(frozenset(), 1); partner = toggle(r1, tuple(range(5)))
    left, right, pi = (1, 2, 3, 4, 0), (0, 2, 4, 3, 1), (4, 0, 2, 3, 1)
    need(directed(r1, left) and not blue(r1, left), "r1 source black word")
    need(directed(partner, right) and not blue(partner, right), "r1 target black word")
    neighbor = toggle(r1, left)
    need(frozenset((pi[u], pi[v]) for u, v in neighbor) == toggle(partner, right),
         "r1 exact common-neighbor isomorphism")
    return dict(exhaustive_blue_middle_orders=dict(cases),
                arbitrary_mirrored_large_middles=larger,
                r2=dict(R_words=small_words, endpoint_sums=small_sums),
                r1=dict(source_R=left, target_R=right, isomorphism=pi))


def expand(arcs, sizes):
    blocks = {}; offset = 0
    for v, size in enumerate(sizes):
        blocks[v] = tuple(range(offset, offset + size)); offset += size
    result = set()
    for block in blocks.values():
        result.update(combinations(block, 2))  # A fixed transitive internal path suffices.
    result.update((x, y) for u, v in arcs for x in blocks[u] for y in blocks[v])
    return frozenset(result), blocks


def expanded_word(word, blocks):
    return tuple(v for u in word for v in blocks[u])


def check_port_criterion():
    checks = 0; hostile = None
    for r in range(3, 7):
        # Cyclic singleton extension of an open transitive middle.
        n = r + 2; mid = tuple(range(1, r + 1)); z = r + 1
        macro = frozenset({(0, v) for v in mid} | {(v, z) for v in mid} |
                          {(z, 0)} | set(combinations(mid, 2)))
        target = tuple(range(n)); tail = mid[1:]
        ledger = [mid + (z, 0), (z,) + tail + (mid[0], 0), (0, z) + tail + (mid[0],)]
        source = macro
        for word in ledger:
            need(directed(source, word) and not blue(source, word), "inherited macro ledger")
            source = kmove(source, word)
        need(source == kmove(macro, target), "inherited macro endpoint")
        for first_size, last_size in product(range(1, 4), repeat=2):
            sizes = [first_size] + [1] * r + [last_size]
            expanded, blocks = expand(macro, sizes)
            state = expanded; uses = Counter(); failures = []
            for step, word in enumerate(ledger, 1):
                guard = all(uses[v, u] % 2 == 0 for u, v in zip(word, word[1:])
                            if sizes[u] > 1 or sizes[v] > 1)
                lifted = expanded_word(word, blocks)
                actual = directed(state, lifted)
                need(guard == actual, "prefix opposite-port parity criterion")
                failures.append(None if actual else {
                    "step": step,
                    "wrongly_directed_joins": [(u, v) for u, v in zip(lifted, lifted[1:])
                                              if (u, v) not in state]})
                state = kmove(state, lifted)  # Formal XOR remains defined even after an illegal word.
                uses.update(zip(word, word[1:]))
            target_uses = Counter(zip(target, target[1:]))
            parity = all(uses[u, v] % 2 == target_uses[u, v]
                         for u in range(n) for v in range(n) if u != v
                         and (sizes[u] > 1 or sizes[v] > 1))
            literal = state == kmove(expanded, expanded_word(target, blocks))
            need(parity == literal, "oriented terminal-port parity criterion")
            checks += 1
            if r == 3 and first_size == last_size == 2:
                need(failures[:2] == [None, None] and failures[2] is not None and not literal,
                     "mirrored endpoint-two hostile")
                hostile = dict(macro_order=n, sizes=sizes,
                               macro_K_words=ledger,
                               expanded_K_words=[expanded_word(w, blocks) for w in ledger],
                               failures=failures,
                               terminal_labelled_target_matches=literal)
    return dict(exact_lift_predicate_controls=checks, mirrored_naive_lift_hostile=hostile)


def last_port_presentations(n, sample=None):
    """Blue standard path, closed endpoint scores(n-3,2), c->x."""
    constraints = {(i, i + 1): 1 for i in range(n - 1)}
    constraints.update({(0, j): int(j < n - 2) for j in range(1, n)})
    constraints[2, n - 2] = 0
    unseen = set(combinations(range(n), 2)); base = set(); free = []
    while unseen:
        u, v = min(unseen); orbit = {(u, v), (n - 1 - v, n - 1 - u)}
        unseen -= orbit
        values = {constraints[p] for p in orbit if p in constraints}
        need(len(values) <= 1, "compatible reflected boundary constraints")
        if values:
            bit = next(iter(values)); base.update((a, b) if bit else (b, a) for a, b in orbit)
        else:
            free.append(orbit)
    assignments = product((0, 1), repeat=len(free))
    if sample is not None:
        rng = random.Random(n * 1000 + 601)
        assignments = [tuple(rng.randrange(2) for _ in free) for _ in range(sample)]
    for bits in assignments:
        arcs = set(base)
        for orbit, bit in zip(free, bits):
            arcs.update((a, b) if bit else (b, a) for a, b in orbit)
        yield frozenset(arcs)


def last_port_ledger(n, cut=None):
    # D starts at position3; U begins at cut and ends at c=n-2.
    if cut is None:
        cut = n - 2
    dblock, ublock = tuple(range(3, cut)), tuple(range(cut, n - 1))
    need(dblock and ublock, "nonempty D and U")
    return [(0, 2) + dblock + ublock + (n - 1, 1),
            dblock + (0, 2) + ublock + (n - 1, 1),
            dblock + (0, 1, 2) + ublock + (n - 1,)]


def check_last_port():
    records = []
    for n, sample in [(7, None), (8, None), (9, None), (10, 16), (11, 16), (16, 16), (25, 16)]:
        counts = Counter(); residual = []
        for arcs in last_port_presentations(n, sample):
            target = tuple(range(n)); d = scores(arcs, n)
            need(blue(arcs, target) and d[0] == n - 3 and d[-1] == 2, "last-port universe")
            words = last_port_ledger(n); state = arcs; flags = []
            for word in words:
                need(directed(state, word), "block ledger always directed")
                flags.append(blue(state, word)); state = kmove(state, word)
            need(state == kmove(arcs, target), "block ledger exact K endpoint")
            need(not flags[0], "first word black for n>=7")
            need(not flags[1] or n in (7, 9), "second word exceptional orders")
            need(not flags[2] or n in (8, 9), "third word exceptional orders")
            counts["presentations"] += 1
            counts["fixed_tail_word_all_black"] += not any(flags)
            if sample is not None:
                need(not any(flags), "uniform last-port theorem for n>=10")
            legal_cut_success = False
            if n <= 9:
                for cut in range(4, n - 1):
                    if (cut, 2) not in arcs:
                        continue
                    state = arcs; ok = True
                    for word in last_port_ledger(n, cut):
                        need(directed(state, word), "guarded variable-cut ledger directed")
                        ok = ok and not blue(state, word)
                        state = kmove(state, word)
                    need(state == kmove(arcs, target), "guarded variable-cut exact endpoint")
                    legal_cut_success |= ok
                counts["some_legal_cut_all_black"] += legal_cut_success
                if not legal_cut_success:
                    residual.append([sum(1 << v for a, v in arcs if a == u) for u in range(n)])
        records.append(dict(n=n, finite_universe="complete fixed-boundary presentations" if sample is None
                            else "sixteen deterministic mirrored-chord samples", counts=dict(counts),
                            variable_cut_failures=residual))
    return records


def main():
    result = dict(status="PASS", source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  family=check_family(), port_criterion=check_port_criterion(),
                  last_port=check_last_port())
    (HERE / "lift_check.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                                          encoding="utf-8", newline="\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
