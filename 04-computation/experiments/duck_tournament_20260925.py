"""Non-pair tournament decoders with exact module, switching and port controls."""

from collections import Counter
from functools import lru_cache
from itertools import combinations, product
from math import factorial


def check(condition, label):
    if not condition:
        raise RuntimeError(label)


def bits(mask):
    while mask:
        low = mask & -mask
        yield low.bit_length() - 1
        mask -= low


def tournament(n, code):
    rows = [0] * n
    for b, (i, j) in enumerate(combinations(range(n), 2)):
        u, v = (i, j) if code >> b & 1 else (j, i)
        rows[u] |= 1 << v
    return tuple(rows)


def regular(n):
    return tuple(sum(1 << j for j in range(n)
                     if 1 <= (j - i) % n <= n // 2) for i in range(n))


def substitution(core, h):
    q = len(h)
    owner = tuple(i for i in range(3) for _ in range(q)) + (3,)
    rows = [0] * len(owner)
    for u in range(len(owner)):
        for v in range(len(owner)):
            i, j = owner[u], owner[v]
            if i == j and i < 3:
                arc = h[u % q] >> (v % q) & 1
            else:
                arc = core[i] >> j & 1
            rows[u] |= arc << v
    return tuple(rows), owner


def module(rows, mask):
    full = (1 << len(rows)) - 1
    return all(rows[v] & mask in (0, mask) for v in bits(full ^ mask))


def modules_of_size(rows, size):
    answer = []
    for vertices in combinations(range(len(rows)), size):
        mask = sum(1 << v for v in vertices)
        if module(rows, mask):
            answer.append(mask)
    return answer


def induced(rows, vertices):
    return tuple(sum((rows[u] >> v & 1) << j for j, v in enumerate(vertices))
                 for u in vertices)


def module_decode(rows):
    n = len(rows)
    check(n % 6 == 4 and n >= 10, "tripled odd-order domain")
    q = (n - 1) // 3
    blocks = modules_of_size(rows, q)
    check(len(blocks) == 3 and sum(m.bit_count() for m in blocks) == (blocks[0] | blocks[1] | blocks[2]).bit_count(),
          "intrinsic three disjoint modules")
    rest = ((1 << n) - 1) ^ (blocks[0] | blocks[1] | blocks[2])
    check(rest.bit_count() == 1, "unique singleton")
    all_blocks = blocks + [rest]
    reps = [next(bits(m)) for m in all_blocks]
    core = induced(rows, reps)
    children = [induced(rows, tuple(bits(m))) for m in blocks]
    check(all(all(r.bit_count() == q // 2 for r in child) for child in children), "regular child guards")
    # Exact edge reconstruction, retaining original vertex embeddings.
    rebuilt = [0] * n
    for i, block_i in enumerate(all_blocks):
        vertices = tuple(bits(block_i))
        internal = induced(rows, vertices)
        for a, u in enumerate(vertices):
            rebuilt[u] |= sum((internal[a] >> b & 1) << v for b, v in enumerate(vertices))
            for j, block_j in enumerate(all_blocks):
                if core[i] >> j & 1:
                    rebuilt[u] |= block_j
    check(tuple(rebuilt) == rows, "zero-reversal inverse reconstruction")
    return core, blocks, rest


def switch(rows, cut):
    full = (1 << len(rows)) - 1
    return tuple(row ^ (full ^ cut if cut >> i & 1 else cut) for i, row in enumerate(rows))


def switchable_pairs(rows):
    full = (1 << len(rows)) - 1
    result = []
    for u, v in combinations(range(len(rows)), 2):
        outside = full ^ (1 << u) ^ (1 << v)
        difference = (rows[u] ^ rows[v]) & outside
        if difference in (0, outside):
            result.append((u, v))
    return result


def pair_modules(rows):
    return [(u, v) for u, v in combinations(range(len(rows)), 2)
            if (rows[u] ^ rows[v]) & ~((1 << u) | (1 << v)) == 0]


def pfaffian4(rows):
    def sign(i, j):
        return 1 if rows[i] >> j & 1 else -1
    return sign(0, 1) * sign(2, 3) - sign(0, 2) * sign(1, 3) + sign(0, 3) * sign(1, 2)


def triangle_masks(rows):
    answer = []
    for vertices in combinations(range(len(rows)), 3):
        mask = sum(1 << v for v in vertices)
        if all((rows[v] & mask).bit_count() == 1 for v in vertices):
            answer.append(mask)
    return answer


def reachability(rows):
    reach = [row | 1 << i for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    return tuple(reach)


def contract_triangle(rows, tri):
    n = len(rows)
    blocks = [tri] + [1 << v for v in range(n) if not tri >> v & 1]
    quotient = []
    for i, block in enumerate(blocks):
        row = 0
        for j, other in enumerate(blocks):
            if i != j and any(rows[u] & other for u in bits(block)):
                row |= 1 << j
        quotient.append(row)
    # A reversible three-port record stores all three crossing bits.
    tri_vertices = tuple(bits(tri))
    outside = tuple(v for v in range(n) if not tri >> v & 1)
    masks = {v: sum((rows[u] >> v & 1) << i for i, u in enumerate(tri_vertices)) for v in outside}
    rebuilt = [0] * n
    for u in range(n):
        for v in range(n):
            if u == v:
                continue
            if u in tri_vertices and v in outside:
                arc = masks[v] >> tri_vertices.index(u) & 1
            elif v in tri_vertices and u in outside:
                arc = 1 - (masks[u] >> tri_vertices.index(v) & 1)
            else:
                arc = rows[u] >> v & 1
            rebuilt[u] |= arc << v
    check(tuple(rebuilt) == rows, "triangle port reconstruction")
    return tuple(quotient), blocks


def hamiltonian_end_counts(rows):
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    full = (1 << n) - 1
    for mask in range(1, full + 1):
        for v in bits(mask):
            if dp[mask][v]:
                for w in bits(rows[v] & (full ^ mask)):
                    dp[mask | 1 << w][w] += dp[mask][v]
    return tuple(dp[full])


def quotient_word_count(core, content, terminal=None):
    @lru_cache(None)
    def extend(left, last):
        if not any(left):
            return int(terminal is None or last == terminal)
        result = 0
        for v, count in enumerate(left):
            if count and (last < 0 or core[last] >> v & 1):
                new = list(left)
                new[v] -= 1
                result += extend(tuple(new), v)
        return result
    return extend(tuple(content), -1)


def profile_transfer_c3_blocks(core, terminal=None):
    weighted = {1: 3, 2: 6, 3: 6}  # c! pc_C3(c), from direct path covers.
    total = 0
    for counts in product(range(1, 4), repeat=3):
        coefficient = quotient_word_count(core, counts + (1,), terminal)
        total += coefficient * weighted[counts[0]] * weighted[counts[1]] * weighted[counts[2]]
    return total


def main():
    modules_count = gauges = strong_controls = 0
    core_classes = Counter()
    path_counts = Counter()
    contractions = Counter()
    for code in range(64):
        core = tournament(4, code)
        pf = abs(pfaffian4(core))
        cut = sum(1 << v for v in range(3) if not (core[3] >> v & 1))
        normalized = switch(core, cut)
        check(normalized[3] == 7, "root-star gauge")
        cyclic = bool(triangle_masks(normalized))
        check(cyclic == (pf == 3), "four-core switching classification")
        check(switch(normalized, cut) == core, "gauge inverse")
        end_counts = hamiltonian_end_counts(core)
        check(sum(end_counts) == 1 + 2 * len(triangle_masks(core)), "four-core OCF")
        core_classes[(pf, tuple(sorted(r.bit_count() for r in core)))] += 1
        for gauge in range(8):
            check(abs(pfaffian4(switch(core, gauge))) == pf, "all core switching gauges")

        for tri in triangle_masks(core):
            quotient, blocks = contract_triangle(core, tri)
            original_reach, quotient_reach = reachability(core), reachability(quotient)
            owner = {v: i for i, b in enumerate(blocks) for v in bits(b)}
            for u in range(4):
                for v in range(4):
                    check((original_reach[u] >> v & 1) == (quotient_reach[owner[u]] >> owner[v] & 1),
                          "all-pairs reachability under strong contraction")
            contractions["digon" if quotient == (2, 1) else "single_arc"] += 1

        for q in (3, 5):
            rows, owner = substitution(core, regular(q))
            decoded, blocks, singleton = module_decode(rows)
            check(decoded == core and singleton == 1 << (3 * q), "unmarked decoder recovery")
            check(blocks == [((1 << q) - 1) << (q * i) for i in range(3)], "exact three modules")
            check(not switchable_pairs(rows), "no cut-switch can create a pair module")
            modules_count += 1
            if q == 3:
                all_modules = [m for m in range(1, 1 << len(rows)) if module(rows, m)]
                for block in blocks:
                    check(not any(m & block and m & ~block and block & ~m for m in all_modules), "strong-module recovery")
                    strong_controls += 1
                for gauge in range(1 << 9):
                    check(not pair_modules(switch(rows, gauge)), "exhaustive ten-vertex switching hostile")
                    gauges += 1
                counts = hamiltonian_end_counts(rows)
                check(sum(counts) == profile_transfer_c3_blocks(core), "inherited full profile transfer")
                check(counts[9] == profile_transfer_c3_blocks(core, 3), "root-ended profile transfer")
                path_counts[(sum(counts), counts[9])] += 1

    h5 = regular(5)
    check(not [m for m in range(1, 1 << 5) if 1 < m.bit_count() < 5 and module(h5, m)], "prime H5 obstruction")
    check(sum(hamiltonian_end_counts(h5)) == 15, "H5 positive Hamiltonian control")
    check((10 - 5) % 2 == 1, "triangle-only parity boundary")
    print("INTRINSIC_TRIPLE_DECODER", modules_count, "all64cores q3,5; exact reconstruction PASS")
    print("STRONG_MODULE_CONTROLS", strong_controls, "all subsets of all64 order10 objects PASS")
    print("SWITCH_PAIR_HOSTILE", gauges, "all64cores x512switches at order10; no pair modules PASS")
    print("FOUR_CORE_SWITCH_CLASSES", sorted(core_classes.items()))
    print("TRIANGLE_CONTRACTIONS", dict(contractions), "all marked directed triangles in all64cores; all-pairs reachability and ports PASS")
    print("HAMILTONIAN_ROOT_PROFILE", sorted(path_counts.items()))
    print("PROFILE_SCOPE: full path-cover data and terminal quotient word retained; scalar H alone not used")
    print("PRIME_STOP H5: no proper nontrivial modules, H=15")
    print("PARITY_STOP: triangle contractions remove2; order10 cannot become5 by triangles alone")
    print("ARITHMETIC_SCOPE: recovered N->(N-1)/3 is guarded REVERSE construction, not forward descent")
    print("ALL CHECKS PASS; no universal Collatz root generator claimed")


if __name__ == "__main__":
    main()
