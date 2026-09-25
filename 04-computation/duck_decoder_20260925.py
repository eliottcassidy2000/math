"""Exact controls for the four-color / ten-letter decoder, standard library only.

Universe: all ten unordered pairs of four XOR colors, words of lengths1..4;
all tournaments of even orders2,4,6; all sixteen signed2x2 blocks.
Positive controls: explicit orbits, direct HP enumeration, reversible blocks.
Hostile controls: first singular HP kernel; loss of block contrasts;
no C3-invariant perfect matching. No asymptotic Collatz claim is tested.
"""
from collections import Counter
from itertools import combinations, combinations_with_replacement, permutations, product


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


COLORS = range(4)
ROT = (0, 2, 3, 1)
PAIRS = tuple(combinations_with_replacement(COLORS, 2))
S = frozenset(p for p in PAIRS if p[0] != p[1] or p == (0, 0))
U = frozenset(PAIRS) - S


def rotate_pair(p):
    return tuple(sorted(ROT[x] for x in p))


def rotate_word(w):
    return tuple(rotate_pair(p) for p in w)


def charge(w):
    answer = 0
    for a, b in w:
        answer ^= a ^ b
    return answer


def k5_edge(p):
    a, b = p
    if a == b:
        return (a, 4)
    return (a, b)


def color_words():
    require(len(S) == 7 and len(U) == 3, "7+3 split")
    require(Counter(a ^ b for a, b in PAIRS) == {0: 4, 1: 2, 2: 2, 3: 2}, "charge multiplicities")
    require(all(ROT[a ^ b] == ROT[a] ^ ROT[b] for a in COLORS for b in COLORS), "linear color action")
    require(len({k5_edge(p) for p in PAIRS}) == 10, "K5 bijection")
    k5_rot = (0, 2, 3, 1, 4)
    require(all(k5_edge(rotate_pair(p)) == tuple(sorted(k5_rot[x] for x in k5_edge(p))) for p in PAIRS), "K5 equivariance")
    print("C3 action: Sym^2_set(F2^2) = 1 fixed letter + 3 free three-cycles; S=7,U=3")
    print("K5 edge action: equivariant bijection, two fixed vertices")
    for k in range(1, 5):
        words = tuple(product(PAIRS, repeat=k))
        deleted = {tuple([p] * k) for p in S}
        remaining = set(words) - deleted
        representatives = set()
        full_representatives = set()
        neutral_representatives = set()
        for w in words:
            rw = rotate_word(w)
            rrw = rotate_word(rw)
            require(rotate_word(rrw) == w, "action order")
            representative = min(w, rw, rrw)
            full_representatives.add(representative)
            if w in remaining:
                require(len({w, rw, rrw}) == 3, "free residual action")
                representatives.add(representative)
                if charge(w) == 0:
                    neutral_representatives.add(representative)
        counts = Counter(map(charge, words))
        require(counts[0] == (10**k + 3 * 2**k) // 4, "neutral Fourier formula")
        require(all(counts[c] == (10**k - 2**k) // 4 for c in (1, 2, 3)), "nonzero Fourier formula")
        require(len(full_representatives) == (10**k + 2) // 3, "Burnside full")
        total = (10**k - 7) // 3
        neutral = (10**k + 3 * 2**k - (4 if k % 2 else 28)) // 12
        nonzero = (10**k - 2**k) // 4 - (2 if k % 2 else 0)
        require((len(representatives), len(neutral_representatives)) == (total, neutral), "residual split")
        require(total == neutral + nonzero, "split conservation")
        print(f"k={k}: residual color-word orbits={total}, neutral={neutral}, nonzero={nonzero}")
    for k in range(1, 10):
        total = (10**k - 7) // 3
        require(total == 10 * ((10**(k-1) - 7) // 3) + 21 if k >= 2 else total == 1, "decimal recursion")
    print("k=9: orbit count333333331 = 17 * 19607843 (counting does not imply primality)")
    require((10**9 - 7) // 3 == 17 * 19607843, "factorization")


def adjacency(n, code):
    edges = tuple(combinations(range(n), 2))
    adj = [0] * n
    for i, (a, b) in enumerate(edges):
        if code >> i & 1:
            adj[a] |= 1 << b
        else:
            adj[b] |= 1 << a
    return edges, adj


def matching_parity(n, code):
    edges, adj = adjacency(n, code)
    edge_bit = {e: 1 << i for i, e in enumerate(edges)}
    dp = [{} for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = (1, 0)
    for mask in range(1, 1 << n):
        is_pair = mask.bit_count() % 2
        for v, (count, kernel) in dp[mask].items():
            todo = adj[v] & ~mask
            while todo:
                bit = todo & -todo
                todo ^= bit
                w = bit.bit_length() - 1
                old_count, old_kernel = dp[mask | bit].get(w, (0, 0))
                added = edge_bit[tuple(sorted((v, w)))] if is_pair and count else 0
                dp[mask | bit][w] = (old_count ^ count, old_kernel ^ kernel ^ added)
    count = kernel = 0
    for c, k in dp[-1].values():
        count ^= c
        kernel ^= k
    rows = [0] * n
    for i, (a, b) in enumerate(edges):
        if kernel >> i & 1:
            rows[a] |= 1 << b
            rows[b] |= 1 << a
    require(count == 1, "odd Hamiltonian-path control")
    require(all(r.bit_count() % 2 for r in rows), "odd-degree parity kernel")
    basis = {}
    for row in rows:
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    return kernel, rows, len(basis)


def direct_matching_parity(n, code):
    edges, adj = adjacency(n, code)
    edge_bit = {e: 1 << i for i, e in enumerate(edges)}
    answer = count = 0
    for path in permutations(range(n)):
        if all(adj[a] >> b & 1 for a, b in zip(path, path[1:])):
            count += 1
            for i in range(0, n, 2):
                answer ^= edge_bit[tuple(sorted(path[i:i+2]))]
    return count, answer


def matchings(vertices):
    if not vertices:
        yield ()
        return
    first, *rest = vertices
    for i, second in enumerate(rest):
        for tail in matchings(rest[:i] + rest[i+1:]):
            yield tuple(sorted(((first, second),) + tail))


def tournament_controls():
    six_kernels = []
    rank_two_degrees = Counter()
    first_rank_two = {}
    for n in (2, 4, 6):
        ranks = Counter()
        first_singular = None
        for code in range(1 << (n * (n-1) // 2)):
            kernel, rows, rank = matching_parity(n, code)
            ranks[rank] += 1
            if n == 6:
                six_kernels.append(kernel)
                if rank == 2:
                    degrees = tuple(sorted(row.bit_count() for row in rows))
                    rank_two_degrees[degrees] += 1
                    first_rank_two.setdefault(degrees, code)
            if n <= 4:
                count, direct = direct_matching_parity(n, code)
                require(count % 2 == 1 and direct == kernel, "independent permutation path")
            if rank < n and first_singular is None:
                first_singular = code
        print(f"HP matching parity: n={n}, ranks={dict(sorted(ranks.items()))}, first singular={first_singular}")
    require(rank_two_degrees == {(1,1,1,1,1,5): 960, (3,3,3,3,3,3): 720}, "rank-two support census")
    print(f"rank2 supports: K1,5=960 (first code72), K3,3=720 (first code83)")
    require(first_rank_two == {(1,1,1,1,1,5): 72, (3,3,3,3,3,3): 83}, "support first witnesses")
    # Independent route: complete each Hamiltonian order to all tournaments
    # containing it, rather than enumerate paths inside a fixed tournament.
    edges = tuple(combinations(range(6), 2))
    edge_bit = {e: 1 << i for i, e in enumerate(edges)}
    completed_counts, completed_kernels = [0] * 32768, [0] * 32768
    for path in permutations(range(6)):
        fixed_mask = fixed_value = paired_edges = 0
        for a, b in zip(path, path[1:]):
            bit = edge_bit[tuple(sorted((a, b)))]
            fixed_mask |= bit
            if a < b:
                fixed_value |= bit
        for i in range(0, 6, 2):
            paired_edges |= edge_bit[tuple(sorted(path[i:i+2]))]
        completions = [fixed_value]
        for i in range(15):
            bit = 1 << i
            if not fixed_mask & bit:
                completions += [code | bit for code in completions]
        for code in completions:
            completed_counts[code] += 1
            completed_kernels[code] ^= paired_edges
    require(completed_kernels == six_kernels, "all-kernel independent path-completion census")
    require(all(c % 2 for c in completed_counts), "all independent path counts odd")
    print("independent path-completion audit: 720 orders x1024 completions, all32768 kernels agree")
    kernel, rows, rank = matching_parity(6, 8)
    count, direct = direct_matching_parity(6, 8)
    require(count == 9 and direct == kernel and rank == 4, "singular hostile")
    for null in ((1 << 0) | (1 << 4), (1 << 2) | (1 << 5)):
        require(all((row & null).bit_count() % 2 == 0 for row in rows), "explicit null vector")
    print("n6 code8: nine HPs; null vectors e0+e4 and e2+e5; universal symplectic halving refuted")
    for n in (4, 10):
        rotate = tuple([1, 2, 0] + list(range(3, n)))
        fixed = 0
        total = 0
        for matching in matchings(list(range(n))):
            total += 1
            rotated = tuple(sorted(tuple(sorted((rotate[a], rotate[b]))) for a, b in matching))
            fixed += rotated == matching
        require(fixed == 0, "no equivariant matching")
        print(f"independent C3 rotation: n={n}, matchings={total}, invariant matchings={fixed}")


def block_controls():
    for a, b, c, d in product((-1, 1), repeat=4):
        s, r, q, h = a+b+c+d, a+b-c-d, a-b+c-d, a-b-c+d
        rebuilt = ((s+r+q+h)//4, (s+r-q-h)//4, (s-r+q-h)//4, (s-r-q+h)//4)
        require(rebuilt == (a, b, c, d), "four-channel inverse")
        require((r == q == h == 0) == (a == b == c == d), "module contrast boundary")
    # Two different balanced blocks have equal mean and equal flip cost2.
    b1, b2 = (1, 1, -1, -1), (1, -1, -1, 1)
    require(sum(b1) == sum(b2) == 0 and b1 != b2, "mean-only hostile")
    print("all16 signed2x2 blocks reconstruct exactly; module boundary = zero three contrasts")


def color_clock_controls():
    # Least-significant color digit first; the action needs carries.
    def step(digits):
        digits = list(digits)
        for i in range(len(digits)):
            digits[i] = (digits[i] + 1) % 3
            if digits[i]:
                break
        return tuple(digits)

    for depth in range(1, 6):
        modulus = 3**depth
        digits = (0,) * depth
        repunit = height = 0
        seen_digits, seen_repunits, seen_heights = set(), set(), set()
        for j in range(modulus):
            require(sum(x * 3**i for i, x in enumerate(digits)) == j, "color digit decoder")
            seen_digits.add(digits)
            seen_repunits.add(repunit)
            seen_heights.add(height)
            digits = step(digits)
            repunit = (10 * repunit + 1) % modulus
            height = (4 * height + 1) % modulus
        require(len(seen_digits) == len(seen_repunits) == len(seen_heights) == modulus, "full residue cycles")
        require(digits == (0,) * depth and repunit == height == 0, "clock return")
    w = (0, 0)
    for _ in range(3):
        w = step(w)
    require(w == (0, 1), "depth2 carry hostile")
    print("color odometer and R4/R10 clocks: full cycles at depths1..5; diagonal color rotation fails at depth2")


if __name__ == "__main__":
    color_words()
    tournament_controls()
    block_controls()
    color_clock_controls()
    print("PASS: exact finite controls; no convergence or graceful/square-sum equivalence claim")
