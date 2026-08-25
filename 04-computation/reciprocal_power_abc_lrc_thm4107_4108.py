from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    out = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def radical(n):
    ans = 1
    for p in factor(abs(n)):
        ans *= p
    return ans


def primitive_terms(a, b):
    aa = a**b
    bb = b**a
    common = gcd(aa, bb)
    return aa // common, bb // common, common


def primitive_power_edge(a, b):
    """True iff a points to b after pairwise gcd reduction."""
    g = gcd(a, b)
    x, y = a // g, b // g
    return x**y > y**x


def exceptional_sorted_pair(a, b):
    require(a < b, "sorted pair required")
    return b % a == 0 or 3 * a == 2 * b


def lrc_observer(speeds, t):
    def circle_distance(x):
        r = x - (x.numerator // x.denominator)
        return min(r, 1 - r)

    return min(circle_distance(v * t) for v in speeds)


def triangle_is_cyclic(a, b, c):
    verts = (a, b, c)
    outdegree = []
    for x in verts:
        outdegree.append(
            sum(x != y and primitive_power_edge(x, y) for y in verts)
        )
    return sorted(outdegree) == [1, 1, 1]


def labelled_tournament_signature(speeds):
    return tuple(
        primitive_power_edge(speeds[i], speeds[j])
        for i, j in combinations(range(len(speeds)), 2)
    )


def triangle_formula(a, b, c):
    require(a < b < c, "sorted triple required")
    e_ab = exceptional_sorted_pair(a, b)
    e_bc = exceptional_sorted_pair(b, c)
    e_ac = exceptional_sorted_pair(a, c)
    return (e_ac and not e_ab and not e_bc) or (
        not e_ac and e_ab and e_bc
    )


def tournament_stats(speeds):
    n = len(speeds)
    edges = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        edges[i][j] = primitive_power_edge(speeds[i], speeds[j])
        edges[j][i] = not edges[i][j]

    c3 = sum(
        triangle_is_cyclic(speeds[i], speeds[j], speeds[k])
        for i, j, k in combinations(range(n), 3)
    )

    @lru_cache(None)
    def paths(mask, last):
        if mask == (1 << last):
            return 1
        previous = mask ^ (1 << last)
        return sum(
            paths(previous, v)
            for v in range(n)
            if previous & (1 << v) and edges[v][last]
        )

    full = (1 << n) - 1
    hamilton_paths = sum(paths(full, last) for last in range(n))
    reversed_edges = [
        (speeds[i], speeds[j])
        for i, j in combinations(range(n), 2)
        if not edges[i][j]
    ]
    # Kosaraju SCC sizes, sorted in condensation order is unnecessary here.
    seen = set()
    finish = []

    def dfs(v):
        seen.add(v)
        for w in range(n):
            if edges[v][w] and w not in seen:
                dfs(w)
        finish.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)
    seen.clear()
    scc_sizes = []

    def reverse_dfs(v):
        seen.add(v)
        size = 1
        for w in range(n):
            if edges[w][v] and w not in seen:
                size += reverse_dfs(w)
        return size

    for v in reversed(finish):
        if v not in seen:
            scc_sizes.append(reverse_dfs(v))
    return reversed_edges, c3, hamilton_paths, tuple(scc_sizes)


# Exact pair classification and triangle formula.
pair_count = 0
for a in range(1, 121):
    for b in range(a + 1, 121):
        pair_count += 1
        require(
            primitive_power_edge(a, b) == (not exceptional_sorted_pair(a, b)),
            f"pair classification failed at {(a, b)}",
        )

triple_count = 0
cycles = []
for a, b, c in combinations(range(1, 61), 3):
    triple_count += 1
    observed = triangle_is_cyclic(a, b, c)
    require(observed == triangle_formula(a, b, c), (a, b, c))
    if observed:
        cycles.append((a, b, c))

cycles_by_max = sorted(cycles, key=lambda t: (t[2], t))
require(cycles_by_max[:2] == [(2, 5, 6), (3, 5, 6)], cycles_by_max[:2])

# Valuation-defect positive/negative parts and primitive additive packets.
packet_count = 0
for a in range(1, 61):
    for b in range(1, 61):
        if a == b:
            continue
        packet_count += 1
        A, B, common = primitive_terms(a, b)
        require(gcd(A, B) == 1, (a, b, A, B))
        require(common == gcd(a**b, b**a), (a, b))
        support = set(factor(a)) | set(factor(b))
        reconstructed_A = 1
        reconstructed_B = 1
        defect_support = 1
        for p in support:
            va = factor(a).get(p, 0)
            vb = factor(b).get(p, 0)
            defect = b * va - a * vb
            if defect > 0:
                reconstructed_A *= p**defect
            elif defect < 0:
                reconstructed_B *= p ** (-defect)
            if defect:
                defect_support *= p
        require((A, B) == (reconstructed_A, reconstructed_B), (a, b))
        require(radical(A * B) == defect_support, (a, b))
        require(gcd(A, A + B) == gcd(B, A + B) == 1, (a, b))

# A common-dilation ray compiles the normalized packet to m-th powers.
ray_count = 0
for r in range(1, 16):
    for s in range(r + 1, 16):
        if gcd(r, s) != 1:
            continue
        for m in range(1, 10):
            ray_count += 1
            a, b = m * r, m * s
            U, V = (m * r) ** s, (m * s) ** r
            d = gcd(U, V)
            X, Y = U // d, V // d
            A, B, _ = primitive_terms(a, b)
            require((A, B) == (X**m, Y**m), (r, s, m))
            raw_edge = a**b > b**a
            threshold_edge = m ** (s - r) * r**s > s**r
            require(raw_edge == threshold_edge, (r, s, m))

# Coprime sum/difference channels are primitive, and their gcd is parity-only.
coprime_channel_count = 0
for a in range(2, 50):
    for b in range(a + 1, 50):
        if gcd(a, b) != 1:
            continue
        coprime_channel_count += 1
        A, B = a**b, b**a
        require(A != B, (a, b))
        plus, minus = A + B, abs(A - B)
        height = max(A, B)
        require(9 * minus >= height, (a, b, minus, height))
        if 9 * minus == height:
            require((a, b) == (2, 3), (a, b))
        require(gcd(A, B) == gcd(A, plus) == gcd(B, plus) == 1, (a, b))
        require(gcd(min(A, B), minus) == gcd(max(A, B), minus) == 1, (a, b))
        require(gcd(plus, minus) == (2 if a % 2 and b % 2 else 1), (a, b))

# Every primitive positive ABC packet induces a transitive primitive-power K3.
abc_count = 0
for A in range(1, 151):
    for B in range(A + 1, 151):
        C = A + B
        if gcd(A, B) != 1:
            continue
        abc_count += 1
        require(not triangle_is_cyclic(A, B, C), (A, B, C))

# LRC straddle ABC normalization, its lift gauge, and an unbounded q/D hostile.
straddle_normalizations = 0
for u in range(1, 20):
    for w in range(1, 20):
        for alpha in range(1, 10):
            for beta in range(1, 10):
                X, Y = u * beta, w * alpha
                if X == Y:
                    continue
                D = abs(X - Y)
                G = gcd(X, Y)
                xx, yy, dd = X // G, Y // G, D // G
                require(gcd(xx, yy) == gcd(xx, dd) == gcd(yy, dd) == 1, (u, w, alpha, beta))
                require(max(xx, yy) == min(xx, yy) + dd, (u, w, alpha, beta))
                straddle_normalizations += 1

for w in range(1, 101):
    u, alpha, beta = w + 1, 1, 1
    q, D = u + w, u * beta - w * alpha
    p = alpha + beta
    require(D == 1 and q == 2 * w + 1, w)
    require(Fraction(u * p, q) - alpha == Fraction(D, q), w)
    require(Fraction(w * p, q) - beta == -Fraction(D, q), w)

for lift in range(0, 101):
    u, w, alpha, beta = 1, 13, 0 + lift, 1 + 13 * lift
    X, Y = u * beta, w * alpha
    require(X - Y == 1, lift)
    require((X, Y) == (1 + 13 * lift, 13 * lift), lift)

# Exact LRC blindness hostiles.
require(primitive_power_edge(1, 2) is False, "{1,2} orientation")
require(primitive_power_edge(2, 3) is False, "{2,3} orientation")
require(lrc_observer((1, 2), Fraction(1, 3)) == Fraction(1, 3), "pair AP")
require(lrc_observer((2, 3), Fraction(1, 5)) == Fraction(2, 5), "pair shifted")

odd_primes = (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
mixed_primes = (2, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43)
for row in (odd_primes, mixed_primes):
    require(
        all(primitive_power_edge(a, b) for a, b in combinations(row, 2)),
        row,
    )
require(lrc_observer(odd_primes, Fraction(1, 2)) == Fraction(1, 2), "odd row")
require(lrc_observer(mixed_primes, Fraction(1, 2)) == 0, "mixed half phase")
require(lrc_observer(mixed_primes, Fraction(1, 4)) == Fraction(1, 4), "mixed quarter")

ap13 = tuple(range(1, 14))
ap_reversals, ap_c3, ap_h, ap_scc = tournament_stats(ap13)
selected_rows = {
    "AP13": ap13,
    "AP12+182": tuple(range(1, 13)) + (182,),
    "AP12+339": tuple(range(1, 13)) + (339,),
    "AP12+5460": tuple(range(1, 13)) + (5460,),
    "26AP12+339": tuple(26 * k for k in range(1, 13)) + (339,),
}
selected_stats = {
    name: tournament_stats(row)[1:] for name, row in selected_rows.items()
}
period27_signature_distance = sum(
    x != y
    for x, y in zip(
        labelled_tournament_signature(selected_rows["AP13"]),
        labelled_tournament_signature(selected_rows["26AP12+339"]),
    )
)
period27_edge_differences = []
for index_pair, x, y in zip(
    combinations(range(13), 2),
    labelled_tournament_signature(selected_rows["AP13"]),
    labelled_tournament_signature(selected_rows["26AP12+339"]),
):
    if x != y:
        i, j = index_pair
        period27_edge_differences.append(
            (
                (selected_rows["AP13"][i], selected_rows["AP13"][j]),
                (selected_rows["26AP12+339"][i], selected_rows["26AP12+339"][j]),
            )
        )

ap_to_period27 = {0: 12, 12: 0, **{i: i for i in range(1, 12)}}
for i, j in combinations(range(13), 2):
    ap_edge = primitive_power_edge(
        selected_rows["AP13"][i], selected_rows["AP13"][j]
    )
    image_i, image_j = ap_to_period27[i], ap_to_period27[j]
    image_edge = primitive_power_edge(
        selected_rows["26AP12+339"][image_i],
        selected_rows["26AP12+339"][image_j],
    )
    require(ap_edge == image_edge, (i, j, image_i, image_j))

print("ABC POWER RECIPROCITY EXACT AUDIT")
print(f"pair_classification={pair_count}")
print(f"triangle_classification={triple_count}")
print(f"cycles={len(cycles)} first_by_max={cycles_by_max[:10]}")
print(f"valuation_packets={packet_count}")
print(f"common_dilation_rays={ray_count}")
print(f"coprime_plus_minus_channels={coprime_channel_count}")
print(f"primitive_abc_packets={abc_count}")
print(f"normalized_straddle_packets={straddle_normalizations}")
print("straddle_hostile=(u,w,alpha,beta)=(n+1,n,1,1), q/D=2n+1 unbounded")
print("AP13_lift_gauge=(13N)+1=(13N)+1 leaves D=1 and the physical time unchanged")
print(f"ap13_reversed_edges={ap_reversals}")
print(f"ap13_c3={ap_c3} ap13_H={ap_h} ap13_SCC={ap_scc}")
print(f"selected_lrc_rows={selected_stats}")
print(f"AP13_period27_labelled_edge_distance={period27_signature_distance}")
print(f"AP13_period27_edge_differences={period27_edge_differences}")
print("AP13_period27_isomorphism=swap endpoint positions 1 and 13")
print("pair_hostile={1,2}:1/3 vs {2,3}:2/5")
print("thirteen_hostile=identical transitive tournament; odd M=1/2, mixed 1/4<=M<1/2")
print("PASS")


