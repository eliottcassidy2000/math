"""Folding a tournament with an anti-involution that fixes one vertex (Theorem Q of the note).

fold(out, n, psi) returns (Z, D1, plus): plus = the out-neighbours of the fixed vertex v0 (one vertex of each pair),
Z = the quotient tournament on the pairs (i -> j iff exactly one of a_i -> a_j, a_i -> psi(a_j) holds, a_i in plus),
D1 = the set of pairs with a_i -> psi(a_i).  paths_ending_in(Z, D1) = number of directed paths of Z ending in D1.
Theorem Q: H(T) = 1 + 2 * paths_ending_in(Z, D1)  (mod 4)."""
def fold(out, n, psi):
    fixed = [a for a in range(n) if psi[a] == a]
    assert len(fixed) == 1 and all(psi[psi[a]] == a for a in range(n))
    arc = lambda a, b: (out[a] >> b) & 1
    assert all(arc(a, b) == arc(psi[b], psi[a]) for a in range(n) for b in range(n) if a != b)
    v0 = fixed[0]
    plus = [a for a in range(n) if a != v0 and arc(v0, a)]
    m = len(plus); assert 2 * m + 1 == n
    Z = [0] * m; D1 = 0
    for i, a in enumerate(plus):
        if arc(a, psi[a]): D1 |= 1 << i
        for j, b in enumerate(plus):
            if i != j and (arc(a, b) ^ arc(a, psi[b])): Z[i] |= 1 << j
    assert all(((Z[i] >> j) & 1) != ((Z[j] >> i) & 1) for i in range(m) for j in range(i + 1, m))
    return Z, D1, plus
def paths_ending_in(Z, D1, by_length=False):
    m = len(Z); dp = [[0] * m for _ in range(1 << m)]
    for v in range(m): dp[1 << v][v] = 1
    tot = [0] * (m + 1)
    for S in range(1, 1 << m):
        k = bin(S).count('1'); row = dp[S]
        for v in range(m):
            c = row[v]
            if not c: continue
            if (D1 >> v) & 1: tot[k] += c
            t = Z[v] & ~S
            while t:
                w = (t & -t).bit_length() - 1; dp[S | (1 << w)][w] += c; t &= t - 1
    return tot if by_length else sum(tot)
def ham(out, n):
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            t = out[v] & ~S
            while t:
                w = (t & -t).bit_length() - 1; dp[S | (1 << w)][w] += c; t &= t - 1
    return sum(dp[full])
