"""Independent check of the two spectrum witnesses quoted in README §7 / THM-4862 (h).

Each witness is a strong tournament given by out-row bitmasks; vertex 0 is the sole-exit vertex,
its unique out-neighbour r is the root; the digit is e_r(T - 0) = number of Hamiltonian paths of
T - 0 ending at r (THM-4862 (b): equivalently the number of Hamiltonian cycles of the entry attachment).
Pure Python subset DP, no dependence on the C generators.
"""
from functools import lru_cache

WITNESSES = {
    "12877 in R2(12)": ([512, 17, 11, 19, 1797, 31, 63, 127, 239, 1518, 2543, 1023], 2),
    "38615 in R1(13)": ([2048, 257, 11, 259, 15, 31, 63, 127, 1781, 255, 2815, 5118, 2047], 1),
}

def is_tournament(rows):
    n = len(rows)
    for i in range(n):
        if rows[i] >> i & 1:
            return False
        for j in range(i + 1, n):
            if (rows[i] >> j & 1) == (rows[j] >> i & 1):
                return False
    return True

def strong(rows):
    n = len(rows)
    def reach(adj):
        seen = 1; stack = [0]
        while stack:
            v = stack.pop(); m = adj[v] & ~seen
            while m:
                b = m & -m; m ^= b; seen |= b; stack.append(b.bit_length() - 1)
        return seen == (1 << n) - 1
    rev = [sum(1 << j for j in range(n) if rows[j] >> i & 1) for i in range(n)]
    return reach(rows) and reach(rev)

def ham_paths_ending(rows, verts):
    """dict end -> number of Hamiltonian paths of the sub-tournament on `verts` ending at end."""
    verts = list(verts); idx = {v: k for k, v in enumerate(verts)}; m = len(verts)
    adj = [sum(1 << idx[u] for u in verts if rows[v] >> u & 1) for v in verts]
    full = (1 << m) - 1
    dp = [[0] * m for _ in range(1 << m)]
    for k in range(m):
        dp[1 << k][k] = 1
    for S in range(1, 1 << m):
        row = dp[S]
        for k in range(m):
            c = row[k]
            if not c:
                continue
            nxt = adj[k] & ~S
            while nxt:
                b = nxt & -nxt; nxt ^= b; j = b.bit_length() - 1
                dp[S | b][j] += c
    return {verts[k]: dp[full][k] for k in range(m)}

ok = True
for name, (rows, want_digit) in WITNESSES.items():
    n = len(rows)
    assert is_tournament(rows), name
    H = sum(ham_paths_ending(rows, range(n)).values())
    outs = [j for j in range(n) if rows[0] >> j & 1]
    assert len(outs) == 1, (name, outs)
    r = outs[0]
    digit = ham_paths_ending(rows, range(1, n))[r]
    st = strong(rows)
    target = int(name.split()[0])
    line = f"{name}: order {n}, H = {H}, strong = {st}, sole exit 0 -> {r}, digit e_r(T - 0) = {digit}"
    good = (H == target and st and digit == want_digit)
    ok &= good
    print(line, "OK" if good else "FAIL")
print("ALL OK" if ok else "SOME FAIL")
