"""aud_block_construction.py -- audit check (fermat_self_duality_20261104_audit_S10) of the construction added to the
note's section 5 (a) while the audit ran: D = AT_m + u (u beaten by m - 1 and m, beating 1..m-2), plus l with u -> l the
only arc leaving D and l beating all of AT_m.  Claimed: order m + 2, no sole exit, hc = 2^(m-3) + 1 (so F_i Hamiltonian
cycles at order 2^i + 5).  Own code: subset DP for Hamiltonian cycles."""

def hc(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    dp[1][0] = 1
    for S in range(1, full + 1, 2):
        for v in range(n):
            c = dp[S][v]
            if c:
                nb = out[v] & ~S
                while nb:
                    w = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    dp[S | 1 << w][w] += c
    return sum(dp[full][v] for v in range(1, n) if out[v] & 1)

ok = True
for m in range(4, 13):
    n = m + 2
    out = [0] * n                      # vertices 0..m-1 = AT_m (labels 1..m), m = u, m + 1 = l
    for i in range(m):
        for j in range(i + 1, m):
            if (i, j) == (0, m - 1):
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
    u, l = m, m + 1
    for v in range(m):
        if v in (m - 2, m - 1):
            out[v] |= 1 << u
        else:
            out[u] |= 1 << v
        out[l] |= 1 << v
    out[u] |= 1 << l
    sole_exits = [v for v in range(n) if bin(out[v]).count("1") == 1]
    h = hc(n, out)
    good = h == 2 ** (m - 3) + 1 and not sole_exits
    ok &= good
    print(f"m = {m:2d}: order {n:2d}, sole exits {sole_exits}, hc = {h} (2^(m-3) + 1 = {2 ** (m - 3) + 1}) {'ok' if good else 'FAIL'}")
print("ALL OK (F_i Hamiltonian cycles at order 2^i + 5 for i = 0..3, m = 4, 5, 7, 11)" if ok else "FAIL")
