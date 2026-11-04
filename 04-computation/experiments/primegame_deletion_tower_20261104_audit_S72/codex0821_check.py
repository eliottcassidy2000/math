# H and hc for the n=7 bitstrings quoted in 07-reflections/tournament-subset-decks-...-codex-20260821.md (3.7), (3.12)
# McKay format: upper triangle in row order, bit 1 means i -> j (i < j).
def adj_from(bits, n=7):
    adj = [0]*n; k = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits[k] == '1': adj[i] |= 1 << j
            else: adj[j] |= 1 << i
            k += 1
    return adj
def H_hc(adj, n=7):
    full = (1 << n) - 1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, full + 1):
        for v in range(n):
            c = dp[S][v]
            if c:
                m = adj[v] & ~S
                while m:
                    w = (m & -m).bit_length() - 1; m &= m - 1
                    dp[S | 1 << w][w] += c
    H = sum(dp[full])
    dp = [[0]*n for _ in range(1 << n)]; dp[1][0] = 1
    for S in range(1, full + 1, 2):
        for v in range(n):
            c = dp[S][v]
            if c:
                m = adj[v] & ~S
                while m:
                    w = (m & -m).bit_length() - 1; m &= m - 1
                    dp[S | 1 << w][w] += c
    hc = sum(dp[full][v] for v in range(1, n) if adj[v] & 1)
    return H, hc
pairs = {"163": "001100000101000110101", "401": "101000100101100110101",
         "19": "000000001001000010101", "64": "000001100001010100101",
         "23": "000000010000100010101", "57": "000001001000010000101"}
r = {k: H_hc(adj_from(v)) for k, v in pairs.items()}
for k, (H, hc) in r.items(): print(f"catalogue {k}: H={H} hc={hc} H-2hc={H-2*hc}")
print("5-deck pair 163/401: H-2hc", r["163"][0]-2*r["163"][1], r["401"][0]-2*r["401"][1])
s1 = sum(r[k][0]-2*r[k][1] for k in ("19","64")); s2 = sum(r[k][0]-2*r[k][1] for k in ("23","57"))
print("6-deck affine relation D(19)+D(64)=D(23)+D(57): sums of H-2hc", s1, s2)
