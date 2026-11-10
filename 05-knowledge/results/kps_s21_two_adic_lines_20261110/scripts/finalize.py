"""finalize.py - kind-pasteur-2026-11-10-S21. Line table from all seed sources:
  seeds_m.txt (direct runs, orders 2..9), dbl_m.txt (doubling vertices, orders 3..8), us_10.txt (unitsets.c on the
  order-9 classes: unit-offset roots of order 10 with root out-degree 2..4).
Writes line_seeds.txt:  c s0 mechanism order code vertex H d E U
  (tail: rooted tournament (code, vertex) with c0 = 1 and U - d = c 2^(s0-1); module: doubling vertex with H - 1 = c 2^s0)
and below_rays.txt: for every odd c <= CMAX and s < s0(c) a certificate that 1 + c 2^s is a path count:
  "census" (a value of H at order <= 9) or "product a b" (a, b values at order <= 9; order join)."""
import sys
CMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 129
S9 = set()
for fn in ["cls_%d.txt" % m for m in range(1, 9)] + ["cls9_9.txt"]:
    for line in open(fn): S9.add(int(line.split()[1]))
def state(code, m, v):
    out = [0] * m; inn = [0] * m; b = 0
    for i in range(m):
        for j in range(i + 1, m):
            if code >> b & 1: out[i] |= 1 << j; inn[j] |= 1 << i
            else: out[j] |= 1 << i; inn[i] |= 1 << j
            b += 1
    full = (1 << m) - 1
    en = [[0] * m for _ in range(full + 1)]; st = [[0] * m for _ in range(full + 1)]
    for x in range(m): en[1 << x][x] = 1; st[1 << x][x] = 1
    for s in range(1, full + 1):
        for x in range(m):
            c = en[s][x]
            if c:
                nx = out[x] & ~s
                while nx:
                    u = (nx & -nx).bit_length() - 1; nx &= nx - 1; en[s | 1 << u][u] += c
            c = st[s][x]
            if c:
                nx = inn[x] & ~s
                while nx:
                    u = (nx & -nx).bit_length() - 1; nx &= nx - 1; st[s | 1 << u][u] += c
    hs = [sum(r) for r in en]; hs[0] = 1
    H = hs[full]; d = hs[full ^ (1 << v)]; E = en[full][v]
    U = sum(st[X][v] * hs[full ^ X] for X in range(1, full + 1) if X >> v & 1)
    return H, d, E, U
best = {}
def offer(c, key, rec):
    if c not in best or key < best[c][0]: best[c] = (key, rec)
for m in range(2, 10):
    for line in open("seeds_%d.txt" % m):
        mm, code, v, H, d, E, U, M = map(int, line.split())
        if M <= 0: continue
        c = M; j = 0
        while c % 2 == 0: c //= 2; j += 1
        offer(c, (j + 1, mm), ("tail", mm, code, v, H, d, E, U))
for m in range(3, 9):
    for line in open("dbl_%d.txt" % m):
        code, H, x, c2, fl = line.split()
        if fl != "1" or int(c2) == 0: continue
        c = int(c2); j = 0
        while c % 2 == 0: c //= 2; j += 1
        offer(c, (j + 1, m), ("module", m, int(code), int(x), int(H), 0, 0, 0))
pend = {}
for line in open("us_10.txt"):
    m, codeB, O, HB, M = map(int, line.split())
    c = M; j = 0
    while c % 2 == 0: c //= 2; j += 1
    key = (j + 1, m)
    if (c not in best or key < best[c][0]) and (c not in pend or key < pend[c][0]): pend[c] = (key, codeB, O)
for c, (key, codeB, O) in pend.items():
    # order-10 code: B on vertices 0..8 (its pairs come first for i<j<=8 in the lexicographic pair order of order 10
    # only after re-indexing), new root v = 9 with N+(v) = O: bit of pair (i,9) is 1 iff i -> 9 iff i not in O.
    mB = 9; A = [[0] * 10 for _ in range(10)]; b = 0
    for i in range(mB):
        for jj in range(i + 1, mB):
            if codeB >> b & 1: A[i][jj] = 1
            else: A[jj][i] = 1
            b += 1
    for i in range(mB):
        if O >> i & 1: A[9][i] = 1
        else: A[i][9] = 1
    code = 0; b = 0
    for i in range(10):
        for jj in range(i + 1, 10):
            if A[i][jj]: code |= 1 << b
            b += 1
    if c <= 4097:
        H, d, E, U = state(code, 10, 9)
        assert H - E + 2 * d - U == 1 and U - d == c * 2 ** (key[0] - 1), (c, H, d, E, U)
        offer(c, key, ("tail", 10, code, 9, H, d, E, U))
    else:
        offer(c, key, ("tail", 10, code, 9, 0, 0, 0, 0))
out = open("line_seeds.txt", "w")
out.write("# c s0 mechanism order code vertex H d E U   (ray: 1 + c 2^s is a path count for all s >= s0; listed for odd c <= 4097)\n")
for c in sorted(best):
    if c > 4097: continue
    key, (mech, m, code, v, H, d, E, U) = best[c]
    out.write("%d %d %s %d %d %d %d %d %d %d\n" % (c, key[0], mech, m, code, v, H, d, E, U))
out.close()
cert = open("below_rays.txt", "w"); cert.write("# c s N certificate   (points 1 + c 2^s below the first ray, odd c <= %d)\n" % CMAX)
bad = []; ncen = 0; nprod = 0
for c in range(1, CMAX + 1, 2):
    if c not in best: bad.append(("no seed", c)); continue
    for s in range(1, best[c][0][0]):
        N = 1 + c * 2 ** s
        if N in (7, 21): cert.write("%d %d %d GAP (THM-1370)\n" % (c, s, N)); continue
        if N in S9: cert.write("%d %d %d census\n" % (c, s, N)); ncen += 1; continue
        f = next(((a, N // a) for a in sorted(S9) if 1 < a < N and N % a == 0 and N // a in S9), None)
        if f: cert.write("%d %d %d product %d %d\n" % (c, s, N, f[0], f[1])); nprod += 1
        else: bad.append((c, s, N))
cert.close()
print("odd c with a seed:", len(best), "| c <= %d without seed or with an uncertified point below the ray:" % CMAX, bad)
print("points below rays for c <= %d: %d by census, %d by product" % (CMAX, ncen, nprod))
print("mechanisms for c <= %d:" % CMAX, {k: sum(1 for c in range(1, CMAX + 1, 2) if c in best and (best[c][1][0], best[c][1][1]) == k) for k in set((best[c][1][0], best[c][1][1]) for c in range(1, CMAX + 1, 2) if c in best)})
print("s0 table c <= 129:", " ".join("%d:%d" % (c, best[c][0][0]) for c in range(1, 130, 2)))
