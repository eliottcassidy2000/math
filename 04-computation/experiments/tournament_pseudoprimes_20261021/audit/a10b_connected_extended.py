# AUDIT a10b (independent auditor, 2026-10-21): connected labelled graphs mod p, extended range.
# Claim tested: cg(m+p-1) = (2/p) cg(m) (mod p) for ALL m >= 2 (fails only at m = 1), every odd prime 5 <= p <= 101,
# 1 <= m <= 8p.  Proof (auditor): shift operator Delta: (a_m) -> (a_{m+1}) maps E_x to x E_x, and the rooted recurrence
# says Delta(cg) = Delta(T) * T^{-1} in the Hurwitz ring; T = iota(t) => Delta(cg) = iota(t' t^{-1}), t' = sum c_x x [x].
from math import comb
def conn_mod(M, p):
    g = [pow(2, j * (j - 1) // 2, p) for j in range(M + 1)]; c = [0] * (M + 1)
    row = [1]  # C(m-1, .)
    for m in range(1, M + 1):
        if m > 1:
            row = [1] + [(row[i - 1] + row[i]) % p for i in range(1, m - 1)] + [1]
        c[m] = (g[m] - sum(row[k - 1] * c[k] * g[m - k] for k in range(1, m))) % p
    return c
allok = True
for p in [q for q in range(5, 102) if all(q % d for d in range(2, int(q ** .5) + 1))]:
    M = 8 * p; cg = conn_mod(M + p, p); e = 1 if pow(2, (p - 1) // 2, p) == 1 else -1
    fails = [m for m in range(1, M + 1) if (cg[m + p - 1] - e * cg[m]) % p]
    allok &= (fails == [1])
    print(f"p={p} eps={e:+d}: m in [1,{M}] where cg(m+p-1) != eps cg(m): {fails}")
print("CONNECTED GRAPHS ARE QUASI-PERIODIC FOR m >= 2 AT EVERY TESTED p" if allok else "PATTERN BROKEN")
