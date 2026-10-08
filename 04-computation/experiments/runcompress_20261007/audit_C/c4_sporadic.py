#!/usr/bin/env python3
"""audit_C item 4 (supplement): equal Terras offset WITHOUT equal odd count.
If sigma_T(M_K) - K = sigma_T(M_K') - K' then the aligned orbits first reach 1 at the same time, hence coincide from the
value 8 on (1 <- 2 <- 4 <- 8 is forced for a first arrival), with debt o(M_K) - o(M_K'). If that debt is nonzero this
is a SPORADIC (nonzero-debt) coincidence: the orbits meet, but it is not a pair-chain absorption (not a deletion
certificate in the sense of HYP-9242). Count such K' for orphans and non-orphans (odd K <= 12800), and verify a few
directly by lockstep."""
from collections import defaultdict
data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split()); data[K] = (o, t)
by_off = defaultdict(list)     # Terras offset sigma_T - K  -> list of K
spor = {}
orphan = {}
for K in range(2, 12801):
    o, t = data[K]
    same_off = by_off[t - K]
    zero = [Kp for Kp in same_off if data[Kp][0] == o]
    nonzero = [Kp for Kp in same_off if data[Kp][0] != o]
    if K % 2 == 1 and K >= 3:
        orphan[K] = not zero
        spor[K] = nonzero
    by_off[t - K].append(K)
orph = [K for K in orphan if orphan[K]]
print("orphans (odd K <= 12800):", len(orph))
hit = [(K, spor[K][-1], data[K][0] - data[spor[K][-1]][0]) for K in orph if spor[K]]
print("orphans with a sporadic equal-time coincidence (equal Terras offset, unequal odd count):", len(hit), hit[:20])
nn = sum(1 for K in orphan if not orphan[K] and spor[K])
print("non-orphans that also have a sporadic coincidence with some other K':", nn)
def lockstep(K, D):
    Kp = K - D
    a = 3 ** D * (1 << Kp) - 1; b = (1 << Kp) - 1; k = D; s = 0
    while a != b:
        if a == 1 or b == 1: return None
        k += (a & 1) - (b & 1)
        a = (3*a + 1) >> 1 if a & 1 else a >> 1
        b = (3*b + 1) >> 1 if b & 1 else b >> 1
        s += 1
    return s, k, a
for K, Kp, dk in hit[:6]:
    print("   direct check", K, Kp, "->", lockstep(K, K - Kp), "(time, debt, value)")
