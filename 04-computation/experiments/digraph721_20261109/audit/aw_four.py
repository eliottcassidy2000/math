"""AUDIT D1 (THM-5055 Remark 7.1): independent re-implementation of the four-move closure (padding, cone over W,
E_2-join, cone over G u K_1, R) from the auditor's exhaustive base (aw_enum tuples, n <= 9), big-int bitsets.
Checks: (a) E(n) = F(n) for n <= NM; (b) the value v_n = 2^(n-1) - 81 lies in S_n (pendant with non-universal
neighbour) exactly for which n; (c) the pendant bound: every value in P_m is <= 3 * 2^(m-2).
usage: python3 aw_four.py base_tuples NB NM"""
import sys
base, NB, NM = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
E = [0] * (NM + 1); W = [0] * (NM + 1); P = [0] * (NM + 1); S = [0] * (NM + 1)
for line in open(base):
    p = line.split('|')[0].split()
    if len(p) < 7: continue
    n, d, U, pend, pnu = map(int, p[:5])
    if n > NB: continue
    E[n] |= 1 << d
    if not U: W[n] |= 1 << d
    if pend: P[n] |= 1 << d
    if pnu: S[n] |= 1 << d
def lim(n): return (1 << (1 << n)) - 1
def sh(x, off, n): return (x << off) & lim(n)
for n in range(2, NM + 1):
    E[n] |= E[n - 1]; W[n] |= E[n - 1]; P[n] |= P[n - 1]; S[n] |= S[n - 1] | P[n - 1]
    if n > NB:
        S[n] |= sh(S[n - 2], 3 << (n - 3), n)
        P[n] |= S[n] | sh(E[n - 2], 1 << (n - 1), n)
        W[n] |= sh(W[n - 2], (3 << (n - 2)) - 2, n) | S[n]
        E[n] |= W[n] | sh(W[n - 1], 1 << (n - 1), n) | sh(E[n - 2], (3 << (n - 2)) - 2, n)
def M(n): return 2**n - n if n % 2 else 2**n - n - 1
tot = 0
for n in range(1, NM + 1):
    Fm = int.from_bytes(b'\xaa' * ((M(n) >> 3) + 1), 'little') & ((1 << (M(n) + 1)) - 1)   # odd k <= M(n)
    for k in (3, 7, 13): Fm &= ~(1 << k)
    if n % 2 == 1 and n >= 5: Fm &= ~(1 << (2**n - n - 2))
    miss = bin(Fm & ~E[n]).count('1'); extra = bin(E[n] & ~Fm).count('1'); tot += miss + extra
    pb = P[n] >> (3 * 2**(n - 2) + 1) if n >= 2 else 0
    v = 2**(n - 1) - 81
    print(f"n={n}: missing {miss} extra {extra}; max P value <= 3*2^(n-2): {pb == 0}; "
          f"v_n = 2^(n-1)-81 = {v}: in S_n {bool(v > 0 and (S[n] >> v) & 1)}, in E_n {bool(v > 0 and (E[n] >> v) & 1)}")
print("total mismatches", tot)
