"""
audit_S25_collatz.py -- independent re-derivation for the audit of THM-4627 / HYP-9258 (session opus-2026-10-13-S25).
Own code; nothing is imported from the session's scripts.

(1) E[s] and the lag correlations of s = (-1)^sigma (true stopping time) at lags 1..8 for n <= 2^NL, with the lag laws
    C(odd h) = -(1 + 2 E[s]),  C(h = 2 mod 4) = +(1 + 2 E[s]),  C(h = 0 mod 4) = 3/4 + (1/4) E[s(n) s(n+h) | n = 3 mod 4]
    checked on the same data (differences are edge effects O(h/N)); also sigma = 1 on even n, sigma = 2 on n = 1 mod 4 (n > 1).
(2) the merge census at lag 1 to N = 2^NM with the rule 'stop when either orbit reaches 1': merged fraction, correlation of
    t = (-1)^tau on all / merged / unmerged pairs, mean merge time.
(3) break density P(tau(n) != tau(n+1)) at 2^16 .. 2^NB, parity-flip share, and the ACTUAL starting integer of the longest run
    below each scale (the session's script printed the run index, not the integer).
(4) periodic (Haar) truncation at depth k = 20 and KP with exact integer comparison 3^a < 2^j: P(kappa > k), E_k[s], C_k(h).
(5) sigma versus kappa on n <= 2^20: how many n have sigma(n) != kappa(n) (the 'finite exceptional set per class').
(6) dyadic-block means of t (cross-check of the parity-bias lane) and the identity E_j = -B_(j-1).
Run: python3 audit_S25_collatz.py [NL=24] [NM=22] [NB=28] [KP=26]
"""
import sys, time
import numpy as np
from numba import njit

NL = int(sys.argv[1]) if len(sys.argv) > 1 else 24
NM = int(sys.argv[2]) if len(sys.argv) > 2 else 22
NB = int(sys.argv[3]) if len(sys.argv) > 3 else 28
KP = int(sys.argv[4]) if len(sys.argv) > 4 else 26
LOG2_3 = 1.5849625007211563


@njit(cache=True)
def sigma_and_tau(N):
    """s[n] = (-1)^sigma(n) with sigma the true stopping time; tau[n] = T-steps to 1 (memoised); n = 1 left at 0."""
    s = np.zeros(N + 1, dtype=np.int8)
    sig = np.zeros(N + 1, dtype=np.int32)
    tau = np.zeros(N + 1, dtype=np.int32)
    for n in range(2, N + 1):
        m = n
        k = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            k += 1
        sig[n] = k
        s[n] = -1 if (k & 1) else 1
        tau[n] = k + tau[m]
    return s, sig, tau


@njit(cache=True)
def corr(x, h, N):
    acc = 0
    for n in range(1, N - h + 1):
        acc += x[n] * x[n + h]
    return acc / (N - h)


@njit(cache=True)
def corr_cond(x, h, N, r, mod):
    acc = 0
    cnt = 0
    for n in range(1, N - h + 1):
        if n % mod == r:
            acc += x[n] * x[n + h]
            cnt += 1
    return acc / cnt


@njit(cache=True)
def mean_cond(x, N, r, mod):
    acc = 0
    cnt = 0
    for n in range(1, N + 1):
        if n % mod == r:
            acc += x[n]
            cnt += 1
    return acc / cnt


@njit(cache=True)
def merge_census(N, h):
    """merge time of the equal-time coincidence of the orbits of n and n+h, or -1 if one orbit reaches 1 first."""
    mt = np.full(N + 1, -1, dtype=np.int32)
    for n in range(1, N + 1):
        x = n
        y = n + h
        k = 0
        while x != y:
            if x == 1 or y == 1:
                k = -1
                break
            x = (3 * x + 1) >> 1 if (x & 1) else (x >> 1)
            y = (3 * y + 1) >> 1 if (y & 1) else (y >> 1)
            k += 1
        mt[n] = k
    return mt


@njit(cache=True)
def tau16(N):
    tau = np.zeros(N + 2, dtype=np.int16)
    for n in range(2, N + 2):
        m = n
        k = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            k += 1
        tau[n] = k + tau[m]
    return tau


@njit(cache=True)
def runs_and_breaks(tau, scales):
    """for each scale M: break density, flip density, longest run (length, start integer) among runs ending at a break n <= M."""
    out = np.zeros((scales.shape[0], 4))
    N = scales[-1]
    best_len = 0
    best_start = 1
    cur_start = 1
    nb = 0
    nf = 0
    si = 0
    for n in range(1, N + 1):
        a = int(tau[n])
        b = int(tau[n + 1])
        if a != b:
            nb += 1
            if (a - b) & 1:
                nf += 1
            L = n - cur_start + 1
            if L > best_len:
                best_len = L
                best_start = cur_start
            cur_start = n + 1
        if n == scales[si]:
            out[si, 0] = nb / n
            out[si, 1] = nf / n
            out[si, 2] = best_len
            out[si, 3] = best_start
            si += 1
            if si >= scales.shape[0]:
                break
    return out


@njit(cache=True)
def kappa_all(k):
    """coefficient stopping time of every residue class mod 2^k, exact integers (0 means > k)."""
    M = 1 << k
    out = np.zeros(M, dtype=np.int8)
    for r in range(M):
        n = r
        p3 = 1
        p2 = 1
        kap = 0
        for j in range(1, k + 1):
            if n & 1:
                n = (3 * n + 1) >> 1
                p3 *= 3
            else:
                n >>= 1
            p2 *= 2
            if p3 < p2:
                kap = j
                break
        out[r] = kap
    return out


@njit(cache=True)
def corr_periodic(s, h):
    M = s.shape[0]
    acc = 0
    for r in range(M):
        acc += s[r] * s[(r + h) & (M - 1)]
    return acc / M


@njit(cache=True)
def sigma_vs_kappa(N):
    """count n <= N with sigma(n) != kappa(n); return count, the largest such n, and the first ten."""
    cnt = 0
    largest = 0
    firsts = np.zeros(10, dtype=np.int64)
    for n in range(2, N + 1):
        m = n
        k = 0
        a = 0
        kap = 0
        sig = 0
        while True:
            if m & 1:
                m = (3 * m + 1) >> 1
                a += 1
            else:
                m >>= 1
            k += 1
            if kap == 0 and a * LOG2_3 < k:
                kap = k
            if sig == 0 and m < n:
                sig = k
            if kap > 0 and sig > 0:
                break
        if kap != sig:
            if cnt < 10:
                firsts[cnt] = n
            cnt += 1
            largest = n
    return cnt, largest, firsts


t0 = time.time()
N = 1 << NL
s, sig, tau = sigma_and_tau(N + 64)
print(f'# (1) stopping-time parity s = (-1)^sigma, n <= 2^{NL} ({time.time()-t0:.1f}s)')
# structural checks
ev = sig[2:N + 1:2]
print(f'   sigma(n) = 1 for every even n <= N: {bool(np.all(ev == 1))}')
m1 = sig[5:N + 1:4]
print(f'   sigma(n) = 2 for every n = 1 mod 4, 5 <= n <= N: {bool(np.all(m1 == 2))}; sigma(1) is undefined (orbit never drops below 1)')
m3 = sig[3:N + 1:4]
print(f'   min sigma on n = 3 mod 4: {int(m3.min())} (should be >= 3)')
Es = float(np.mean(s[1:N + 1]))
m_odd = float(mean_cond(s, N, 1, 2))
m_3 = float(mean_cond(s, N, 3, 4))
print(f'   E[s] = {Es:+.6f};  1 + 2E[s] = {1 + 2 * Es:+.6f};  E[s | odd] = {m_odd:+.6f};  E[s | 3 mod 4] = {m_3:+.6f};  2 m_odd - 1 = {2 * m_odd - 1:+.6f}')
print('   h   C(h)        law value     C - law     (law: -(1+2E[s]) odd; +(1+2E[s]) 2 mod 4; 3/4 + E[ss|3 mod 4]/4 for 0 mod 4)')
for h in range(1, 9):
    c = float(corr(s, h, N))
    if h % 2 == 1:
        law = -(1 + 2 * Es)
    elif h % 4 == 2:
        law = 1 + 2 * Es
    else:
        law = 0.75 + 0.25 * float(corr_cond(s, h, N, 3, 4))
    print(f'   {h}  {c:+.6f}   {law:+.6f}   {c - law:+.2e}')
for NBk in (20, 22):
    Nb = 1 << NBk
    print(f'   N=2^{NBk}: E[s] = {float(np.mean(s[1:Nb+1])):+.6f}  C(1) = {float(corr(s, 1, Nb)):+.6f}  C(2) = {float(corr(s, 2, Nb)):+.6f}  C(4) = {float(corr(s, 4, Nb)):+.6f}')

t_par = np.where(tau % 2 == 1, -1, 1).astype(np.int8)
t_par[0] = 0
print(f'\n# (2) merge census at lag 1, n <= 2^{NM}, stop when either orbit reaches 1')
Nm = 1 << NM
t0 = time.time()
mt = merge_census(Nm, 1)
merged = mt[1:Nm + 1] >= 0
prod = (t_par[1:Nm + 1].astype(np.int32) * t_par[2:Nm + 2].astype(np.int32))
pm = float(merged.mean())
print(f'   P(merge before 1) = {pm:.5f};  corr(all) = {float(prod.mean()):+.5f};  corr(merged) = {float(prod[merged].mean()):+.5f};  '
      f'corr(unmerged) = {float(prod[~merged].mean()):+.5f};  tau equal on merged: {float((tau[1:Nm+1] == tau[2:Nm+2])[merged].mean()):.4f};  '
      f'mean merge time {float(mt[1:Nm+1][merged].mean()):.2f}  ({time.time()-t0:.1f}s)')
for NBk in (20,):
    Nb = 1 << NBk
    print(f'   N=2^{NBk}: P(merge) = {float(merged[:Nb].mean()):.5f}  corr = {float(prod[:Nb].mean()):+.5f}')
eq_unm = float((tau[1:Nm + 1] == tau[2:Nm + 2])[~merged].mean())
print(f'   fraction of UNMERGED pairs with equal tau: {eq_unm:.6f}  (so break density = 1 - P(merge) up to this)')
print(f'   break density at 2^{NM} from tau: {float((tau[1:Nm+1] != tau[2:Nm+2]).mean()):.5f};  1 - P(merge) = {1-pm:.5f}')
del s, sig, tau, t_par, mt, merged, prod

print(f'\n# (3) break density, parity flips, longest run with its ACTUAL start, scales 2^16..2^{NB} (int16 tau)')
t0 = time.time()
tb = tau16(1 << NB)
scales = np.array([1 << L for L in range(16, NB + 1, 2)], dtype=np.int64)
res = runs_and_breaks(tb, scales)
print(f'   tau to 2^{NB} in {time.time()-t0:.1f}s; max tau {int(tb.max())}')
print('   scale   break density  flip density  flip/break  mean run   longest run  starts at n')
for i, L in enumerate(range(16, NB + 1, 2)):
    bd, fd, bl, bs = res[i]
    print(f'   2^{L:2d}   {bd:.6f}       {fd:.6f}      {fd/bd:.4f}     {1/bd:.4f}    {int(bl):4d}        {int(bs)}')
bl, bs = int(res[-1, 2]), int(res[-1, 3])
vals = tb[bs:bs + bl + 1]
print(f'   check: tau constant on [{bs}, {bs + bl - 1}]: {bool(np.all(vals[:-1] == vals[0]))}, tau = {int(vals[0])}, tau({bs - 1}) = {int(tb[bs - 1])}, tau({bs + bl}) = {int(vals[-1])}')
print('\n# (6) dyadic-block means B_j of t = (-1)^tau over [2^(j-1), 2^j), even half E_j, odd half O_j, check E_j + B_(j-1) = 0')
tp = np.where(tb % 2 == 1, -1, 1).astype(np.int8)
tp[0] = 0
prevB = None
for j in range(16, NB + 1):
    lo, hi = 1 << (j - 1), 1 << j
    blk = tp[lo:hi]
    B = float(blk.mean())
    E = float(blk[0::2].mean())
    O = float(blk[1::2].mean())
    chk = f'{E + prevB:+.2e}' if prevB is not None else 'n/a'
    print(f'   2^{j:2d}: B = {B:+.6f}  E = {E:+.6f}  O = {O:+.6f}  E_j + B_(j-1) = {chk}  cumulative mean to 2^{j}: {float(tp[1:hi].mean()):+.6f}')
    prevB = B
del tb, tp

print(f'\n# (4) periodic (Haar) truncation of s at depth k, exact integer test 3^a < 2^j')
for k in (20, KP):
    t0 = time.time()
    kap = kappa_all(k)
    undet = float(np.mean(kap == 0))
    sk = np.where(kap == 0, 0, np.where(kap % 2 == 1, -1, 1)).astype(np.int8)
    Ek = float(np.mean(sk))
    cs = [float(corr_periodic(sk, h)) for h in (1, 2, 4, 8, 16, 32)]
    print(f'   k={k}: P(kappa > k) = {undet:.6f}  E_k[s] = {Ek:+.6f}  C_k(1,2,4,8,16,32) = ' + ' '.join(f'{c:+.5f}' for c in cs) + f'  ({time.time()-t0:.1f}s)')
    print(f'         lag-law check on the truncated periodic data: -(1+2E_k[s]) = {-(1+2*Ek):+.5f} vs C_k(1) = {cs[0]:+.5f} (the truncation zeros break the exact law; difference = mass of undetermined odd classes)')

print('\n# (5) sigma versus kappa on n <= 2^20')
cnt, largest, firsts = sigma_vs_kappa(1 << 20)
print(f'   #{{n <= 2^20 : sigma(n) != kappa(n)}} = {cnt}; largest such n = {largest}; first ten: {list(firsts[:min(cnt,10)])}')
print('done')
