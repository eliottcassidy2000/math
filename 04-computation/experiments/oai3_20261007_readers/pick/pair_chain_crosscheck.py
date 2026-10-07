"""Cross-check: chain merge time == direct equal-time, equal-odd-count merge of T-orbits of p and p-1
(the lag-1 Mersenne switch U^i(p) = U^(i+1)(q), q = (p-2)/3), on random 400-bit p = 1 mod 16,
and exhaustive check of P(merged by n=20) over all p mod 2^20, p = 1 mod 16."""
import random
from pair_chain_verify import step, T

def chain_merge_time(p, nmax):
    v = p - 1; k, num = 0, 1
    for n in range(nmax):
        k, num, _ = step(k, num, v % 2)
        v = T(v)
        if (k, num) == (0, 0):
            return n + 1
    return None

def direct_merge_time(p, nmax):
    u, v = p, p - 1; cu = cv = 0
    for n in range(nmax):
        cu += u % 2; cv += v % 2
        u, v = T(u), T(v)
        if u == v and cu == cv:
            return n + 1
    return None

rng = random.Random(3)
agree = 0; tot = 0; merged = 0
for _ in range(3000):
    p = (rng.getrandbits(400) << 4) | 1
    a = chain_merge_time(p, 300); b = direct_merge_time(p, 300)
    tot += 1; agree += (a == b); merged += (a is not None)
print(f"chain vs direct merge time agree on {agree}/{tot} random 404-bit p (merged by 300: {merged})")
# lag-1 Syracuse check: U^i(p) = U^(i+1)(q) with q = (p-2)/3 when merged
def U(x):
    x = 3 * x + 1
    while x % 2 == 0: x //= 2
    return x
ok = 0; cnt = 0
for _ in range(300):
    p = (rng.getrandbits(400) << 4) | 1
    if (p - 2) % 3: continue
    q = (p - 2) // 3
    if q % 2 == 0: continue
    t = chain_merge_time(p, 300)
    if t is None: continue
    cnt += 1
    # find i with U^i(p) == U^(i+1)(q)
    a, b = p, U(q); found = False
    for i in range(400):
        if a == b: found = True; break
        a, b = U(a), U(b)
    ok += found
print(f"lag-1 Syracuse merge U^i(p)=U^(i+1)(q) found for {ok}/{cnt} chain-merged samples")
# exhaustive n = 20
n = 20; hits = 0; total = 0
for p in range(1, 2 ** n, 16):
    total += 1
    t = chain_merge_time(p + 2 ** 30 * 0, n)  # p mod 2^n representative; chain uses first n parities only
    hits += (t is not None)
print(f"exhaustive: P(merged by n=20) = {hits}/{total} = {hits/total:.6f}  (DP gave 0.121429)")
