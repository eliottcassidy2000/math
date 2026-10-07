# THM-4580 (1),(2),(6) brute-force checks
# (1) lift lemma, k <= 6: for n >= k+1, 2^(n+jT_k), j=0..4 share last k digits; digit k+1 runs over the 5 digits of one parity
ok1 = True
for k in range(1, 7):
    T = 4 * 5**(k-1); M = 10**(k+1)
    # LTE: v5(2^T - 1) == k
    v = 2**T - 1; e = 0
    while v % 5 == 0: v //= 5; e += 1
    ok1 &= (e == k)
    for n in range(k + 1, k + 1 + T):
        xs = [pow(2, n + j*T, M) for j in range(5)]
        ok1 &= len(set(x % 10**k for x in xs)) == 1
        ds = sorted(x // 10**k for x in xs)
        ok1 &= ds in ([0, 2, 4, 6, 8], [1, 3, 5, 7, 9])
        # parity of digit k+1 equals parity of q = (2^n mod 10^k)/2^k
        q = (pow(2, n, 10**k)) // 2**k
        ok1 &= (ds[0] % 2) == (q % 2) and (pow(2, n, 10**k) % 2**k == 0)
    # n = k (boundary): does the lemma fail? (the THM requires n >= k+1)
print("(1) lift lemma + LTE for k<=6:", ok1)
# (2) tails of 2^n (n>=k) = k-digit strings divisible by 2^k, prime to 5; Z_k = #zeroless multiples of 2^k
ok2 = True
for k in range(1, 7):
    T = 4 * 5**(k-1)
    tails = set(pow(2, n, 10**k) for n in range(k, k + T))
    target = set(r for r in range(10**k) if r % 2**k == 0 and r % 5 != 0)
    ok2 &= (tails == target) and len(tails) == T
    zl = sum(1 for r in tails if '0' not in str(r).zfill(k))
    zl2 = sum(1 for r in range(0, 10**k, 2**k) if '0' not in str(r).zfill(k))
    ok2 &= (zl == zl2)
print("(2) bijection for k<=6:", ok2)
# (6) doubling criterion for zeroless x < 10^6: 2x zeroless iff no block 5[1-4] and not ending in 5
ok6 = True
for x in range(1, 10**6):
    s = str(x)
    if '0' in s: continue
    crit = not any(s[i] == '5' and s[i+1] in '1234' for i in range(len(s)-1)) and s[-1] != '5'
    ok6 &= (('0' not in str(2*x)) == crit)
print("(6) doubling criterion for zeroless x < 1e6:", ok6)
