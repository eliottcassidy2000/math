import random, time
W = 300; M = 10**W
t0 = time.time()
# (1) all n in [0, 10^6]: exact for n < 1000, then doubling mod 10^300 (2^n >= 10^300 for n >= 997)
zeroless = []
for n in range(0, 1000):
    if '0' not in str(2**n): zeroless.append(n)
print("zeroless 2^n for n < 1000:", len(zeroless), "largest", zeroless[-1], zeroless == [0,1,2,3,4,5,6,7,8,9,13,14,15,16,18,19,24,25,27,28,31,32,33,34,35,36,37,39,49,51,67,72,76,77,81,86])
x = pow(2, 957, M)
best = -1; recs = []; maxpos = 0; fails = []
for n in range(957, 10**6 + 1):
    s = str(x).zfill(W)
    k = s.rfind('0')
    zl = W - 1 - k if k >= 0 else W   # number of zeroless trailing digits
    if zl >= 251: fails.append(n)
    if zl > best: best = zl; recs.append((n, zl))
    x = (2 * x) % M
print("n in [957,1e6]: zero within last 251 digits fails at", fails[:5], "; records (n, zeroless suffix):", recs)
print("check 2^1e6 state:", x == pow(2, 10**6 + 1, M))
# (2) random sample in [957, 1.1e11)
rng = random.Random(4580)
bad = []; mx = (0, 0)
for _ in range(3000):
    n = rng.randrange(957, 11 * 10**10)
    s = str(pow(2, n, M)).zfill(W)
    zl = W - 1 - s.rfind('0')
    if zl >= 251: bad.append(n)
    if zl > mx[0]: mx = (zl, n)
# extra sample right below 1.1e11 and around the records
for n in list(range(11 * 10**10 - 2000, 11 * 10**10)) + list(range(109171987836 - 500, 109171987836 + 500)):
    s = str(pow(2, n, M)).zfill(W)
    zl = W - 1 - s.rfind('0')
    if zl >= 251: bad.append(n)
print("random/targeted sample: violations (>=251 zeroless trailing digits):", bad, " max zeroless suffix seen", mx)
print("R1 final low 18 digits check: 2^1e10 mod 10^18 =", pow(2, 10**10, 10**18), "== 374549681787109376:", pow(2, 10**10, 10**18) == 374549681787109376)
print("elapsed %.1fs" % (time.time() - t0))
