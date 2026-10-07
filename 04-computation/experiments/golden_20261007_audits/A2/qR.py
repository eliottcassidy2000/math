import random, sys, math
def stepE(k, E, b):
    s = E & 1
    if k >= 0:
        if s == 0 and b == 0: return k, E >> 1
        if s == 0: return k, (3 * E + 1 - 3 ** k) >> 1
        if b == 0: return k + 1, (3 * E + 1) >> 1
        if k >= 1: return k - 1, (E - 3 ** (k - 1)) >> 1
        return -1, (3 * E - 1) >> 1
    h = -k
    if s == 0 and b == 0: return k, E >> 1
    if s == 0: return k, (3 * E + 3 ** h - 1) >> 1
    if b == 0: return k + 1, (E + 3 ** (h - 1)) >> 1
    return k - 1, (3 * E - 1) >> 1
R = int(sys.argv[1]); N = int(sys.argv[2]); TMAX = int(sys.argv[3]); random.seed(int(sys.argv[4]))
Ts = [100, 1000, 3000, 10000, 30000, 100000]
Ts = [t for t in Ts if t <= TMAX]
surv = {t: 0 for t in Ts}
for i in range(N):
    states = {(0, r) for r in range(1, R + 1)}
    tau = None
    for t in range(TMAX):
        b = random.getrandbits(1)
        states = {stepE(k, E, b) for (k, E) in states}
        if (0, 0) in states: tau = t + 1; break
    for T in Ts:
        if tau is None or tau > T: surv[T] += 1
for T in Ts:
    q = surv[T] / N
    print(f"R={R} T={T:6d} q={q:.5f} sqrtT*q={math.sqrt(T)*q:6.2f} +- {math.sqrt(T)*math.sqrt(q*(1-q)/N):.2f}  R*sqrtT*q={R*math.sqrt(T)*q:6.2f}")
