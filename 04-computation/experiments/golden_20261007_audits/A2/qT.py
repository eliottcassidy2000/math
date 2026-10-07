# Monte Carlo of q(T) = P(chain from (0,1) not absorbed by T), integer representation E = 3^max(0,-k) e
import random
random.seed(7)
def run(T):
    k, E = 0, 1          # e = E / 3^max(0,-k)
    for t in range(T):
        b = random.getrandbits(1)
        s = E & 1        # parity of e equals parity of E (3^h odd)
        if k >= 0:
            if s == 0 and b == 0: E //= 2
            elif s == 0 and b == 1: E = (3 * E + 1 - 3 ** k) // 2
            elif s == 1 and b == 0: E = (3 * E + 1) // 2; k += 1
            else:
                if k >= 1: E = (E - 3 ** (k - 1)) // 2; k -= 1
                else:      E = (3 * E - 1) // 2; k = -1          # k=0 -> -1: E' = 3 e' = (3e - 1)/2
        else:
            h = -k
            if s == 0 and b == 0: E //= 2
            elif s == 0 and b == 1: E = (3 * E + 3 ** h - 1) // 2
            elif s == 1 and b == 0: E = (E + 3 ** (h - 1)) // 2; k += 1
            else: E = (3 * E - 1) // 2; k -= 1
        if k == 0 and E == 0: return t + 1
    return None
N = 40000; Ts = [20, 48, 72, 96, 134, 200]
cnt = {T: 0 for T in Ts}
for i in range(N):
    tau = run(max(Ts))
    for T in Ts:
        if tau is None or tau > T: cnt[T] += 1
for T in Ts:
    q = cnt[T] / N
    print(f"q({T}) = {q:.4f} +- {((q*(1-q)/N)**0.5):.4f}   0.666*q = {0.666*q:.3f}")
