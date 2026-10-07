# SRW persistence: sqrt(T) * P(S_n > -j for n <= T) -> j sqrt(2/pi); reflection: P = P(-j < S_T <= j)
from math import lgamma, exp, log, sqrt, pi
def pmf(T, k):
    if (T + k) % 2 or abs(k) > T: return 0.0
    m = (T + k) // 2
    return exp(lgamma(T + 1) - lgamma(m + 1) - lgamma(T - m + 1) - T * log(2))
for j in (1, 2):
    for T in (10**4, 10**6, 10**8):
        p = sum(pmf(T, k) for k in range(-j + 1, j + 1))
        print(f"j={j} T={T}: sqrt(T)*P = {sqrt(T)*p:.5f}   j*sqrt(2/pi) = {j*sqrt(2/pi):.5f}")
