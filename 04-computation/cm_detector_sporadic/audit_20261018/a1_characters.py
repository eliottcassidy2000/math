"""Audit check 1 (THM-4708 part 1): product of all (2,2K+1) Virasoro characters.

Characters from the Rocha-Caridi formula (NOT from Andrews-Gordon products):
  chi_{r,s}^{(p,p')} = q^{h_{r,s} - c/24} / prod(1-q^n) * sum_{k in Z} ( q^{pp'k^2 + k(p'r - ps)} - q^{pp'k^2 + k(p'r + ps) + rs} ),
  h_{r,s} = ((p'r - ps)^2 - (p' - p)^2) / (4pp'),  c = 1 - 6 (p'-p)^2/(pp').
Here (p,p') = (2, 2K+1), r = 1, s = 1..K (s and 2K+1-s are identified).
Target: prod_s chi_{1,s} = (eta((2K+1)tau)/eta(tau))^{K-1}, i.e. leading power K(K-1)/12 and the integer series
prod_n ((1-q^{(2K+1)n})/(1-q^n))^{K-1}.
"""
from fractions import Fraction as Fr
import sys

ORDER = int(sys.argv[1]) if len(sys.argv) > 1 else 300


def mul(a, b, N):
    c = [0] * N
    for i, x in enumerate(a[:N]):
        if x:
            for j, y in enumerate(b[:N - i]):
                c[i + j] += x * y
    return c


def inv_phi(N):  # 1/prod(1-q^n) = partition numbers
    p = [0] * N
    p[0] = 1
    for n in range(1, N):
        for k in range(n, N):
            p[k] += p[k - n]
    return p


def phi_pow(step, e, N):  # prod_n (1-q^{step n})^e, e integer (may be negative)
    s = [0] * N
    s[0] = 1
    for n in range(1, N):
        if step * n >= N:
            break
        d = step * n
        if e >= 0:
            for _ in range(e):
                for k in range(N - 1, d - 1, -1):
                    s[k] -= s[k - d]
        else:
            for _ in range(-e):
                for k in range(d, N):
                    s[k] += s[k - d]
    return s


def rocha_caridi(p, pp, r, s, N):
    h = Fr((pp * r - p * s) ** 2 - (pp - p) ** 2, 4 * p * pp)
    c = 1 - Fr(6 * (pp - p) ** 2, p * pp)
    theta = [0] * N
    K = 0
    while True:
        done = True
        for k in (K, -K) if K else (0,):
            e1 = p * pp * k * k + k * (pp * r - p * s)
            e2 = p * pp * k * k + k * (pp * r + p * s) + r * s
            if 0 <= e1 < N:
                theta[e1] += 1; done = False
            if 0 <= e2 < N:
                theta[e2] -= 1; done = False
            if e1 < 0 or e2 < 0:
                raise ValueError('negative exponent')
        if done and K > 2:
            break
        K += 1
    return h - c / 24, mul(theta, inv_phi(N), N), h, c


ok_all = True
for K in range(2, 9):
    p, pp = 2, 2 * K + 1
    N = ORDER
    lead = Fr(0)
    prod = [1] + [0] * (N - 1)
    hs = []
    for s in range(1, K + 1):
        e, ser, h, c = rocha_caridi(p, pp, 1, s, N)
        hs.append(h)
        lead += e
        prod = mul(prod, ser, N)
    target_lead = Fr((K - 1) * (pp - 1), 24)
    target = mul(phi_pow(pp, K - 1, N), phi_pow(1, -(K - 1), N), N)
    ok = (lead == target_lead) and prod == target and lead == Fr(K * (K - 1), 12)
    # also: Andrews-Gordon product form of each character (independent cross-check of the proof's input)
    ag_ok = True
    for s in range(1, K + 1):
        e, ser, h, c = rocha_caridi(p, pp, 1, s, N)
        agp = [1] + [0] * (N - 1)
        for n in range(1, N):
            if n % pp not in (0, s, pp - s):
                for k in range(n, N):
                    agp[k] += agp[k - n]
        ag_ok &= (agp == ser)
    ok_all &= ok and ag_ok
    print(f'K={K} (2,{pp}): c={c}, h_1s={[str(x) for x in hs]}, leading exponent {lead} '
          f'(target (K-1)(2K)/24 = {target_lead}); series equal to O(q^{N}): {prod == target}; '
          f'Andrews-Gordon products = Rocha-Caridi for all s: {ag_ok}')
print('ALL OK' if ok_all else 'FAILURE')
