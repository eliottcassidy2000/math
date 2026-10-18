#!/usr/bin/env python3
"""Extended check of the CM harmonic log law (L_p) at the three open reflection classes (S43, HYP-9310).

(L_p): 6 sum_{k<=n} C(2k,k)^3 m^-k (H_2k - H_k) == q_p(m) sum_{k<=n} C(2k,k)^3 m^-k (mod p), n = (p-1)/2.
By the reflection lemma of THM-4685 the law at m and at 4096/m are equivalent, and it is proved at m = +-64; the open
classes are represented by m = 1 (Q(sqrt-7)), m = -8 (Q(i)), m = 16 (Q(sqrt-3)).  Also checks the class-6 congruence (K).
Inverses mod p by the recurrence inv[i] = -(p // i) inv[p % i].
Usage: python cm_harmonic_log_law_extended_20261018.py [PMIN] [PMAX]
Session opus-2026-10-17-S43.
"""
import os, sys, json


def primes_between(a, b):
    s = bytearray([1]) * (b + 1)
    s[0] = s[1] = 0
    for i in range(2, int(b ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(max(a, 5), b + 1) if s[i]]


def main():
    PMIN = int(sys.argv[1]) if len(sys.argv) > 1 else 5
    PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    ms = [1, -8, 16]
    fails = {str(m): [] for m in ms}
    kfail = []
    count = 0
    for p in primes_between(PMIN, PMAX):
        n = (p - 1) // 2
        inv = [0, 1] + [0] * (p - 2)
        for i in range(2, p):
            inv[i] = (-(p // i) * inv[p % i]) % p
        H = [0] * (2 * n + 1)
        for i in range(1, 2 * n + 1):
            H[i] = (H[i - 1] + inv[i]) % p
        cb = [1] * (n + 1)
        for k in range(1, n + 1):
            cb[k] = cb[k - 1] * (2 * k) % p * (2 * k - 1) % p * inv[k] % p * inv[k] % p
        c3 = [pow(c, 3, p) for c in cb]
        dH = [(H[2 * k] - H[k]) % p for k in range(n + 1)]
        for m in ms:
            if m % p == 0:
                continue
            im = pow(m % p, p - 2, p)
            F = 0
            h = 0
            w = 1
            for k in range(n + 1):
                t = c3[k] * w % p
                F += t
                h += t * dH[k]
                w = w * im % p
            q = (pow(m % (p * p), p - 1, p * p) - 1) // p % p
            if (6 * h - q * F) % p:
                fails[str(m)].append(p)
        # (K): sum_{j<=n/2} a_j^3 a_2j (H_2j - 2H_j - 2H_4j + 2H_(n-j) + H_(n-2j)) == 0, a_k = C(2k,k)/4^k
        inv4 = inv[4]
        a = [cb[k] * pow(inv4, k, p) % p for k in range(n + 1)]
        sK = 0
        for j in range(n // 2 + 1):
            br = (H[2 * j] - 2 * H[j] - 2 * H[4 * j] + 2 * H[n - j] + H[n - 2 * j]) % p
            sK += pow(a[j], 3, p) * a[2 * j] % p * br
        if sK % p:
            kfail.append(p)
        count += 1
    print('primes %d..%d (%d primes):' % (PMIN, PMAX, count))
    for m in ms:
        print('  (L_p) at m = %3d: failures %s' % (m, fails[str(m)][:10] if fails[str(m)] else 'none'))
    print('  (K) class-6 congruence: failures %s' % (kfail[:10] if kfail else 'none'))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(pmin=PMIN, pmax=PMAX, primes=count, fails=fails, K_fail=kfail),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'cm_harmonic_log_law_extended_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
