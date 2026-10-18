"""Audit falsification attempt for the general weight dichotomy of HYP-9330 (title): (B) in weight 1 and
(C) [rad F_p | V_p] in weight 2, on families NOT tested by the author.
weight 1 (2nd order): C(2k,k)^2, (3k)!/k!^3, (4k)!/((2k)! k!^2), (6k)!/((3k)!(2k)! k!)
weight 2 (3rd order): C(2k,k)^3, C(2k,k)^2 C(3k,k), C(2k,k)^2 C(4k,2k), C(2k,k) C(3k,k) C(6k,3k),
                      Almkvist-Zudilin sum_k (-1)^k 3^{n-3k} C(n,3k) C(n+k,n) (3k)!/k!^3,
                      Cooper s18 (Beukers App. C 9.13 formula), Franel (x) central C(2n,n) sum_k C(n,k)^3.
usage: python a6c_more_families.py PLO PHI"""
import sys, os, time
from math import comb, factorial
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import primes_in
from a6_polys_modp import norm, pmod, deg, pgcd, pderiv, radical, sqfree_decomposition


def fam(N):
    F = [sum(comb(n, k) ** 3 for k in range(n + 1)) for n in range(N)]
    def s18(n):
        tot = 0
        for k in range(n // 2 + 1):
            a = comb(2 * n - 3 * k - 1, n) if 2 * n - 3 * k - 1 >= 0 else 0
            b = comb(2 * n - 3 * k, n) if 2 * n - 3 * k >= 0 else 0
            tot += (-1) ** k * comb(n, k) * comb(2 * k, k) * comb(2 * n - 2 * k, n - k) * (a + b)
        return tot
    def az(n):
        return sum((-1) ** k * 3 ** (n - 3 * k) * comb(n, 3 * k) * comb(n + k, n) * factorial(3 * k) // factorial(k) ** 3
                   for k in range(n // 3 + 1))
    return {
        ('w1', 'C(2k,k)^2'): [comb(2 * k, k) ** 2 for k in range(N)],
        ('w1', '(3k)!/k!^3'): [factorial(3 * k) // factorial(k) ** 3 for k in range(N)],
        ('w1', '(4k)!/((2k)!k!^2)'): [factorial(4 * k) // (factorial(2 * k) * factorial(k) ** 2) for k in range(N)],
        ('w1', '(6k)!/((3k)!(2k)!k!)'): [factorial(6 * k) // (factorial(3 * k) * factorial(2 * k) * factorial(k)) for k in range(N)],
        ('w2', 'C(2k,k)^3'): [comb(2 * k, k) ** 3 for k in range(N)],
        ('w2', 'C(2k,k)^2C(3k,k)'): [comb(2 * k, k) ** 2 * comb(3 * k, k) for k in range(N)],
        ('w2', 'C(2k,k)^2C(4k,2k)'): [comb(2 * k, k) ** 2 * comb(4 * k, 2 * k) for k in range(N)],
        ('w2', 'C(2k,k)C(3k,k)C(6k,3k)'): [comb(2 * k, k) * comb(3 * k, k) * comb(6 * k, 3 * k) for k in range(N)],
        ('w2', 'Almkvist-Zudilin'): [az(n) for n in range(N)],
        ('w2', 'Cooper s18'): [s18(n) for n in range(N)],
        ('w2', 'Franel(x)central'): [comb(2 * n, n) * F[n] for n in range(N)],
    }


if __name__ == '__main__':
    plo, phi = int(sys.argv[1]), int(sys.argv[2])
    primes = primes_in(plo, phi)
    t0 = time.time()
    fams = fam(2 * max(primes))
    print(f'sequences built in {time.time() - t0:.0f}s; heads:', {k[1]: v[:5] for k, v in fams.items()})
    for (w, name), u in fams.items():
        nl, badB, badC, simple = [], [], [], []
        for p in primes:
            if any((u[p + j] - u[1] * u[j]) % p for j in range(p)):
                nl.append(p)
                continue
            F = norm([u[k] % p for k in range(p)])
            V = norm([((u[p + j] - u[1] * u[j]) // p) % p for j in range(p)])
            if deg(F) <= 0:
                continue
            if w == 'w1':
                sqf = deg(pgcd(F, pderiv(F, p), p)) == 0
                cop = (not V) and False or (deg(pgcd(F, V, p)) == 0 if V else False)
                if not (sqf and cop):
                    badB.append((p, 'sqf' if not sqf else 'gcd'))
            else:
                R = radical(F, p)
                if V and pmod(V, R, p):
                    badC.append(p)
                if deg(sqfree_decomposition(F, p).get(1, [1])) > 0:
                    simple.append(p)
        if w == 'w1':
            print(f'{w} {name:24s}: Lucas fails at {nl}; (B) fails at {badB}', flush=True)
        else:
            print(f'{w} {name:24s}: Lucas fails at {nl}; rad(F_p)|V_p fails at {badC}; simple roots at {len(simple)} primes', flush=True)
