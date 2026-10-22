# Exact linear-forest numbers N_L of every doubly regular 23-tournament class (drt23_classes.json): the forward-step
# polynomial F(T, x) is evaluated mod three primes at x = 0..22 (forest_numbers.c), interpolated and lifted by CRT
# (0 <= f_k <= 23!), and (n - L)! N_L = sum_k f_k [t^L] (1 + t)^k (1 - t)^(n-1-k) (THM-4757 section 5).
# Checks: f_0 = f_(n-1) = H (hk64), sum_k f_k = n!, f palindromic, N_L = 0 for odd L, sum_L N_L/(n)_L = H 2^(n-1)/n!.
# Output: drt23_forest_exact.json (f_k and N_L as decimal strings), drt23_forest_exact.out.  mac-mini-2026-10-22-loop4.
import json, subprocess, sys
from math import factorial, comb
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from sympy import isprime

PRIMES = [2147483647, 2147483629, 2147483587]
assert all(isprime(p) for p in PRIMES)
cls = json.load(open('drt23_classes.json'))
n = 23

def run(job):
    k, p = job
    r = subprocess.run(['./forest_numbers', str(p)], input=cls[k]['code'] + '\n', capture_output=True, text=True, check=True)
    vals = {}
    for line in r.stdout.split('\n'):
        if line.strip(): x, v = map(int, line.split()); vals[x] = v
    return k, p, vals

def interpolate_mod(ys, p):
    """coefficients c_0..c_(m-1) mod p of the polynomial with values ys[x] at x = 0..m-1 (Newton form)"""
    m = len(ys)
    dd = list(ys)
    for j in range(1, m):
        for i in range(m - 1, j - 1, -1):
            dd[i] = (dd[i] - dd[i - 1]) * pow(j, -1, p) % p      # nodes are 0..m-1: x_i - x_(i-j) = j
    coef = [0] * m                                                 # expand sum dd[j] prod_(i<j) (x - i)
    basis = [1]
    for j in range(m):
        for i, b in enumerate(basis): coef[i] = (coef[i] + dd[j] * b) % p
        nb = [0] * (len(basis) + 1)
        for i, b in enumerate(basis):
            nb[i + 1] = (nb[i + 1] + b) % p; nb[i] = (nb[i] - j * b) % p
        basis = nb
    return coef

def crt(rs, ps):
    x, M = 0, 1
    for r, p in zip(rs, ps):
        t = (r - x) * pow(M, -1, p) % p
        x += M * t; M *= p
    return x, M

if __name__ == "__main__":
    jobs = [(k, p) for k in range(len(cls)) for p in PRIMES]
    with ThreadPoolExecutor(max_workers=2) as ex:
        results = list(ex.map(run, jobs))
    vals = {(k, p): v for k, p, v in results}
    out = []
    for c in cls:
        k = c['class']
        co = [interpolate_mod([vals[(k, p)][x] for x in range(n)], p) for p in PRIMES]
        f = []
        for i in range(n):
            x, M = crt([co[j][i] for j in range(3)], PRIMES)
            assert x <= factorial(n) < M
            f.append(x)
        assert f[0] == c['H'] and f[n - 1] == c['H'], 'F(T,0) and top coefficient must equal H'
        assert sum(f) == factorial(n) and f == f[::-1]
        N = []
        for L in range(n):
            s = 0
            for kk, fk in enumerate(f):
                # [t^L] (1+t)^kk (1-t)^(n-1-kk)
                s += fk * sum(comb(kk, a) * comb(n - 1 - kk, L - a) * (-1) ** (L - a) for a in range(max(0, L - (n - 1 - kk)), min(kk, L) + 1))
            assert s % factorial(n - L) == 0
            N.append(s // factorial(n - L))
        assert all(N[L] == 0 for L in range(1, n, 2)) and N[0] == 1
        kap = sum(Fraction(N[L], factorial(n) // factorial(n - L)) for L in range(n))
        assert kap == Fraction(c['H'] * 2 ** (n - 1), factorial(n))
        print(f"{c['label']:7s}: F(T,x) exact (f_0 = H = {f[0]}, sum = 23!, palindromic); N_L for even L: "
              + ' '.join(f"{N[L]}" for L in range(0, 9, 2)) + ' ...; forest identity sum N_L/(n)_L = H 2^22/23! exact', flush=True)
        out.append({'class': k, 'label': c['label'], 'H': str(c['H']), 'f': [str(x) for x in f], 'N': [str(x) for x in N]})
    json.dump(out, open('drt23_forest_exact.json', 'w'), indent=1)
