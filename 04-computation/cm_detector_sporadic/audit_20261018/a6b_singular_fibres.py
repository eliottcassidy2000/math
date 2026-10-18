"""Audit check 6b: per singular fibre, at which primes is it a simple root of F_p, and is that set the inert set of
an imaginary quadratic discriminant?  (HYP-9330 (C): 'simple roots exactly when p is inert in their CM field'.)"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from audit_seqs import seq, primes_in
from a6_polys_modp import norm, pmod, deg, sqfree_decomposition, kron

FACTORS = {  # irreducible-over-Q singular factors in x = 1/m
    'Apery': {'x^2-34x+1': [1, -34, 1]},
    'Domb': {'1-4x (m=4)': [1, -4], '1-16x (m=16)': [1, -16]},
    's7': {'1+x (m=-1)': [1, 1], '1-27x (m=27)': [1, -27]},
    's10': {'1+4x (m=-4)': [1, 4], '1-16x (m=16)': [1, -16]},
    'A274786': {'1-44x-16x^2': [1, -44, -16]},
}
DS = [-3, -4, -7, -8, -11, -12, -15, -16, -19, -20, -24, -27, -28, -32, -35, -36, -40, -48, -51, -52, -60, -64, -72,
      -75, -84, -88, -91, -99, -100, -112, -115, -120, -123, -147, -148]

if __name__ == '__main__':
    plo, phi = int(sys.argv[1]), int(sys.argv[2])
    primes = primes_in(plo, phi)
    for fam, facs in FACTORS.items():
        u = seq(fam, 2 * max(primes))
        simple = {p: sqfree_decomposition(norm([u[k] % p for k in range(p)]), p).get(1, [1]) for p in primes}
        for lab, f in facs.items():
            hit = []
            for p in primes:
                S = simple[p]
                fp = norm([c % p for c in f])
                if deg(fp) <= 0:
                    continue
                # does the factor (mod p) share a root with S?  i.e. gcd(S, f mod p) nontrivial
                a, b = S[:], fp[:]
                while b:
                    a, b = b, pmod(a, b, p)
                if deg(a) > 0:
                    hit.append(p)
            match = [D for D in DS if [p for p in primes if kron(D, p) == -1 and p not in ()] == hit]
            # allow primes dividing D: compare on primes not dividing D
            match2 = [D for D in DS if [p for p in primes if D % p and kron(D, p) == -1] == [p for p in hit if D % p]]
            print(f'{fam:8s} {lab:14s}: simple root at {len(hit)}/{len(primes)} primes; '
                  f'= inert set of D in {match}  (ignoring p|D: {match2})')
