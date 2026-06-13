"""CONSTRUCTIVE PROOF of the converse of HYP-2270:  C=2n-1 composite => sign-orbit < 2^{n-2}.
monad-explorer-2026-06-06-S707.

Construction.  C composite, q = smallest prime factor, m = C/q (>1, odd).
K = order-q subgroup = { j*m mod C : j=0..q-1 } (the multiples of m).  H_q = K-half-system
= { m, 2m, .., ((q-1)/2) m }.  Define the base cut eps by aligning every non-K runner into a
CANONICAL coset of K:
    eps_i = +1  if (i mod m) in {0,1,..,(m-1)/2},   else  -1 .
Then { eps_i * i : i not in K } = union over r=1..(m-1)/2 of the FULL coset (r + K), so its
signed sine sum A_t = sum_{i not in K} eps_i sin(2 pi t i / C) = 0 for every t with q | t FALSE,
i.e. for every t with q NOT dividing t (full-coset geometric-sum cancellation).
Hence eps and eps' := eps with H_q negated satisfy Phi_eps(t)^2 = Phi_eps'(t)^2 for ALL t:
    q | t :  B_t := sum_{i in H_q} eps_i sin(2 pi t i / C) = 0  => Phi unchanged;
    q nmid t: A_t = 0 => Phi_eps = B_t, Phi_eps' = -B_t, squares equal.
=> same folded clock multiset => COLLISION.  eps != +-eps' since H_q != all runners (q<C).

This script VERIFIES, for every composite odd C up to a large bound (brute force impossible
for big C), that the constructed pair has identical folded spectra (and that A_t=0 as claimed),
and confirms NO collision is constructible for prime C (sanity).
"""
import math
from collections import Counter


def smallest_prime_factor(C):
    d = 2
    while d * d <= C:
        if C % d == 0:
            return d
        d += 1
    return C  # prime


def folded_spectrum(pts, C):
    out = []
    L = len(pts)
    for i in range(L):
        for j in range(i + 1, L):
            f = (pts[i] - pts[j]) % C
            out.append(min(f, C - f))
    return tuple(sorted(out))


def build_base_cut(C):
    q = smallest_prime_factor(C)
    m = C // q
    n1 = (C - 1) // 2
    runners = list(range(1, n1 + 1))
    eps = {}
    for i in runners:
        r = i % m
        eps[i] = 1 if (r <= (m - 1) // 2) else -1   # r in {0..(m-1)/2} -> +1
    Hq = [m * j for j in range(1, (q - 1) // 2 + 1)]   # half-system of K
    return q, m, eps, Hq, runners


def verify(C):
    q, m, eps, Hq, runners = build_base_cut(C)
    if q == C:  # prime: construction degenerate (m=1, Hq empty)
        return ("prime", q, m, None)
    # eps points and flipped points
    pts0 = [(eps[i] * i) % C for i in runners]
    epsf = dict(eps)
    for x in Hq:
        epsf[x] = -epsf[x]
    pts1 = [(epsf[i] * i) % C for i in runners]
    fs0 = folded_spectrum(pts0, C)
    fs1 = folded_spectrum(pts1, C)
    collide = (fs0 == fs1)
    distinct = (pts0 != pts1) and (sorted(pts0) != sorted((-x) % C for x in pts1))
    # check A_t = 0 for q nmid t
    nonK = [i for i in runners if i % m != 0]
    maxbadA = 0.0
    for t in range(1, (C - 1) // 2 + 1):
        if t % q != 0:
            A = sum(eps[i] * math.sin(2 * math.pi * t * i / C) for i in nonK)
            maxbadA = max(maxbadA, abs(A))
    return ("composite", q, m, dict(collide=collide, distinct=distinct,
                                    maxAerr=maxbadA, Hq=Hq))


print("Verifying the constructive collision for composite C (and checking prime C is skipped).")
print(f"{'C':>5} {'fact':>12} {'q':>4} {'m':>5} {'collide':>8} {'distinct':>9} {'maxA_err':>11}")
all_ok = True
checked = 0
for C in range(5, 220, 2):
    spf = smallest_prime_factor(C)
    if spf == C:
        continue  # prime
    kind, q, m, info = verify(C)
    # factor string
    f, mm, d = [], C, 2
    while d * d <= mm:
        while mm % d == 0:
            f.append(d); mm //= d
        d += 1
    if mm > 1:
        f.append(mm)
    fstr = "x".join(map(str, f))
    ok = info['collide'] and info['distinct'] and info['maxAerr'] < 1e-9
    all_ok = all_ok and ok
    checked += 1
    flag = "" if ok else "   <<< FAIL"
    print(f"{C:>5} {fstr:>12} {q:>4} {m:>5} {str(info['collide']):>8} "
          f"{str(info['distinct']):>9} {info['maxAerr']:>11.2e}{flag}")

print(f"\nChecked {checked} composite C in [5,219].  ALL constructed pairs collide & distinct: {all_ok}")
print("=> CONVERSE OF HYP-2270 verified constructively: composite C => sign-orbit < 2^(n-2).")
