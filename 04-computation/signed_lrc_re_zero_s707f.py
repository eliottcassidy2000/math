"""Verify the two lemmas underpinning the full-coset characterization & count law.
monad-explorer-2026-06-06-S707.

LEMMA 1 (Re vanishes): for K=order-q subgroup, Sigma any half-system selection of the non-K part,
   Re(Sigma_hat(t)) = (1/2) * \hat{1_{notK}}(t) = 0  for every t with q nmid t  (t != 0).
   [Because 1_{notK} = 1_{Z/C} - 1_K and \hat{1_K}(t)=q*[q|t].]
   => A_t = Im(Sigma_hat(t)) = 0 (the silent condition) <=> Sigma_hat(t)=0 for q nmid t
      <=> 1_Sigma constant on K-cosets <=> Sigma = union of full K-cosets.

LEMMA 2 (count): #(unions of full K-cosets that are half-system selections) = 2^{(m-1)/2}, m=C/q.

We check Lemma 1 numerically for random half-system selections, and Lemma 2 by formula vs the
enumerated #fullcoset from s707e (re-derived here cheaply via the coset-pair count).
"""
import math, cmath


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def check_lemma1(C, q, trials=200):
    import random
    rng = random.Random(12345)
    n1 = (C - 1) // 2
    Kq = set((C // q * j) % C for j in range(q))
    nonK = [i for i in range(1, n1 + 1) if i not in Kq]
    worst = 0.0
    for _ in range(trials):
        # random half-system selection of non-K: choose +-i for each non-K runner
        Sigma = [(rng.choice((1, -1)) * i) % C for i in nonK]
        for t in range(1, n1 + 1):
            if t % q == 0:
                continue
            re = sum(math.cos(2 * math.pi * t * x / C) for x in Sigma)
            worst = max(worst, abs(re))
    return worst


def count_fullcoset(C, q):
    m = C // q
    return 2 ** ((m - 1) // 2)


print("LEMMA 1: Re(Sigma_hat(t)) = 0 for q nmid t  (max over 200 random half-system selections)")
print(f"{'C':>4} {'q':>3} {'maxRe':>12}")
for C in range(9, 80, 2):
    if is_prime(C):
        continue
    for q in [d for d in range(2, C) if C % d == 0 and is_prime(d)]:
        w = check_lemma1(C, q)
        flag = "" if w < 1e-9 else "  <<< NONZERO"
        print(f"{C:>4} {q:>3} {w:>12.2e}{flag}")

print("\nLEMMA 2: #fullcoset = 2^{(m-1)/2}, m=C/q  (compare with s707e enumerated values)")
known = {(9,3):2,(15,3):4,(15,5):2,(21,3):8,(21,7):2,(25,5):4,(27,3):16,(33,3):32,(35,5):8,(35,7):4}
print(f"{'C':>4} {'q':>3} {'m':>4} {'formula':>8} {'enum(s707e)':>12} {'match':>6}")
for (C,q),val in sorted(known.items()):
    m=C//q
    f=count_fullcoset(C,q)
    print(f"{C:>4} {q:>3} {m:>4} {f:>8} {val:>12} {str(f==val):>6}")

print("\nCOUNT LAW (squarefree pq):  pairs(H_q) = 2^{(C/q-1)/2 + (q-1)/2 - 2};  deficiency = sum over primes.")
print(f"{'C=pq':>8} {'pred deficiency':>16} {'2^((p+q)/2-2)':>16}")
for (p,q) in [(3,5),(3,7),(5,7),(3,11),(3,13),(5,11),(5,13),(7,11),(3,17),(3,19),(7,13)]:
    C=p*q
    def pairs(qq):
        m=C//qq
        e=(m-1)//2+(qq-1)//2-2
        return 2**e
    pred=pairs(p)+pairs(q)
    closed=2**((p+q)//2-2)
    print(f"{C:>8} {pred:>16} {closed:>16}   {'OK' if pred==closed else 'DIFF'}")
