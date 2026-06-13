"""COUNT LAW for signed-LRC collisions: characterize the exact silent-cut set & test the
squarefree deficiency conjecture 2^{(p+q)/2-2}.   monad-explorer-2026-06-06-S707.

For a prime q|C with subgroup K_q=<C/q>, flipping H_q (its half-system) is silent for cut eps iff
A_t := sum_{i not in K_q} eps_i sin(2 pi t i/C) = 0 for all t with q nmid t.
We:
  (1) enumerate the NON-K_q signings (much smaller than the full cube) and count solutions;
  (2) test whether the solution set = exactly {signings whose signed point set is a union of full
      K_q-cosets}  (the 'full-coset' characterization behind THM-417);
  (3) assemble deficiency = sum over primes q|C of (#silent pairs for H_q)  [valid for squarefree C,
      where classes are size 2 and each comes from a single H_q] and compare with brute force.
"""
import math
from itertools import product
from collections import defaultdict


def is_prime(C):
    return C > 1 and all(C % d for d in range(2, int(C**0.5) + 1))


def prime_factors(C):
    f, m, d = set(), C, 2
    while d * d <= m:
        while m % d == 0:
            f.add(d); m //= d
        d += 1
    if m > 1:
        f.add(m)
    return sorted(f)


def count_silent_for_q(C, q):
    """Return (#solutions on non-K signings, #that are full-coset unions, total non-K size)."""
    n1 = (C - 1) // 2
    Kq = set((C // q * j) % C for j in range(q))          # order-q subgroup
    nonK = [i for i in range(1, n1 + 1) if i not in Kq]    # runners not in K_q
    free_t = [t for t in range(1, n1 + 1) if t % q != 0]   # frequencies q nmid t
    m = C // q
    # precompute sin table
    sin = {(t, i): math.sin(2 * math.pi * t * i / C) for t in free_t for i in nonK}
    sols = 0
    fullcoset = 0
    for bits in product((1, -1), repeat=len(nonK)):
        ok = True
        for t in free_t:
            A = 0.0
            for idx, i in enumerate(nonK):
                A += bits[idx] * sin[(t, i)]
            if abs(A) > 1e-7:
                ok = False
                break
        if ok:
            sols += 1
            # full-coset test: signed points, grouped by residue mod m, each residue class either
            # entirely present as a full coset
            pts = set((bits[idx] * nonK[idx]) % C for idx in range(len(nonK)))
            # for each canonical residue r in 1..(m-1)/2, the union of pts must be a full coset r+K or (m-r)+K
            is_fc = True
            byres = defaultdict(set)
            for x in pts:
                byres[x % m].add(x)
            for r, s in byres.items():
                full = set((r + (C // q) * j) % C for j in range(q))
                if s != full:
                    is_fc = False
                    break
            if is_fc:
                fullcoset += 1
    return sols, fullcoset, len(nonK)


print("Characterizing silent cuts per prime q|C, and assembling the deficiency.")
print(f"{'C':>4} {'fact':>10} {'q':>4} {'nonK':>5} {'#sol(nonK)':>11} {'#fullcoset':>11} "
      f"{'pairs(Hq)':>10}")
defic_pred = {}
for C in range(9, 40, 2):
    if is_prime(C):
        continue
    pf = prime_factors(C)
    fstr = "x".join(map(str, []))
    # factorization string
    f2, m2, d = [], C, 2
    while d * d <= m2:
        while m2 % d == 0:
            f2.append(d); m2 //= d
        d += 1
    if m2 > 1:
        f2.append(m2)
    fstr = "x".join(map(str, f2))
    total_pairs = 0
    for q in pf:
        if (C - 1) // 2 - (q - 1) // 2 > 17:
            print(f"{C:>4} {fstr:>10} {q:>4} {(C-1)//2-(q-1)//2:>5} {'(skip>2^17)':>11}")
            total_pairs = None
            continue
        sols, fc, nk = count_silent_for_q(C, q)
        # silent cuts for H_q = sols(nonK) * 2^{|H_q|}(free signs on H_q) ; pairs = /2 ; mod global flip /...
        Hq = (q - 1) // 2
        cuts = sols * (2 ** Hq)        # full sign vectors in {+-1}^{n-1}, before global-flip quotient
        pairs = cuts // 2 // 2         # /2 pairing (eps, eps+Hq); /2 global-flip quotient
        if total_pairs is not None:
            total_pairs += pairs
        print(f"{C:>4} {fstr:>10} {q:>4} {nk:>5} {sols:>11} {fc:>11} {pairs:>10}")
    defic_pred[C] = total_pairs
    if total_pairs is not None:
        # squarefree two-prime conjecture
        conj = ""
        if len(pf) == 2 and all(C % (p*p) for p in pf):
            p_, q_ = pf
            conj = f"  2^((p+q)/2-2)=2^{(p_+q_)//2-2}={2**((p_+q_)//2-2)}"
        print(f"{C:>4} {fstr:>10}  ==> predicted deficiency = {total_pairs}{conj}")
print("\n(compare 'predicted deficiency' with brute-force deficiency from s707/s705 tables)")
