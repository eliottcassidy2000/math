"""clock_arrows_circulant_20261026.py -- opus-2026-10-26-S62

The owner's clock-arrow reading of an integer sequence (a_k): its PATTERN modulo m is the set S_m = {a_k mod m} of
residues attained (0 removed), and the arrows on the clock Z/m are i -> i + s for s in S_m: the circulant digraph
C(m, S_m).  It is an oriented graph (a "partial tournament": every pair joined by at most one arc) iff S_m and -S_m are
disjoint, and a tournament iff moreover S_m u -S_m = Z/m \ {0}, i.e. S_m is a skew half-set (m odd).  Squares modulo a
prime p = 3 mod 4 give the Paley tournament; powers of 2 modulo p give the subgroup <2>, a tournament iff -1 is not in
<2> and <2> has index 2 (THM-126, THM-4644, THM-4529 for the prior circulant results).
This script tabulates, for classical sequences and moduli m <= 120, which (sequence, m) give oriented graphs and
which give tournaments, and the arc density |S_m|/(m-1) otherwise.
Usage: python3 clock_arrows_circulant_20261026.py
"""
from sympy import isprime, primerange, factorint

def seq_squares(L): return [k * k for k in range(1, L)]
def seq_cubes(L): return [k ** 3 for k in range(1, L)]
def seq_pow2(L): return [2 ** k for k in range(0, L)]
def seq_pow3(L): return [3 ** k for k in range(0, L)]
def seq_fib(L):
    a, b, out = 1, 1, []
    for _ in range(L): out.append(a); a, b = b, a + b
    return out
def seq_primes(L): return list(primerange(2, L))
def seq_tri(L): return [k * (k + 1) // 2 for k in range(1, L)]
def seq_fact(L):
    f, out = 1, []
    for k in range(1, L): f *= k; out.append(f)
    return out
def seq_collatz_odd(L):                       # the odd Collatz successors (3k+1)/2^v of odd k
    out = []
    for k in range(1, L, 2):
        x = 3 * k + 1
        while x % 2 == 0: x //= 2
        out.append(x)
    return out

SEQS = {"squares": seq_squares, "cubes": seq_cubes, "2^k": seq_pow2, "3^k": seq_pow3, "Fibonacci": seq_fib,
        "primes": seq_primes, "triangular": seq_tri, "factorials": seq_fact, "Collatz odd successors": seq_collatz_odd}

def pattern(seq, m):
    S = set(a % m for a in seq) - {0}
    return S

def main():
    out = []
    for name, fn in SEQS.items():
        seq = fn(3000)
        oriented, tourn = [], []
        for m in range(3, 121):
            S = pattern(seq, m)
            if not S: continue
            neg = {(-s) % m for s in S}
            if S & neg: continue
            oriented.append(m)
            if len(S) == (m - 1) // 2 and m % 2 == 1: tourn.append(m)
        out.append(f"{name}: oriented moduli (S and -S disjoint) {oriented}; tournament moduli (skew half-set) {tourn}")
    # squares: tournaments exactly at primes p = 3 mod 4 (Paley) -- check and also prime powers
    sq = seq_squares(3000); sqt = [m for m in range(3, 121, 2) if len(pattern(sq, m)) == (m - 1) // 2 and not (pattern(sq, m) & {(-s) % m for s in pattern(sq, m)})]
    out.append(f"squares: tournament moduli {sqt} versus primes = 3 mod 4 below 121: {[p for p in primerange(3, 121) if p % 4 == 3]}")
    # powers of 2: tournament moduli versus the criterion ord_p(2) = (p-1)/2 and p = 7 mod 8
    p2 = seq_pow2(3000)
    crit = []
    for p in primerange(3, 121):
        o = 1; x = 2 % p
        while x != 1: x = x * 2 % p; o += 1
        if o == (p - 1) // 2 and p % 8 == 7: crit.append(p)
    p2t = [m for m in range(3, 121, 2) if len(pattern(p2, m)) == (m - 1) // 2 and not (pattern(p2, m) & {(-s) % m for s in pattern(p2, m)})]
    out.append(f"2^k: tournament moduli {p2t} versus primes with ord_p(2) = (p-1)/2 and p = 7 mod 8: {crit}")
    print("\n".join(out))

if __name__ == "__main__":
    main()
