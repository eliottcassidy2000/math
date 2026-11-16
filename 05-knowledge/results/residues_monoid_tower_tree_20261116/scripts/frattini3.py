#!/usr/bin/env python3
"""frattini3.py: generation of (Z/m)^x for every modulus m with prime factors in P (any exponents) by realised counts
coprime to m, via the Frattini criterion (components: p = 2: bits x = 3 mod 4 and x in {3,5} mod 8; q | p-1: the q-th
power residue class mod p; q = p: (x^(p-1)-1)/p mod p).  Cases: (A') P = {2,3,5}, generators = first 32 odd numbers
minus 7, 21, coprime to 30; (B') P = primes <= 63 except 7, generators = odd numbers <= 4000 (not 7, 21) whose prime
factors are 7 or >= 67; (B'') P = all primes <= 63, generators with all prime factors >= 67.  All generators are
Hamiltonian-path counts (every odd number <= 2^32 except 7, 21: THM-4875 (D); below 609 already at order <= 8, THM-1370)."""
def primes_upto(n): return [p for p in range(2, n+1) if all(p % d for d in range(2, int(p**0.5)+1))]
P_all = primes_upto(63)
def factors_ok(x, allowed_small):
    y = x
    for p in primes_upto(63):
        while y % p == 0:
            if p not in allowed_small: return False
            y //= p
    return True
def dlog_q(y, p, q):
    g = next(pow(c, (p-1)//q, p) for c in range(2, p) if pow(c, (p-1)//q, p) != 1)
    return next(k for k in range(q) if pow(g, k, p) == y)
def comp(x, q, P):
    v = []
    for p in P:
        if p == 2:
            if q == 2: v += [1 if x % 4 == 3 else 0, 1 if x % 8 in (3, 5) else 0]
        else:
            if (p-1) % q == 0: v.append(dlog_q(pow(x, (p-1)//q, p), p, q))
            if q == p: v.append(((pow(x, p-1, p*p) - 1)//p) % p)
    return v
def rank_mod_q(rows, q):
    rows = [r[:] for r in rows]; rank = 0
    for col in range(len(rows[0]) if rows else 0):
        piv = next((i for i in range(rank, len(rows)) if rows[i][col] % q), None)
        if piv is None: continue
        rows[rank], rows[piv] = rows[piv], rows[rank]; inv = pow(rows[rank][col], -1, q)
        rows[rank] = [(v*inv) % q for v in rows[rank]]
        for i in range(len(rows)):
            if i != rank and rows[i][col] % q:
                f = rows[i][col]; rows[i] = [(a - f*b) % q for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank
cases = [
  ("(A') P = {2,3,5}; generators = first 32 odd numbers minus 7, 21, coprime to 30", [2,3,5], [x for x in range(3, 64, 2) if x not in (7, 21) and x % 3 and x % 5]),
  ("(B') P = primes <= 63 except 7; generators = odd <= 4000 with prime factors 7 or >= 67", [p for p in P_all if p != 7], [x for x in range(3, 4001, 2) if x not in (7, 21) and factors_ok(x, {7})]),
  ("(B'') P = all primes <= 63; generators = odd <= 4000 with all prime factors >= 67", P_all, [x for x in range(3, 4001, 2) if factors_ok(x, set())]),
]
for name, P, G in cases:
    qs = sorted(set(q for p in P for q in primes_upto(63) if (p == 2 and q == 2) or (p > 2 and ((p-1) % q == 0 or q == p))))
    ok = True; detail = []
    for q in qs:
        vecs = [comp(x, q, P) for x in G]; need = len(vecs[0]); r = rank_mod_q(vecs, q)
        detail.append(f"q={q}:{r}/{need}"); ok = ok and r == need
    print(name, f"({len(G)} generators, e.g. {G[:8]}):", ' '.join(detail), "->", "GENERATES (Z/m)^x for all exponents" if ok else "FAILS")
def closure(gens, m, odd_only):
    seen = {1 % m}; frontier = [1 % m]
    while frontier:
        nxt = []
        for x in frontier:
            for g in gens:
                y = (x*g) % m
                if y not in seen: seen.add(y); nxt.append(y)
        frontier = nxt
    return [r for r in range(m) if (r % 2 == 1 or not odd_only) and r not in seen]
S = [x for x in range(3, 64, 2) if x not in (7, 21)]
bad = [(a, b, c) for a in range(0, 8) for b in range(0, 4) for c in range(0, 3) if 1 < 2**a*3**b*5**c <= 200000 and closure(S, 2**a*3**b*5**c, a >= 1)]
print("(A) replay: products of the first 32 odd numbers minus 7, 21 miss an odd residue for moduli 2^a 3^b 5^c <= 2e5 (a<8,b<4,c<3):", bad if bad else "none")
L3 = [1 + 3*2**s for s in range(4, 40) if s % 3 == 1]
G7 = S + [x for x in range(65, 2000, 2) if factors_ok(x, {7})] + [49, 343] + L3
bad2 = [(m, closure(G7, m, m % 2 == 0)[:4]) for m in [2**3*3*5*7, 2**3*5*7*11, 3*5*7*11*13, 2**4*7*13, 49*3*5*8, 2**3*9*5*7*11, 343*2, 49*11*13, 7*11*13*17] if closure(G7, m, m % 2 == 0)]
print("(C) replay with 7 in the modulus (S + units with factors 7 or >= 67 up to 2000 + 49, 343 + L_3 multiples of 7):", bad2 if bad2 else "none missed")
Pn7 = [p for p in P_all if p not in (2, 7)]
found = {s0: next(s for s in range(4, 200000) if s % 3 == 1 and s % 7 == s0 and (1 + 3*2**s) % 49 != 0 and all((1 + 3*2**s) % p for p in Pn7)) for s0 in range(7)}
print("(C') for each s mod 7 the least s = 1 mod 3 with 1 + 3*2^s a unit modulo every prime <= 63 other than 7 and not divisible by 49:", found)
print("(C'') u' = (1 + 3*2^s)/7 mod 7^2 as s runs over s = 4 mod 21 avoided? classes mod 49 hit by s in [4, 4+21*49) with s = 1 mod 3:", len(set(((1+3*2**s)//7) % 49 for s in range(4, 4+21*49) if s % 3 == 1 and (1+3*2**s) % 7 == 0)), "of 49")
