"""Single primes: for every odd prime p <= B, counts generate (Z/p^2)^x (hence every residue mod every p^k is a count).
Generators: odd primes q <= 20405 other than p, and 7 (= 35/5 in the group).  For each prime l | p(p-1) we need a
generator outside the index-l subgroup: an l-th power non-residue mod p (l | p-1), or a non-Wieferich base (l = p).
Report the largest generator ever needed."""
import sys
def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** .5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]
def factor(n):
    f = []; d = 2
    while d * d <= n:
        if n % d == 0:
            f.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: f.append(n)
    return f
B = int(sys.argv[1])
gens = [q for q in primes_upto(20405) if q > 2]
worst = (0, None, None); fails = []
for p in primes_upto(B):
    if p == 2: continue
    for l in factor(p - 1) + [p]:
        need = None
        for q in gens:
            if q == p: continue
            if l == p:
                if pow(q, p - 1, p * p) != 1: need = q; break
            else:
                if pow(q, (p - 1) // l, p) != 1: need = q; break
        if need is None: fails.append((p, l))
        elif need > worst[0]: worst = (need, p, l)
print(f"odd primes p <= {B}: failures {fails[:5]} (count {len(fails)}); largest generator needed: q = {worst[0]} (at p = {worst[1]}, l = {worst[2]})")
