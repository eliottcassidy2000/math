# taxicab_split_lemma_kpc1_verify.py
# ADVERSARIAL VERIFIER (kind-pasteur-2026-06-10-S1) for THM-461 claim A3 (split lemma).
# Checks, for ALL primitive pairs 1 <= y <= x <= 2000 (gcd(x,y)=1), with q = x^2-xy+y^2:
#   (0) q odd
#   (1) every prime p != 3 dividing q satisfies p == 1 (mod 3)
#   (2) v3(q) <= 1, with v3(q) = 1  <=>  3 | (x+y)
# Method: smallest-prime-factor sieve (different from any direct trial division route),
# pair count cross-checked against an independent Euler-phi sieve.
import math

X = 2000
M = X * X  # q = x^2 - xy + y^2 <= x^2 for 1 <= y <= x
print(f"sieving spf up to {M} ...")
spf = list(range(M + 1))
for i in range(2, math.isqrt(M) + 1):
    if spf[i] == i:
        for j in range(i * i, M + 1, i):
            if spf[j] == j:
                spf[j] = i

pairs = 0
viol_odd = viol_inert = viol_v3 = viol_v3iff = 0
g = math.gcd
for x in range(1, X + 1):
    for y in range(1, x + 1):
        if g(x, y) != 1:
            continue
        pairs += 1
        q = x * x - x * y + y * y
        if q % 2 == 0:
            viol_odd += 1
        v3 = 0
        qq = q
        while qq % 3 == 0:
            qq //= 3
            v3 += 1
        if v3 > 1:
            viol_v3 += 1
        if (v3 == 1) != ((x + y) % 3 == 0):
            viol_v3iff += 1
        t = qq
        while t > 1:
            p = spf[t]
            if p % 3 != 1:
                viol_inert += 1
                print(f"  INERT VIOLATION: x={x} y={y} q={q} prime={p}")
                break
            while t % p == 0:
                t //= p

print(f"primitive pairs 1<=y<=x<={X}: {pairs}")
# independent count via totient sieve
phi = list(range(X + 1))
for i in range(2, X + 1):
    if phi[i] == i:  # i prime
        for j in range(i, X + 1, i):
            phi[j] -= phi[j] // i
print(f"sum_(x<=X) phi(x) (independent pair count): {sum(phi[1:])}")
print(f"violations: q even: {viol_odd}, inert/3-adic prime p!=3 with p%3!=1: {viol_inert}, "
      f"v3(q)>1: {viol_v3}, (v3=1 <=> 3|x+y) failures: {viol_v3iff}")
print("DONE taxicab_split_lemma_kpc1_verify")
