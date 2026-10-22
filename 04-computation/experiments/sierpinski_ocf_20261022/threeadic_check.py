# Brute-force checks of the 3-adic Sierpinski tournament S_k = Cay(Z/3^k, {x : unit part of x = 1 mod 3}):
#  (a) S_k is a tournament; (b) x -> 3x+b embeds S_(k-1) onto residue class b (self-similarity);
#  (c) class b => class b+1 (so S_k = C_3[S_(k-1)^3]); (d) multiplication by lambda is an automorphism if lambda = 1
#  mod 3 and an anti-automorphism (reverses every arc) if lambda = 2 mod 3, so halving reverses all arcs;
#  (e) along two integer Collatz trajectories with the same parity vector, the S_infinity arc between them flips once
#  per halving (checked on Z/3^k for random pairs).  mac-mini-2026-10-22-sierp.
import random
def unit_mod3(z, N):
    z %= N
    if z == 0: return None
    while z % 3 == 0: z //= 3
    return z % 3
def arc(x, y, N): return unit_mod3(y - x, N) == 1
ok = True
for k in range(1, 6):
    N = 3 ** k
    for x in range(N):
        for y in range(N):
            if x != y and arc(x, y, N) == arc(y, x, N): ok = False
    if k >= 2:
        M = N // 3
        for b in range(3):
            for x in range(M):
                for y in range(M):
                    if x != y and arc(x, y, M) != arc(3 * x + b, 3 * y + b, N): ok = False
        for x in range(N):
            for y in range(N):
                if (y - x) % 3 == 1 and not arc(x, y, N): ok = False
    for lam in range(1, N):
        if lam % 3 == 0: continue
        for _ in range(200):
            x, y = random.randrange(N), random.randrange(N)
            if x == y: continue
            same = arc(x, y, N) == arc(lam * x % N, lam * y % N, N)
            if same != (lam % 3 == 1): ok = False
print("checks (a)-(d) for k <= 5:", ok)
# (e) Collatz: C(x) = 3x+1 (odd) or x/2 (even); integers x, y with x = y mod 2^t share t-step parity vectors
ok2 = True; tested = 0
for _ in range(2000):
    t = random.randint(1, 18); k = 12; N = 3 ** k
    x = random.randrange(1, 10 ** 9); y = x + 2 ** t * random.randrange(1, 10 ** 6)
    a, b, halv = x, y, 0
    for _ in range(t):
        if a % 2: a, b = 3 * a + 1, 3 * b + 1
        else: a, b = a // 2, b // 2; halv += 1
    if (a % 2) != (b % 2): pass
    expected = arc(x, y, N) if halv % 2 == 0 else (not arc(x, y, N))
    if arc(a, b, N) != expected: ok2 = False
    tested += 1
print(f"check (e): arc flips once per halving along shared parity vectors ({tested} random pairs, Z/3^12):", ok2)
