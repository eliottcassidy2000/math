"""(1) survival from several initial states; (2) return-map drift E|e_next|^theta vs |e|^theta from (0,E)."""
import random, math, sys
from pair_chain_verify import step

def survival(k0, num0, N, cps, rng):
    alive = {c: 0 for c in cps}; nmax = max(cps)
    for _ in range(N):
        k, num = k0, num0; ci = 0
        for n in range(nmax):
            k, num, _ = step(k, num, rng.getrandbits(1))
            if k == 0 and num == 0: break
            while ci < len(cps) and cps[ci] == n + 1:
                alive[cps[ci]] += 1; ci += 1
    return {c: alive[c] / N for c in cps}

def next_return(E, rng, theta=0.5, maxsteps=10**6):
    k, num = 0, E
    left = False
    for n in range(maxsteps):
        k2, num2, _ = step(k, num, rng.getrandbits(1))
        if k != 0 and k2 == 0:
            return abs(num2) ** theta
        k, num = k2, num2
    return None

rng = random.Random(99)
cps = [100, 400, 1600, 6400]
for (k0, num0, lab) in [(1, 0, "u=3v"), (0, 2**20, "u=v+2^20"), (0, 1, "u=v+1"), (5, 7, "(5,7)"), (-3, 1, "u=v/27+1/27"), (9, 3**9 - 1, "lag-9 Mersenne")]:
    S = survival(k0, num0, 600, cps, rng)
    print(f"{lab:16s} S(n): " + "  ".join(f"{c}:{S[c]:.3f}(sqrt n S={math.sqrt(c)*S[c]:.1f})" for c in cps))
theta = 0.5
s = 2 ** theta * (1 - math.sqrt(1 - 0.75 ** theta)); lam = 2 ** -theta * s
print(f"theta={theta}: s={s:.4f}, lambda=2^-theta s={lam:.4f}")
for E in [10**3, 10**5, 10**7, 10**9, 10**12, 10**15, 3**30, 2**40]:
    vals = [next_return(E, rng, theta) for _ in range(3000)]
    vals = [v for v in vals if v is not None]
    m = sum(vals) / len(vals)
    print(f"  E={E:.3e}: E|e_next|^1/2 = {m:12.2f}   |E|^1/2 = {E**0.5:12.2f}   ratio = {m/E**0.5:.4f}   (returns {len(vals)})")
