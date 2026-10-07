"""Autocorrelation of the flip (volatility) process sigma_n = e_n mod 2 of the lag-1 pair chain,
among paths alive on the whole window [W0, W1).  Also the cross moment E[dk_m dk_(m+j)] (should be 0)
and the 'leverage' correlation E[dk_m sigma_(m+j)].  usage: N W0 W1 seed"""
import random, sys, math
from pair_chain_verify import step

def main(N, W0, W1, seed, J=(1, 2, 3, 4, 8, 16, 32, 64)):
    rng = random.Random(seed)
    paths = 0
    sums = {j: [0.0, 0.0, 0.0] for j in J}   # sigma acf numerators, dk dk, dk sigma
    cnts = {j: 0 for j in J}
    s1 = 0; s2 = 0; cnt = 0
    while paths < N:
        k, num = 0, 1; sig = []; dks = []; alive = True
        for n in range(W1):
            beta = 0 if n < 4 else rng.getrandbits(1)
            if n >= W0:
                sig.append(num & 1)
            k2, num2, _ = step(k, num, beta)
            if n >= W0:
                dks.append(k2 - k)
            k, num = k2, num2
            if (k, num) == (0, 0):
                alive = False; break
        if not alive:
            continue
        paths += 1
        L = len(sig)
        for x in sig:
            s1 += x; s2 += x * x; cnt += 1
        for j in J:
            for m in range(L - j):
                sums[j][0] += sig[m] * sig[m + j]
                sums[j][1] += dks[m] * dks[m + j]
                sums[j][2] += dks[m] * sig[m + j]
            cnts[j] += L - j
    mu = s1 / cnt; var = s2 / cnt - mu * mu
    print(f"paths alive on [{W0},{W1}): {N}; flip density {mu:.4f}, var {var:.4f}")
    for j in J:
        c = sums[j][0] / cnts[j] - mu * mu
        print(f"  lag {j:3d}: corr(sigma_m, sigma_m+j) = {c/var:+.4f}   E[dk_m dk_m+j] = {sums[j][1]/cnts[j]:+.4f}   E[dk_m sigma_m+j] = {sums[j][2]/cnts[j]:+.4f}")
    print(f"  (s.e. of a correlation ~ {1/math.sqrt(cnts[J[0]]):.4f} if independent)")

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
