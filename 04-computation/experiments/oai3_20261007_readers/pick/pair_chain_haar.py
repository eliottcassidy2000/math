"""Haar-model Monte Carlo of the pair chain u = 3^k v + e (Terras T-steps).
Lag-1 Mersenne switch = pair (p, p-1), p Haar on 1 + 16 Z_2  ->  start (k,e) = (0,1),
first four driving bits 0 (v = p-1 = 16 m), then i.i.d. fair bits.
Lag-D start: (D, 3^D - 1) with v = q_D Haar on its coset (not used by default).
Outputs: survival S(n), sqrt(n) S(n), flip density among survivors, E k_n, E k_n^2,
|f| = |e|/3^|k| at returns of k to 0.
usage: python3 pair_chain_haar.py N NMAX seed
"""
import random, sys, math
from pair_chain_verify import step

def run(N=2000, NMAX=20000, seed=11, k0=0, num0=1, forced_zero=4):
    rng = random.Random(seed)
    checkpoints = [27, 50, 100, 200, 400, 800, 1600, 3200, 6400, 12800, 19000, 25600, 40000]
    checkpoints = [c for c in checkpoints if c <= NMAX]
    alive_at = {c: 0 for c in checkpoints}
    # flip density among survivors in windows
    windows = [(0, 100), (100, 1000), (1000, 5000), (5000, 20000), (20000, 40000)]
    wflip = {w: [0, 0] for w in windows}
    k_sum = {c: 0 for c in checkpoints}; k2_sum = {c: 0 for c in checkpoints}
    ret_f = []          # |f| at returns to k = 0 (survivors), after n >= 200
    merge_times = []
    for s in range(N):
        k, num = k0, num0
        merged_at = None
        bits = rng.getrandbits(64); nb = 64
        ci = 0
        for n in range(NMAX):
            if n < forced_zero:
                beta = 0
            else:
                if nb == 0:
                    bits = rng.getrandbits(64); nb = 64
                beta = bits & 1; bits >>= 1; nb -= 1
            flip = num & 1
            for w in windows:
                if w[0] <= n < w[1]:
                    wflip[w][0] += flip; wflip[w][1] += 1
            k_prev = k
            k, num, _ = step(k, num, beta)
            if k == 0 and num == 0:
                merged_at = n + 1
                break
            if k == 0 and k_prev != 0 and n >= 200:
                ret_f.append(abs(num))
            while ci < len(checkpoints) and checkpoints[ci] == n + 1:
                alive_at[checkpoints[ci]] += 1
                k_sum[checkpoints[ci]] += k; k2_sum[checkpoints[ci]] += k * k
                ci += 1
        if merged_at is not None:
            merge_times.append(merged_at)
    return checkpoints, alive_at, k_sum, k2_sum, wflip, ret_f, merge_times

if __name__ == "__main__":
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    NMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 20000
    seed = int(sys.argv[3]) if len(sys.argv) > 3 else 11
    cps, alive, ks, k2s, wflip, ret_f, mt = run(N, NMAX, seed)
    print(f"N={N} NMAX={NMAX} seed={seed}; start (k,e)=(0,1), 4 forced zero bits")
    print(" n      S(n)    sqrt(n)S(n)   E[k_n] (all, merged=0)   E[k_n^2]/n")
    for c in cps:
        S = alive[c] / N
        Ek = ks[c] / N
        Ek2 = k2s[c] / N
        print(f"{c:6d}  {S:.4f}   {math.sqrt(c)*S:7.2f}     {Ek:+.3f}                 {Ek2/c:.4f}")
    print("flip density P(e odd) among survivors by window:")
    for w, (a, b) in wflip.items():
        if b:
            print(f"  steps {w}: {a/b:.4f}  ({b} samples)")
    if ret_f:
        ret_f.sort()
        q = lambda x: ret_f[min(len(ret_f)-1, int(x*len(ret_f)))]
        print(f"|e| at returns of k to 0 (n>=200, survivors): count {len(ret_f)}; quantiles 50% {q(.5)}, 90% {q(.9)}, 99% {q(.99)}; P(|e|<=10) = {sum(1 for x in ret_f if x<=10)/len(ret_f):.3f}")
    mt.sort()
    if mt:
        print(f"merged: {len(mt)}/{N}; median merge time {mt[len(mt)//2]}")
