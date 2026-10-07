"""Larger Haar Monte Carlo of the lag-1 pair chain with diagnostics.
usage: python3 pair_chain_mc2.py N NMAX seed
"""
import random, sys, math
from pair_chain_verify import step

def main(N, NMAX, seed):
    rng = random.Random(seed)
    cps = [c for c in [27, 100, 400, 1600, 3200, 6400, 12800, 19000, 25600, 51200] if c <= NMAX]
    alive = {c: 0 for c in cps}
    win = [(0, 100), (100, 1000), (1000, 10000), (10000, 10**9)]
    wfl = {w: [0, 0] for w in win}
    epochs = [(200, 1000), (1000, 5000), (5000, 25000), (25000, 10**9)]
    ret = {ep: [] for ep in epochs}
    flips_total = 0; steps_alive = 0
    merge_from = {}
    for s in range(N):
        k, num = 0, 1
        bits = 0; nb = 0; ci = 0
        for n in range(NMAX):
            if n < 4:
                beta = 0
            else:
                if nb == 0:
                    bits = rng.getrandbits(64); nb = 64
                beta = bits & 1; bits >>= 1; nb -= 1
            fl = num & 1
            flips_total += fl; steps_alive += 1
            for w in win:
                if w[0] <= n < w[1]:
                    wfl[w][0] += fl; wfl[w][1] += 1; break
            kp = k; prev = (k, num)
            k, num, _ = step(k, num, beta)
            if k == 0 and num == 0:
                key = (prev[0], prev[1], beta)
                merge_from[key] = merge_from.get(key, 0) + 1
                break
            if k == 0 and kp != 0:
                for ep in epochs:
                    if ep[0] <= n < ep[1]:
                        ret[ep].append(abs(num)); break
            while ci < len(cps) and cps[ci] == n + 1:
                alive[cps[ci]] += 1; ci += 1
    print(f"N={N} NMAX={NMAX} seed={seed}")
    for c in cps:
        S = alive[c] / N; se = math.sqrt(S * (1 - S) / N)
        print(f"  n={c:6d}  S={S:.4f} (se {se:.4f})  sqrt(n) S = {math.sqrt(c)*S:6.2f} (se {math.sqrt(c)*se:.2f})")
    print("  flip density among survivors:", {f"{w[0]}-{w[1] if w[1]<10**9 else 'inf'}": round(a / b, 4) for w, (a, b) in wfl.items() if b})
    print(f"  overall flip density (alive steps) {flips_total/steps_alive:.4f} over {steps_alive} steps")
    for ep, L in ret.items():
        if len(L) > 20:
            L.sort(); q = lambda x: L[min(len(L)-1, int(x * len(L)))]
            tail = {t: sum(1 for x in L if x > t) / len(L) for t in (10, 100, 1000)}
            print(f"  |e| at k-returns to 0, n in {ep}: count {len(L)}, median {q(.5)}, q90 {q(.9)}, q99 {q(.99)}, P(>10,>100,>1000) = {tail[10]:.3f},{tail[100]:.3f},{tail[1000]:.4f}")
    print("  merge entry states (k, e_num, beta):", sorted(merge_from.items(), key=lambda x: -x[1])[:4])

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]))
