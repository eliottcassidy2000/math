"""Any-lag Mersenne switch in the Haar model via pair chains driven by one tape.
For odd D <= DMAX: u = q_D, v = p, (k,e) = (-D, 3^-D - 1), i.e. (k, num) = (-D, 1 - 3^D).
All chains are driven by the parities of T^n(p); p = 1 mod 16 gives forced bits 1,0,1,0.
q(T) = P(no chain absorbed by T); q1(T) for D = 1 alone.
usage: python3 pair_chain_anylag.py N NMAX DMAX seed
"""
import random, sys, math
from pair_chain_verify import step

def main(N, NMAX, DMAX, seed):
    rng = random.Random(seed)
    cps = [c for c in [27, 100, 400, 1600, 3200, 6400, 12800, 19000, 25600] if c <= NMAX]
    alive_any = {c: 0 for c in cps}; alive1 = {c: 0 for c in cps}
    first_lag = {}
    Ds = list(range(1, DMAX + 1, 2))
    forced = [1, 0, 1, 0]
    for s in range(N):
        st = {D: (-D, 1 - 3 ** D) for D in Ds}
        any_alive = True; lag1_alive = True
        bits = 0; nb = 0
        ci = 0
        for n in range(NMAX):
            if n < 4:
                beta = forced[n]
            else:
                if nb == 0:
                    bits = rng.getrandbits(64); nb = 64
                beta = bits & 1; bits >>= 1; nb -= 1
            for D in list(st.keys()):
                k, num = st[D]
                k, num, _ = step(k, num, beta)
                if k == 0 and num == 0:
                    del st[D]
                    if any_alive:
                        first_lag[D] = first_lag.get(D, 0) + 1
                    any_alive = False
                    if D == 1:
                        lag1_alive = False
                else:
                    st[D] = (k, num)
            while ci < len(cps) and cps[ci] == n + 1:
                alive_any[cps[ci]] += any_alive; alive1[cps[ci]] += lag1_alive; ci += 1
            if not lag1_alive and not any_alive:
                # still need lag-1 survival only; both dead -> remaining checkpoints contribute 0
                break
            if not any_alive:
                # keep only lag-1 chain to finish q1 bookkeeping
                st = {D: v for D, v in st.items() if D == 1}
    print(f"N={N} NMAX={NMAX} odd D<={DMAX} seed={seed}")
    pts = []
    for c in cps:
        q = alive_any[c] / N; q1 = alive1[c] / N
        pts.append((c, q))
        print(f"  T={c:6d}  q(T)={q:.4f}  q1(T)={q1:.4f}  sqrtT q1={math.sqrt(c)*q1:6.2f}")
    xs = [math.log(c) for c, q in pts if c >= 400 and q > 0]; ys = [math.log(q) for c, q in pts if c >= 400 and q > 0]
    if len(xs) >= 3:
        mx = sum(xs)/len(xs); my = sum(ys)/len(ys)
        sl = sum((x-mx)*(y-my) for x, y in zip(xs, ys)) / sum((x-mx)**2 for x in xs)
        print(f"  least-squares exponent alpha on T>=400: {-sl:.3f}")
    tot = sum(first_lag.values())
    print("  first-switch lag shares:", {D: round(c/tot, 3) for D, c in sorted(first_lag.items())[:6]})

if __name__ == "__main__":
    main(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]))
