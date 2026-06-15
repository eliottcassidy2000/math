#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The "0+0=1" world of H: the achievable odd-H spectrum, its gaps, and a falsifiable
2-adic pattern for the forbidden values.  kind-pasteur-2026-06-14-S5.

H = #directed Hamiltonian paths = I(Omega,2) = sum_k alpha_k 2^k, alpha_0=1 (Redei =
the vacuum unit, the "0+0=1" ground floor; THM-466). Forbidden values (THM-029/075):
H=7 and H=21 never occur. Observation: 7 = (111)_2 = 1+2+4, 21 = (111)_4 = 1+4+16.
=> the 1+2^k+4^k pattern predicts the NEXT forbidden value 73 = (111)_8 = 1+8+64.
TEST it: is 73 (or 273=(111)_16) achievable?  If yes, pattern dead (honest negative);
if no after hard search, a striking 2-adic structure for the forbidden set.

Method: exhaustive achievable-H sets for n<=6 (exact); heavy sampling n=7,8; and a
TARGETED hill-climb trying to realize each special value at growing n.
"""
import sys, itertools, random
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def H_count(out, n):
    """#directed Hamiltonian paths via Held-Karp. out[v] = bitmask of beaten vertices."""
    full = (1<<n)-1
    # dp[mask] = list over last-vertex of path-counts; use dict of arrays
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for mask in range(1<<n):
        row = dp[mask]
        for v in range(n):
            c = row[v]
            if not c: continue
            # extend to w not in mask with v->w (w in out[v])
            avail = out[v] & ~mask
            w = 0; a = avail
            while a:
                if a & 1:
                    dp[mask|(1<<w)][w] += c
                a >>= 1; w += 1
    return sum(dp[full])

def rand_out(n, rng):
    out = [0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.random() < 0.5: out[i] |= 1<<j
            else: out[j] |= 1<<i
    return out

def out_from_bits(bits, pairs, n):
    out=[0]*n
    for k,(i,j) in enumerate(pairs):
        if (bits>>k)&1: out[i]|=1<<j
        else: out[j]|=1<<i
    return out

def exhaustive_H(n):
    pairs=list(itertools.combinations(range(n),2)); m=len(pairs)
    S=set()
    for bits in range(1<<m):
        S.add(H_count(out_from_bits(bits,pairs,n), n))
    return S

def sample_H(n, k, rng):
    S=set()
    for _ in range(k):
        S.add(H_count(rand_out(n,rng), n))
    return S

def gaps(S):
    mx=max(S); odd=[h for h in range(1,mx+1,2)]
    return [h for h in odd if h not in S]

def neighbors_flip(out, n, rng):
    """flip one arc."""
    i,j = rng.randrange(n), rng.randrange(n)
    while i==j: i,j=rng.randrange(n),rng.randrange(n)
    o=[x for x in out]
    if o[i]&(1<<j): o[i]&=~(1<<j); o[j]|=1<<i
    else: o[j]&=~(1<<i); o[i]|=1<<j
    return o

def target_search(target, n, rng, restarts=400, steps=300):
    """hill-climb / random search for a tournament on n vertices with H=target."""
    for _ in range(restarts):
        out=rand_out(n,rng)
        h=H_count(out,n)
        for _ in range(steps):
            if h==target: return out
            o2=neighbors_flip(out,n,rng); h2=H_count(o2,n)
            # accept moves that get closer to target (|h-target| not increasing) + noise
            if abs(h2-target)<=abs(h-target) or rng.random()<0.3:
                out,h=o2,h2
        if h==target: return out
    return None

def main():
    rng=random.Random(614)
    print("=== achievable odd-H spectrum (exhaustive n<=6, sampled n=7,8) ===",flush=True)
    for n in (4,5,6):
        S=exhaustive_H(n); g=gaps(S)
        print(f"   n={n} exhaustive: maxH={max(S)}, odd-gaps in [1,max]={g}",flush=True)
    for n,k in ((7,200000),(8,400000)):
        S=sample_H(n,k,rng); g=gaps(S)
        print(f"   n={n} sampled({k}): maxH={max(S)}, odd-gaps in [1,max]={g if len(g)<40 else str(g[:40])+'...'}",flush=True)

    print("\n=== the 1+2^k+4^k = (111)_{2^k} pattern: {7,21,73,273,1057} ===",flush=True)
    specials=[1+2**k+4**k for k in range(1,6)]
    print(f"   values: {specials}",flush=True)
    for val in specials:
        found_n=None
        for n in range(4,11):
            # only search if val plausibly <= max H at this n (maxH grows fast)
            res=target_search(val,n,rng,restarts=250,steps=250)
            if res is not None:
                found_n=n; break
        print(f"   H={val:5d} (=1+2^{specials.index(val)+1}+4^{specials.index(val)+1}): "
              + (f"ACHIEVABLE at n={found_n}" if found_n else "NOT FOUND (n<=10, hard search) -- candidate forbidden"),flush=True)

if __name__=="__main__":
    main()
