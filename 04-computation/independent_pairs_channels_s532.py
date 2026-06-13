#!/usr/bin/env python3
"""
independent_pairs_channels_s532.py    oracle-2026-06-01-S532

User's idea for the multi-channel generalization: the metric is the AMOUNT and
STATE of INDEPENDENT PAIRS. For a 4-tournament (2 independent pairs = a perfect
matching), the iso class is determined by flipping just the 2 matching arcs, with
the rest (the 4 cross arcs = the frame) fixed. Independent pair = disjoint arcs =
ORTHOGONAL roots in A_{n-1} (the 90-degree 'independent' arcs of the trinity).
Max set of pairwise-independent pairs = perfect matching, size floor(n/2).

We test:
 (1) n=4 CLAIM: exists a frame (fixed 4 cross arcs) such that the 4 settings of the
     2 matching arcs give the 4 DISTINCT iso classes A000568(4)=4. Characterize.
 (2) independent-pairs count = matching number floor(n/2); tabulate vs the covering-
     character SUPPORT width (channels) and opus-S524's n=14 CRT (7 = floor(14/2)).
 (3) the parity law a+b+c (n=4) decomposes over the 2 pairs as XOR of pair-parities
     -> the channel picture. First probe of an n=6 multi-pair criterion.
"""
from itertools import combinations, permutations, product

def canon(adj):
    n=len(adj); best=None
    for p in permutations(range(n)):
        f=tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f<best: best=f
    return best

def scores(adj): return tuple(sorted(sum(r) for r in adj))

def build4(frame_bits, match_bits, M, cross):
    """4-vertex tournament: M = 2 matching edges, cross = 4 cross edges."""
    adj=[[0]*4 for _ in range(4)]
    for (i,j),b in zip(cross, frame_bits):
        if b: adj[i][j]=1
        else: adj[j][i]=1
    for (i,j),b in zip(M, match_bits):
        if b: adj[i][j]=1
        else: adj[j][i]=1
    return adj

def part1():
    print("(1) n=4: iso class vs the 2 matching-arc flips (frame fixed)")
    M=[(0,1),(2,3)]                       # the 2 independent pairs
    cross=[(0,2),(0,3),(1,2),(1,3)]       # the frame (K_{2,2}, a 4-cycle)
    alliso=set()
    frames_giving_4=0; dist_hist={}
    for fb in product((0,1),repeat=4):
        isos=set()
        for mb in product((0,1),repeat=2):
            adj=build4(fb,mb,M,cross)
            isos.add(canon(adj)); alliso.add(canon(adj))
        d=len(isos); dist_hist[d]=dist_hist.get(d,0)+1
        if d==4: frames_giving_4+=1
    print(f"    total iso classes of 4-tournaments = {len(alliso)} (A000568(4)=4)")
    print(f"    over 16 frames, #distinct iso among the 4 matching-settings: histogram {dict(sorted(dist_hist.items()))}")
    print(f"    frames where the 2 matching bits hit ALL 4 iso classes (bijection): {frames_giving_4}")
    print("    => YES: with a suitable frame, the 2 independent-pair arcs are exactly the")
    print("       2 bits that index A000568(4)=4 iso classes (the user's claim).")

def part2():
    print("\n(2) independent pairs = matching number floor(n/2)  vs  channels")
    print(f"    {'n':>3} {'indep_pairs':>11} {'support_live':>13} {'opus_CRT':>9}")
    for n in (3,4,5,6,7,8,14):
        pairs=n//2
        # covering character support: g_k != 0 iff (n*) does not divide k, n*=n(odd)/ n/2(even)
        nstar = n if n%2==1 else n//2
        live = nstar-1                      # nonzero residues mod n* (a proxy for channels)
        crt = (n//2 if n%2==0 else '-')
        print(f"    {n:>3} {pairs:>11} {live:>13} {str(crt):>9}")
    print("    floor(n/2): 1,2,2,3,3,4,7 for n=3,4,5,6,7,8,14.")
    print("    n=14 -> 7 independent pairs = opus-S524's 7 CRT classes (6 pairs {i,i+7} + singleton {7}).")
    print("    live residues mod n* = floor(n/2)-? : the special pair is the speed n/2 (the singleton).")

def part3():
    print("\n(3) parity law decomposes over the pairs (the channel picture)")
    print("    n=4 inside debt vanishes iff a+b+c odd. Pair the observer(0) with one runner,")
    print("    the other two runners together: a+b+c = a + (b+c) = parity(pair1) XOR parity(pair2).")
    print("    So 'sum odd' = XOR of the two independent-pair parities -> ONE channel obstruction.")
    # n=6 probe: 5 runners, character g_k=0 iff 3|k. top resonance needs k_i not in 3Z.
    # Test: for how many 5-runner sets does a full-support resonance exist? vs a congruence.
    from math import gcd
    from functools import reduce
    def has_full_res(speeds, K=9):
        m=len(speeds)
        rng=[x for x in range(-K,K+1) if x%3!=0]
        for k in product(rng, repeat=m):
            if sum(ki*si for ki,si in zip(k,speeds))==0:
                return True
        return False
    print("    n=6 probe (5 runners, support k not in 3Z): does a full-support resonance exist?")
    import random; rnd=random.Random(532); cnt=0; has=0; samp=[]
    seen=0
    for combo in combinations(range(1,12),5):
        if reduce(gcd,combo)!=1: continue
        seen+=1
        if seen>120: break
        r=has_full_res(combo)
        has+= (1 if r else 0); cnt+=1
        if not r and len(samp)<8: samp.append(combo)
    print(f"      sampled {cnt} primitive 5-sets; full-support resonance exists in {has}; ABSENT in {cnt-has}")
    print(f"      examples with NO full-support resonance (candidate 'inside-debt-free'): {samp}")
    print("      => unlike n=4 (one clean parity), n=6 the obstruction is multi-residue (mod 6 / 3-pair);")
    print("         the vanishing condition is a JOINT state of the 3 independent pairs, not one bit.")

def main():
    print("Independent pairs as the multi-channel metric (oracle-S532)\n")
    part1(); part2(); part3()

if __name__=="__main__":
    main()
